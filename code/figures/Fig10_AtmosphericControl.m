%% Repository-relative setup
thisDir = fileparts(mfilename('fullpath'));
codeDir = fileparts(thisDir);
addpath(codeDir,fullfile(codeDir,'functions'));
cfg = project_config();
sys = river_shelf_systems();
close all;

%% ========================================================================
% Figure 10
% Atmospheric modulation of HSI-SIC relationships
%
% Updated for annual_river_shelf_metrics_1979_2024.mat
%
% (a) Raw Pearson correlation vs higher-order partial correlation
%     controlling simultaneously for V10, SAT, and AO
%
% (b) DeltaR = R_partial - R_raw
%     Bootstrap 95% confidence interval
%
% (c) HSI-SIC correlations under positive and negative AO phases
%
% Period: 1979-2024, JAS
%
% IMPORTANT:
% - Uses final HSI_JAS_dt and SIC_JAS_dt directly from annual_river_shelf_metrics_1979_2024.mat
% - Uses ordinary two-sided significance tests (p < 0.05)
% - NO lag-1 effective-sample-size correction
% - AO-phase correlations are recalculated automatically
%% ========================================================================




%% ========================================================================
% 1. LOAD FINAL HSI / SIC DATA
%% ========================================================================

dataDir = cfg.dataDir;

D = load( ...
    fullfile(dataDir,'annual_river_shelf_metrics_1979_2024.mat'));


% ------------------------------------------------------------
% Required variables
% ------------------------------------------------------------

HSI_dt = double(D.HSI_JAS_dt);
SIC_dt = double(D.SIC_JAS_dt);

HSI_z  = double(D.HSI_JAS_z);
SIC_z  = double(D.SIC_JAS_z);

years = double(D.years(:));


% ------------------------------------------------------------
% Convert to 46 × 6
% ------------------------------------------------------------

if size(HSI_dt,1) == 6
    HSI_dt = HSI_dt';
end

if size(SIC_dt,1) == 6
    SIC_dt = SIC_dt';
end

if size(HSI_z,1) == 6
    HSI_z = HSI_z';
end

if size(SIC_z,1) == 6
    SIC_z = SIC_z';
end


assert(isequal(size(HSI_dt),[46 6]), ...
    'HSI_JAS_dt must be 46 × 6.');

assert(isequal(size(SIC_dt),[46 6]), ...
    'SIC_JAS_dt must be 46 × 6.');

assert(isequal(size(HSI_z),[46 6]), ...
    'HSI_JAS_z must be 46 × 6.');

assert(isequal(size(SIC_z),[46 6]), ...
    'SIC_JAS_z must be 46 × 6.');


nYear = length(years);
nSystem = 6;


%% ========================================================================
% 2. LOAD ATMOSPHERIC DATA
%% ========================================================================

windFile = ...
    fullfile(cfg.dataDir,'wind_speed_JAS_1979_2024.mat');

satFile = ...
    fullfile(cfg.dataDir,'surface_air_temperature_JAS_1979_2024.mat');

aoFile = ...
    fullfile(cfg.dataDir,'arctic_oscillation_monthly_1979_2024.mat');


W = load(windFile);
S = load(satFile);
A = load(aoFile);


%% ------------------------------------------------------------------------
% Wind-speed variable compatibility
%
% Prefer V10.
% If the old file still stores the wind-speed matrix as U10,
% use it but rename internally as V10.
%% ------------------------------------------------------------------------

if isfield(W,'V10')

    V10 = double(W.V10);

elseif isfield(W,'U10')

    V10 = double(W.U10);

    fprintf(['NOTE: wind.mat contains U10 rather than V10.\n' ...
             'U10 is used here as the saved 10-m wind-speed series.\n']);

else

    error('Neither V10 nor U10 is found in wind.mat.');

end


%% ------------------------------------------------------------------------
% SAT
%% ------------------------------------------------------------------------

if isfield(S,'SAT')

    SAT = double(S.SAT);

else

    error('SAT variable is not found in SAT.mat.');

end


%% ------------------------------------------------------------------------
% Orientation -> 46 × 6
%% ------------------------------------------------------------------------

if size(V10,1) == 6 && size(V10,2) == 46
    V10 = V10';
end

if size(SAT,1) == 6 && size(SAT,2) == 46
    SAT = SAT';
end


assert(isequal(size(V10),[46 6]), ...
    'V10 must be 46 × 6.');

assert(isequal(size(SAT),[46 6]), ...
    'SAT must be 46 × 6.');


%% ========================================================================
% 3. AO MONTHLY -> JAS ANNUAL MEAN
%% ========================================================================

if ~isfield(A,'AO_index')
    error('AO_index variable is not found in AO_index.mat.');
end

AO_index = double(A.AO_index(:));


% Your original AO file was treated as Jan 1979-Dec 2024
if numel(AO_index) ~= 552

    error(['AO_index should contain exactly 552 monthly values ' ...
           '(Jan 1979-Dec 2024). Current length = %d.'], ...
           numel(AO_index));

end


AO = nan(nYear,1);


for iy = 1:nYear

    idxYear = ...
        (iy-1)*12 + (1:12);

    idxJAS = ...
        idxYear([7 8 9]);

    AO(iy) = ...
        mean(AO_index(idxJAS),'omitnan');

end


%% ========================================================================
% 4. DETREND ATMOSPHERIC VARIABLES
%
% HSI_dt and SIC_dt are already detrended in annual metrics file.
%% ========================================================================

V10_dt = nan(size(V10));
SAT_dt = nan(size(SAT));


for i = 1:nSystem

    V10_dt(:,i) = ...
        local_detrend(V10(:,i),years);

    SAT_dt(:,i) = ...
        local_detrend(SAT(:,i),years);

end


AO_dt = ...
    local_detrend(AO,years);


%% ========================================================================
% 5. BASIN NAMES
%% ========================================================================

basins = { ...
    'Mackenzie–Beaufort', ...
    'Yukon–Chukchi', ...
    'Kolyma–East Siberian', ...
    'Lena–Laptev', ...
    'Yenisei–Kara', ...
    'Ob–Kara'};


shortNames = { ...
    'Mackenzie-Beaufort', ...
    'Yukon-Chukchi', ...
    'Kolyma-East Siberian', ...
    'Lena-Laptev', ...
    'Yenisei-Kara', ...
    'Ob-Kara'};


%% ========================================================================
% 6. RAW AND HIGHER-ORDER PARTIAL CORRELATIONS
%
% Ordinary two-sided significance tests.
% NO effective-sample-size correction.
%% ========================================================================

T = table;

T.Basin = basins';

T.N = nan(nSystem,1);

T.R_raw = nan(nSystem,1);
T.p_raw = nan(nSystem,1);

T.R_partial = nan(nSystem,1);
T.p_partial = nan(nSystem,1);

T.DeltaR = nan(nSystem,1);

T.CI_low = nan(nSystem,1);
T.CI_high = nan(nSystem,1);


% Number of bootstrap resamples
BOOT_N = 3000;


% Reproducibility
rng(20260913);


for i = 1:nSystem

    x = HSI_dt(:,i);
    y = SIC_dt(:,i);

    z1 = V10_dt(:,i);
    z2 = SAT_dt(:,i);
    z3 = AO_dt;


    valid = ...
        isfinite(x) & ...
        isfinite(y) & ...
        isfinite(z1) & ...
        isfinite(z2) & ...
        isfinite(z3);


    x = x(valid);
    y = y(valid);

    Z = [ ...
        z1(valid), ...
        z2(valid), ...
        z3(valid)];


    T.N(i) = length(x);


    %% --------------------------------------------------------
    % Raw Pearson correlation
    %% --------------------------------------------------------

    [T.R_raw(i),T.p_raw(i)] = ...
        local_pearson(x,y);


    %% --------------------------------------------------------
    % Higher-order partial correlation
    %
    % Simultaneously controls:
    % V10, SAT, AO
    %
    % Uses precision-matrix formulation consistent with Eq. (2)
    %% --------------------------------------------------------

    [T.R_partial(i),T.p_partial(i)] = ...
        local_partialcorr_precision(x,y,Z);


    %% --------------------------------------------------------
    % Delta R
    %% --------------------------------------------------------

    T.DeltaR(i) = ...
        T.R_partial(i) - ...
        T.R_raw(i);


    %% --------------------------------------------------------
    % Bootstrap 95% CI for DeltaR
    %% --------------------------------------------------------

    dR_boot = ...
        nan(BOOT_N,1);


    n = length(x);


    for ib = 1:BOOT_N

        idx = ...
            randi(n,n,1);


        xb = ...
            x(idx);

        yb = ...
            y(idx);

        Zb = ...
            Z(idx,:);


        rRaw_b = ...
            local_pearson_r(xb,yb);


        rPar_b = ...
            local_partialcorr_r_precision( ...
            xb,yb,Zb);


        dR_boot(ib) = ...
            rPar_b-rRaw_b;

    end


    dR_boot = ...
        dR_boot(isfinite(dR_boot));


    if length(dR_boot) > 100

        ci = ...
            prctile( ...
            dR_boot, ...
            [2.5 97.5]);


        T.CI_low(i) = ...
            ci(1);

        T.CI_high(i) = ...
            ci(2);

    end

end


%% ========================================================================
% 7. DISPLAY FINAL VALUES
%% ========================================================================

fprintf('\n============================================================\n');
fprintf('FINAL FIGURE 10 RESULTS\n');
fprintf('============================================================\n');

disp(T);


%% ========================================================================
% 8. CHECK RAW CORRELATIONS AGAINST FINAL annual_river_shelf_metrics_1979_2024.mat
%
% These should reproduce the HSI-SIC correlations used in Figure 8.
%% ========================================================================

fprintf('\nRaw HSI-SIC correlations from final annual series:\n');

for i = 1:nSystem

    fprintf('%-24s R = %6.3f, p = %.5f\n', ...
        basins{i}, ...
        T.R_raw(i), ...
        T.p_raw(i));

end


%% ========================================================================
% 9. AO-PHASE HSI-SIC CORRELATIONS
%
% HSI and SIC are already detrended and standardized within each basin.
%
% Positive/negative AO years are selected using the JAS AO index.
% Standardized basin-year anomalies are then pooled across the six
% river-shelf systems.
%% ========================================================================

AO_pos = ...
    AO > 0;

AO_neg = ...
    AO < 0;


X_pos = ...
    HSI_z(AO_pos,:);

Y_pos = ...
    SIC_z(AO_pos,:);


X_neg = ...
    HSI_z(AO_neg,:);

Y_neg = ...
    SIC_z(AO_neg,:);


X_pos = X_pos(:);
Y_pos = Y_pos(:);

X_neg = X_neg(:);
Y_neg = Y_neg(:);


[R_pos,p_pos] = ...
    local_pearson( ...
    X_pos,Y_pos);


[R_neg,p_neg] = ...
    local_pearson( ...
    X_neg,Y_neg);


fprintf('\nAO-phase pooled correlations:\n');

fprintf('AO > 0: R = %.3f, p = %.5f, N = %d\n', ...
    R_pos,p_pos, ...
    sum(isfinite(X_pos) & isfinite(Y_pos)));

fprintf('AO < 0: R = %.3f, p = %.5f, N = %d\n', ...
    R_neg,p_neg, ...
    sum(isfinite(X_neg) & isfinite(Y_neg)));


%% ========================================================================
% 10. SAVE NUMERICAL RESULTS
%% ========================================================================

outDir = cfg.outputDir;


if ~exist(outDir,'dir')
    mkdir(outDir);
end


%% ========================================================================
% 11. FIGURE SETTINGS
%% ========================================================================

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 25 7.6]);


tiledlayout( ...
    1,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');


fontAxis   =5.8;
fontLabel  = 8.2;
fontLegend = 7.0;
fontPanel  = 10;


cRaw = ...
    [0.15 0.45 0.85];

cPartial = ...
    [0.95 0.45 0.10];

cDelta = ...
    [0.20 0.45 0.75];

cAOpos = ...
    [0.55 0.68 0.90];

cAOneg = ...
    [0.10 0.35 0.80];


%% ========================================================================
% 12. PANEL (a)
% Raw vs partial correlations
%% ========================================================================

ax1 = nexttile;

hold(ax1,'on');


barData = ...
    [T.R_raw,T.R_partial];


bh = bar( ...
    ax1, ...
    barData, ...
    'grouped');


bh(1).FaceColor = cRaw;
bh(1).EdgeColor = 'none';

bh(2).FaceColor = cPartial;
bh(2).EdgeColor = 'none';


yline( ...
    ax1, ...
    0, ...
    '-', ...
    'Color',[0.45 0.45 0.45], ...
    'LineWidth',0.8);


xlim(ax1,[0.5 6.5]);
ylim(ax1,[-0.40 1.00]);

yticks(ax1,-0.4:0.2:1.0);


set(ax1, ...
    'XTick',1:6, ...
    'XTickLabel',shortNames, ...
    'FontName','Times New Roman', ...
    'FontSize',fontAxis+0.5, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'Box','off', ...
    'Layer','top');


% basin names bold + small
ax1.XAxis.FontWeight = 'bold';

xtickangle(ax1,40);


ylabel( ...
    ax1, ...
    'Correlation R', ...
    'FontName','Times New Roman', ...
    'FontSize',fontLabel, ...
    'FontWeight','bold');


grid(ax1,'on');

ax1.GridColor = ...
    [0.75 0.75 0.75];

ax1.GridAlpha = ...
    0.30;


lgd = legend( ...
    ax1, ...
    {'R_{raw}','R_{partial}'}, ...
    'Location','southwest', ...
    'Box','off', ...
    'FontName','Times New Roman', ...
    'FontSize',fontLegend, ...
    'Interpreter','tex');


lgd.ItemTokenSize = ...
    [12 6];


%% ------------------------------------------------------------------------
% Significance markers
%% ------------------------------------------------------------------------

drawnow;


xRaw = ...
    bh(1).XEndPoints;

xPar = ...
    bh(2).XEndPoints;


starOffset = ...
    0.04;


for i = 1:nSystem

    %% Raw

    if isfinite(T.p_raw(i)) && ...
       T.p_raw(i) < 0.05

        if T.R_raw(i) >= 0

            yy = ...
                T.R_raw(i)+starOffset;

            va = ...
                'bottom';

        else

            yy = ...
                T.R_raw(i)-starOffset;

            va = ...
                'top';

        end


        text( ...
            ax1, ...
            xRaw(i), ...
            yy, ...
            '*', ...
            'Color',cRaw, ...
            'FontName','Times New Roman', ...
            'FontSize',fontPanel, ...
            'FontWeight','bold', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment',va, ...
            'Clipping','off');

    end


    %% Partial

    if isfinite(T.p_partial(i)) && ...
       T.p_partial(i) < 0.05

        if T.R_partial(i) >= 0

            yy = ...
                T.R_partial(i)+starOffset;

            va = ...
                'bottom';

        else

            yy = ...
                T.R_partial(i)-starOffset;

            va = ...
                'top';

        end


        text( ...
            ax1, ...
            xPar(i), ...
            yy, ...
            '*', ...
            'Color',cPartial, ...
            'FontName','Times New Roman', ...
            'FontSize',fontPanel, ...
            'FontWeight','bold', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment',va, ...
            'Clipping','off');

    end

end


text( ...
    ax1, ...
    0.03,0.98, ...
    '(a)', ...
    'Units','normalized', ...
    'FontName','Times New Roman', ...
    'FontSize',fontPanel, ...
    'FontWeight','bold', ...
    'VerticalAlignment','top');


%% ========================================================================
% 13. PANEL (b)
% Delta R + bootstrap 95% CI
%% ========================================================================

ax2 = nexttile;

hold(ax2,'on');


xpos = ...
    1:nSystem;


dR = ...
    T.DeltaR;


errLow = ...
    dR-T.CI_low;


errHigh = ...
    T.CI_high-dR;


yline( ...
    ax2, ...
    0, ...
    '--', ...
    'Color',[0.45 0.45 0.45], ...
    'LineWidth',0.9);


errorbar( ...
    ax2, ...
    xpos, ...
    dR, ...
    errLow, ...
    errHigh, ...
    'o', ...
    'Color',[0.20 0.20 0.20], ...
    'MarkerFaceColor',cDelta, ...
    'MarkerEdgeColor',[0.10 0.10 0.10], ...
    'MarkerSize',4.8, ...
    'LineWidth',1.0, ...
    'CapSize',6);


xlim(ax2,[0.5 6.5]);


%% automatic y range

ymin = ...
    min(T.CI_low,[],'omitnan');

ymax = ...
    max(T.CI_high,[],'omitnan');


if ~isfinite(ymin) || ...
   ~isfinite(ymax)

    ymin = -0.30;
    ymax = 0.20;

end


ypad = 0.02;


ymin = ...
    floor((ymin-ypad)/0.05)*0.05;

ymax = ...
    ceil((ymax+ypad)/0.05)*0.05;


ylim(ax2,[ymin ymax]);


set(ax2, ...
    'XTick',xpos, ...
    'XTickLabel',shortNames, ...
    'FontName','Times New Roman', ...
    'FontSize',fontAxis+0.5, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'Box','off', ...
    'Layer','top');


ax2.XAxis.FontWeight = ...
    'bold';


xtickangle(ax2,40);


ylabel( ...
    ax2, ...
    '\DeltaR', ...
    'FontName','Times New Roman', ...
    'FontSize',fontLabel, ...
    'FontWeight','bold', ...
    'Interpreter','tex');


grid(ax2,'on');

ax2.GridColor = ...
    [0.75 0.75 0.75];

ax2.GridAlpha = ...
    0.30;


text( ...
    ax2, ...
    0.03,0.98, ...
    '(b)', ...
    'Units','normalized', ...
    'FontName','Times New Roman', ...
    'FontSize',fontPanel, ...
    'FontWeight','bold', ...
    'VerticalAlignment','top');


%% ========================================================================
% 14. PANEL (c)
% Positive vs negative AO phases
%% ========================================================================

ax3 = nexttile;

hold(ax3,'on');


bar( ...
    ax3, ...
    1,R_pos,0.55, ...
    'FaceColor',cAOpos, ...
    'EdgeColor','none');


bar( ...
    ax3, ...
    2,R_neg,0.55, ...
    'FaceColor',cAOneg, ...
    'EdgeColor','none');


yline( ...
    ax3, ...
    0, ...
    '--', ...
    'Color',[0.40 0.40 0.40], ...
    'LineWidth',0.8);


%% ------------------------------------------------------------------------
% Values
%% ------------------------------------------------------------------------

text( ...
    ax3, ...
    1,R_pos+0.025, ...
    sprintf('%.2f',R_pos), ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontName','Times New Roman', ...
    'FontSize',fontAxis, ...
    'FontWeight','bold');


text( ...
    ax3, ...
    2,R_neg+0.025, ...
    sprintf('%.2f',R_neg), ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontName','Times New Roman', ...
    'FontSize',fontAxis, ...
    'FontWeight','bold');


%% ------------------------------------------------------------------------
% Significance
%% ------------------------------------------------------------------------

sPos = ...
    sig_marker(p_pos);

sNeg = ...
    sig_marker(p_neg);


if ~isempty(sPos)

    text( ...
        ax3, ...
        1,R_pos+0.075, ...
        sPos, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontName','Times New Roman', ...
        'FontSize',fontPanel, ...
        'FontWeight','bold');

end


if ~isempty(sNeg)

    text( ...
        ax3, ...
        2,R_neg+0.075, ...
        sNeg, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontName','Times New Roman', ...
        'FontSize',fontPanel, ...
        'FontWeight','bold');

end


set(ax3, ...
    'XTick',[1 2], ...
    'XTickLabel',{'AO > 0','AO < 0'}, ...
    'FontName','Times New Roman', ...
    'FontSize',fontAxis, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'Box','off', ...
    'Layer','top');


ax3.XAxis.FontWeight = ...
    'bold';


ylabel( ...
    ax3, ...
    'Correlation R(HSI, SIC)', ...
    'FontName','Times New Roman', ...
    'FontSize',fontLabel, ...
    'FontWeight','bold');


xlim(ax3,[0.4 2.6]);


% Dynamic y limit
rMax = max([R_pos,R_neg]);

upperY = ...
    max(0.8,ceil((rMax+0.15)*10)/10);

ylim(ax3,[0 upperY]);


grid(ax3,'on');

ax3.GridColor = ...
    [0.75 0.75 0.75];

ax3.GridAlpha = ...
    0.30;


text( ...
    ax3, ...
    0.03,0.98, ...
    '(c)', ...
    'Units','normalized', ...
    'FontName','Times New Roman', ...
    'FontSize',fontPanel, ...
    'FontWeight','bold', ...
    'VerticalAlignment','top');


%% ========================================================================
% 15. EXPORT
%% ========================================================================
drawnow;

set(fig,'InvertHardcopy','off');

exportgraphics( ...
    gcf, ...
    fullfile(cfg.outputDir,'Fig10_AtmosphericControl.png'), ...
    'Resolution',600);


%% ========================================================================
% LOCAL FUNCTIONS
%% ========================================================================


function y = local_detrend(x,years)

    x = double(x(:));
    years = double(years(:));

    y = nan(size(x));


    valid = ...
        isfinite(x) & ...
        isfinite(years);


    if nnz(valid) < 3
        return
    end


    p = ...
        polyfit( ...
        years(valid), ...
        x(valid), ...
        1);


    trend = ...
        polyval( ...
        p, ...
        years(valid));


    y(valid) = ...
        x(valid)-trend;

end


function [r,p] = local_pearson(x,y)
% Ordinary Pearson correlation and ordinary two-sided p-value.

    x = double(x(:));
    y = double(y(:));


    valid = ...
        isfinite(x) & ...
        isfinite(y);


    x = x(valid);
    y = y(valid);


    r = NaN;
    p = NaN;


    if length(x) < 3 || ...
       std(x) == 0 || ...
       std(y) == 0

        return

    end


    [r,p] = ...
        corr( ...
        x,y, ...
        'Rows','complete');

end


function r = local_pearson_r(x,y)

    [r,~] = ...
        local_pearson(x,y);

end


function [r,p] = local_partialcorr_precision(x,y,Z)
% ========================================================================
% Higher-order partial correlation using the precision matrix:
%
% r12.3...k = -P12 / sqrt(P11*P22)
%
% Ordinary two-sided significance test:
% df = n - k - 2
% ========================================================================

    x = double(x(:));
    y = double(y(:));
    Z = double(Z);


    valid = ...
        isfinite(x) & ...
        isfinite(y) & ...
        all(isfinite(Z),2);


    x = x(valid);
    y = y(valid);
    Z = Z(valid,:);


    n = length(x);
    k = size(Z,2);


    r = NaN;
    p = NaN;


    if n <= k+2
        return
    end


    X = ...
        [x y Z];


    R = ...
        corrcoef(X);


    if any(~isfinite(R(:))) || ...
       rcond(R) < 1e-12

        return

    end


    P = ...
        inv(R);


    r = ...
        -P(1,2) / ...
        sqrt(P(1,1)*P(2,2));


    % Numerical protection
    r = ...
        max(-1,min(1,r));


    %% --------------------------------------------------------
    % Ordinary two-sided significance test
    %% --------------------------------------------------------

    df = ...
        n-k-2;


    if df <= 0 || abs(r) >= 1
        return
    end


    tValue = ...
        r * ...
        sqrt(df/(1-r^2));


    p = ...
        2 * ...
        (1-tcdf(abs(tValue),df));

end


function r = local_partialcorr_r_precision(x,y,Z)

    [r,~] = ...
        local_partialcorr_precision( ...
        x,y,Z);

end


function s = sig_marker(p)

    if ~isfinite(p)

        s = '';

    elseif p < 0.001

        s = '***';

    elseif p < 0.01

        s = '**';

    elseif p < 0.05

        s = '*';

    else

        s = '';

    end

end
