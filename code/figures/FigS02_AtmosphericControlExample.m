%% Repository-relative setup
thisDir = fileparts(mfilename('fullpath'));
codeDir = fileparts(thisDir);
addpath(codeDir,fullfile(codeDir,'functions'));
cfg = project_config();
sys = river_shelf_systems();
close all;

%% ========================================================================
% Figure S2
% Atmospheric-control analysis for Mackenzie-Beaufort
%
% (a) Detrended, standardized, and 5-year smoothed
%     HSI, SIC, and AO anomalies
%
% (b) Partial-residual relationship between HSI and SIC
%     after simultaneously controlling for V10, SAT, and AO
%
% Period: 1979-2024, JAS
%
% IMPORTANT:
% This script is consistent with the FINAL Figure 10 calculation.
%% ========================================================================




%% ========================================================================
% 1. PATHS
%% ========================================================================

HSI_file = ...
    fullfile(cfg.dataDir,'annual_river_shelf_metrics_1979_2024.mat');

wind_file = ...
    fullfile(cfg.dataDir,'wind_speed_JAS_1979_2024.mat');

SAT_file = ...
    fullfile(cfg.dataDir,'surface_air_temperature_JAS_1979_2024.mat');

AO_file = ...
    fullfile(cfg.dataDir,'arctic_oscillation_monthly_1979_2024.mat');

outDir = cfg.outputDir;

if ~exist(outDir,'dir')
    mkdir(outDir);
end


%% ========================================================================
% 2. LOAD FINAL HSI AND SIC
%% ========================================================================

D = load(HSI_file);

if ~isfield(D,'HSI_JAS_dt')
    error('Annual metrics file does not contain HSI_JAS_dt.');
end

if ~isfield(D,'SIC_JAS_dt')
    error('Annual metrics file does not contain SIC_JAS_dt.');
end

if ~isfield(D,'years')
    error('Annual metrics file does not contain years.');
end


HSI_dt = double(D.HSI_JAS_dt);
SIC_dt = double(D.SIC_JAS_dt);
years  = double(D.years(:));


% Convert to 46 × 6
if size(HSI_dt,1) == 6
    HSI_dt = HSI_dt';
end

if size(SIC_dt,1) == 6
    SIC_dt = SIC_dt';
end


assert(isequal(size(HSI_dt),[46 6]), ...
    'HSI_JAS_dt must be 46 × 6.');

assert(isequal(size(SIC_dt),[46 6]), ...
    'SIC_JAS_dt must be 46 × 6.');


nYear = length(years);


%% ========================================================================
% 3. LOAD V10 AND SAT
%% ========================================================================

W = load(wind_file);
S = load(SAT_file);


% ------------------------------------------------------------------------
% Wind-speed compatibility
% ------------------------------------------------------------------------

if isfield(W,'V10')

    V10 = double(W.V10);

elseif isfield(W,'U10')

    % Previous files used U10 as the saved wind-speed matrix
    V10 = double(W.U10);

    fprintf('NOTE: wind.mat contains U10; using it as V10 series.\n');

else

    error('Neither V10 nor U10 exists in wind.mat.');

end


if ~isfield(S,'SAT')
    error('SAT variable not found in SAT.mat.');
end

SAT = double(S.SAT);


if size(V10,1) == 6
    V10 = V10';
end

if size(SAT,1) == 6
    SAT = SAT';
end


assert(isequal(size(V10),[46 6]), ...
    'V10 must be 46 × 6.');

assert(isequal(size(SAT),[46 6]), ...
    'SAT must be 46 × 6.');


%% ========================================================================
% 4. AO MONTHLY -> JAS ANNUAL MEAN
%% ========================================================================

A = load(AO_file);

if ~isfield(A,'AO_index')
    error('AO_index not found in AO_index.mat.');
end

AO_index = double(A.AO_index(:));


if numel(AO_index) ~= 552
    error('AO_index should contain 552 monthly values for 1979-2024.');
end


AO = nan(nYear,1);


for iy = 1:nYear

    idxJAS = ...
        (iy-1)*12 + [7 8 9];

    AO(iy) = ...
        mean(AO_index(idxJAS),'omitnan');

end


%% ========================================================================
% 5. DETREND ATMOSPHERIC VARIABLES
%
% HSI_dt and SIC_dt are already detrended in annual metrics file.
%% ========================================================================

V10_dt = nan(size(V10));
SAT_dt = nan(size(SAT));


for ib = 1:6

    V10_dt(:,ib) = ...
        local_detrend(V10(:,ib),years);

    SAT_dt(:,ib) = ...
        local_detrend(SAT(:,ib),years);

end


AO_dt = ...
    local_detrend(AO,years);


%% ========================================================================
% 6. MACKENZIE-BEAUFORT SERIES
%% ========================================================================

ib = 1;      % Mackenzie-Beaufort


x = HSI_dt(:,ib);
y = SIC_dt(:,ib);

v = V10_dt(:,ib);
s = SAT_dt(:,ib);
a = AO_dt;


valid = ...
    isfinite(x) & ...
    isfinite(y) & ...
    isfinite(v) & ...
    isfinite(s) & ...
    isfinite(a);


x = x(valid);
y = y(valid);

v = v(valid);
s = s(valid);
a = a(valid);

yearsUse = years(valid);


Z = [v s a];


%% ========================================================================
% 7. HIGHER-ORDER PARTIAL CORRELATION
%
% Same precision-matrix definition used for Figure 10.
%% ========================================================================

[rPartial,pPartial] = ...
    local_partialcorr_precision(x,y,Z);


fprintf('\n============================================\n');
fprintf('Mackenzie-Beaufort partial correlation\n');
fprintf('============================================\n');
fprintf('R_partial = %.5f\n',rPartial);
fprintf('p_partial = %.6g\n',pPartial);


%% ========================================================================
% 8. PARTIAL RESIDUALS
%
% Regress both HSI and SIC separately against:
% V10 + SAT + AO
%
% Correlation between residuals = higher-order partial correlation.
%% ========================================================================

Xdesign = ...
    [ones(length(x),1) Z];


betaX = ...
    Xdesign \ x;

betaY = ...
    Xdesign \ y;


resHSI = ...
    x - Xdesign*betaX;

resSIC = ...
    y - Xdesign*betaY;


%% ------------------------------------------------------------------------
% Standardize residuals only for visual presentation
%
% This does NOT change their correlation.
%% ------------------------------------------------------------------------

resHSI_z = ...
    local_zscore(resHSI);

resSIC_z = ...
    local_zscore(resSIC);


rResidual = ...
    corr( ...
    resHSI_z, ...
    resSIC_z, ...
    'Rows','complete');


fprintf('Correlation of partial residuals = %.5f\n',rResidual);


if abs(rResidual-rPartial) > 1e-10

    warning('Residual correlation does not equal partial correlation.');

end


%% ========================================================================
% 9. PANEL (a): STANDARDIZED + 5-YEAR SMOOTHED SERIES
%% ========================================================================

% ------------------------------------------------------------------------
% Standardize detrended HSI / SIC / AO
%% ------------------------------------------------------------------------

HSI_z = ...
    local_zscore(HSI_dt(:,ib));

SIC_z = ...
    local_zscore(SIC_dt(:,ib));

AO_z = ...
    local_zscore(AO_dt);


%% ------------------------------------------------------------------------
% 5-year moving mean
%% ------------------------------------------------------------------------

smoothWindow = 5;


HSI_sm = ...
    movmean( ...
    HSI_z, ...
    smoothWindow, ...
    'omitnan', ...
    'Endpoints','shrink');


SIC_sm = ...
    movmean( ...
    SIC_z, ...
    smoothWindow, ...
    'omitnan', ...
    'Endpoints','shrink');


AO_sm = ...
    movmean( ...
    AO_z, ...
    smoothWindow, ...
    'omitnan', ...
    'Endpoints','shrink');


%% ========================================================================
% 10. OPTIONAL: RECOMPUTE TABLE S1
%
% Pooled atmospheric controls across all six systems.
%
% AO is repeated for each river-shelf system because the same AO series
% enters each basin-specific partial-correlation analysis.
%% ========================================================================

V_pool = ...
    V10_dt(:);

SAT_pool = ...
    SAT_dt(:);

AO_pool = ...
    repmat(AO_dt,6,1);


validPool = ...
    isfinite(V_pool) & ...
    isfinite(SAT_pool) & ...
    isfinite(AO_pool);


Controls = [ ...
    V_pool(validPool), ...
    SAT_pool(validPool), ...
    AO_pool(validPool)];


R_controls = ...
    corr( ...
    Controls, ...
    'Rows','complete');


fprintf('\n============================================\n');
fprintf('Table S1: pooled atmospheric correlations\n');
fprintf('============================================\n');

disp(array2table( ...
    R_controls, ...
    'VariableNames',{'V10','SAT','AO'}, ...
    'RowNames',{'V10','SAT','AO'}));


%% ========================================================================
% 11. FIGURE
%% ========================================================================

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 16.5 7.2]);


tl = tiledlayout( ...
    1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


fontAxis  = 7.2;
fontLabel = 8.0;
fontPanel = 10;
fontLegend = 7.0;


%% ------------------------------------------------------------------------
% Colors
%% ------------------------------------------------------------------------

cHSI = [0.05 0.42 0.88];
cSIC = [0.95 0.35 0.08];
cAO  = [0.35 0.35 0.35];

cScatter = [0.25 0.48 0.88];


%% ========================================================================
% 12. PANEL (a)
%% ========================================================================

ax1 = nexttile(tl,1);

hold(ax1,'on');


h1 = plot( ...
    ax1, ...
    years, ...
    HSI_sm, ...
    '-o', ...
    'Color',cHSI, ...
    'MarkerFaceColor','w', ...
    'MarkerEdgeColor',cHSI, ...
    'MarkerSize',2.8, ...
    'LineWidth',1.05);


h2 = plot( ...
    ax1, ...
    years, ...
    SIC_sm, ...
    '-s', ...
    'Color',cSIC, ...
    'MarkerFaceColor','w', ...
    'MarkerEdgeColor',cSIC, ...
    'MarkerSize',2.5, ...
    'LineWidth',1.00);


h3 = plot( ...
    ax1, ...
    years, ...
    AO_sm, ...
    '-', ...
    'Color',cAO, ...
    'LineWidth',0.95);


yline( ...
    ax1, ...
    0, ...
    '-', ...
    'Color',[0.75 0.75 0.75], ...
    'LineWidth',0.6);


xlim(ax1,[1979 2024]);

xticks(ax1,1980:10:2020);


% Use symmetric y range
allSeries = [ ...
    HSI_sm(:); ...
    SIC_sm(:); ...
    AO_sm(:)];


ymax = ...
    max(abs(allSeries),[],'omitnan');

ymax = ...
    ceil((ymax+0.05)*2)/2;


ylim(ax1,[-ymax ymax]);


xlabel( ...
    ax1, ...
    'Year', ...
    'FontName','Times New Roman', ...
    'FontSize',fontLabel, ...
    'FontWeight','bold');


ylabel( ...
    ax1, ...
    'Standardized anomaly', ...
    'FontName','Times New Roman', ...
    'FontSize',fontLabel, ...
    'FontWeight','bold');


grid(ax1,'on');

ax1.GridColor = ...
    [0.80 0.80 0.80];

ax1.GridAlpha = ...
    0.25;


set( ...
    ax1, ...
    'FontName','Times New Roman', ...
    'FontSize',fontAxis, ...
    'FontWeight','bold', ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'Box','off', ...
    'Layer','top');


h=legend( ...
    ax1, ...
    [h1 h2 h3], ...
    {'HSI','SIC','AO'}, ...
    'Location','southwest', ...
    'Orientation','horizontal', ...
    'Box','off', ...
    'FontName','Times New Roman', ...
    'FontSize',fontLegend, ...
    'FontWeight','bold');
h.ItemTokenSize=[7,2];

text( ...
    ax1, ...
    0.03, ...
    0.95, ...
    '(a)', ...
    'Units','normalized', ...
    'FontName','Times New Roman', ...
    'FontSize',fontPanel, ...
    'FontWeight','bold', ...
    'VerticalAlignment','top');


title( ...
    ax1, ...
    'Mackenzie-Beaufort', ...
    'FontName','Times New Roman', ...
    'FontSize',8.5, ...
    'FontWeight','bold');


%% ========================================================================
% 13. PANEL (b): PARTIAL RESIDUALS
%% ========================================================================

ax2 = nexttile(tl,2);

hold(ax2,'on');


%% ------------------------------------------------------------------------
% OLS line on standardized residuals
%% ------------------------------------------------------------------------

pp = ...
    polyfit( ...
    resHSI_z, ...
    resSIC_z, ...
    1);


xmin = ...
    min(resHSI_z);

xmax = ...
    max(resHSI_z);


xpad = ...
    0.08*(xmax-xmin);


xx = ...
    linspace( ...
    xmin-xpad, ...
    xmax+xpad, ...
    100);


yy = ...
    polyval(pp,xx);


%% ------------------------------------------------------------------------
% Scatter
%% ------------------------------------------------------------------------

scatter( ...
    ax2, ...
    resHSI_z, ...
    resSIC_z, ...
    18, ...
    'MarkerFaceColor',cScatter, ...
    'MarkerEdgeColor','none');


%% ------------------------------------------------------------------------
% Regression line
%% ------------------------------------------------------------------------

plot( ...
    ax2, ...
    xx, ...
    yy, ...
    'k-', ...
    'LineWidth',1.15);


%% ------------------------------------------------------------------------
% Zero reference lines
%% ------------------------------------------------------------------------

xline( ...
    ax2, ...
    0, ...
    '-', ...
    'Color',[0.82 0.82 0.82], ...
    'LineWidth',0.5);


yline( ...
    ax2, ...
    0, ...
    '-', ...
    'Color',[0.82 0.82 0.82], ...
    'LineWidth',0.5);


%% ------------------------------------------------------------------------
% Axes
%% ------------------------------------------------------------------------

xlabel( ...
    ax2, ...
    'HSI residual', ...
    'FontName','Times New Roman', ...
    'FontSize',fontLabel, ...
    'FontWeight','bold');


ylabel( ...
    ax2, ...
    'SIC residual', ...
    'FontName','Times New Roman', ...
    'FontSize',fontLabel, ...
    'FontWeight','bold');


grid(ax2,'on');

ax2.GridColor = ...
    [0.80 0.80 0.80];

ax2.GridAlpha = ...
    0.25;


set( ...
    ax2, ...
    'FontName','Times New Roman', ...
    'FontSize',fontAxis, ...
    'FontWeight','bold', ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'Box','off', ...
    'Layer','top');


%% ------------------------------------------------------------------------
% Correlation annotation
%% ------------------------------------------------------------------------

starText = ...
    sig_marker(pPartial);


if pPartial < 0.001

    statText = ...
        sprintf( ...
        'r = %.2f ***', ...
        rPartial);

elseif pPartial < 0.01

    statText = ...
        sprintf( ...
        'r = %.2f **', ...
        rPartial);

elseif pPartial < 0.05

    statText = ...
        sprintf( ...
        'r = %.2f *', ...
        rPartial);

else

    statText = ...
        sprintf( ...
        'r = %.2f', ...
        rPartial);

end


text( ...
    ax2, ...
    0.94, ...
    0.08, ...
    statText, ...
    'Units','normalized', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','bottom', ...
    'FontName','Times New Roman', ...
    'FontSize',8.0, ...
    'FontWeight','bold');


text( ...
    ax2, ...
    0.03, ...
    0.95, ...
    '(b)', ...
    'Units','normalized', ...
    'FontName','Times New Roman', ...
    'FontSize',fontPanel, ...
    'FontWeight','bold', ...
    'VerticalAlignment','top');


title( ...
    ax2, ...
    'Partial residuals', ...
    'FontName','Times New Roman', ...
    'FontSize',8.5, ...
    'FontWeight','bold');


%% ========================================================================
% 14. EXPORT
%% ========================================================================

drawnow;

set(fig,'InvertHardcopy','off');

exportgraphics( ...
    gcf, ...
    fullfile(cfg.outputDir,'FigS02_AtmosphericControlExample.png'), ...
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


function z = local_zscore(x)

    x = double(x(:));

    z = nan(size(x));


    valid = ...
        isfinite(x);


    if nnz(valid) < 2
        return
    end


    mu = ...
        mean(x(valid));


    sigma = ...
        std(x(valid));


    if sigma == 0
        return
    end


    z(valid) = ...
        (x(valid)-mu)/sigma;

end


function [r,p] = local_partialcorr_precision(x,y,Z)

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


    r = ...
        max(-1,min(1,r));


    df = ...
        n-k-2;


    if df <= 0 || ...
       abs(r) >= 1

        return
    end


    tValue = ...
        r * ...
        sqrt(df/(1-r^2));


    p = ...
        2 * ...
        (1-tcdf(abs(tValue),df));

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
