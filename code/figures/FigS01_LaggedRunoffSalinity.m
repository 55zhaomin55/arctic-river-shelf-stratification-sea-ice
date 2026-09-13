%% Repository-relative setup
thisDir = fileparts(mfilename('fullpath'));
codeDir = fileparts(thisDir);
addpath(codeDir,fullfile(codeDir,'functions'));
cfg = project_config();
sys = river_shelf_systems();
close all;

discharge = load_monthly_discharge(cfg.dataDir);


%% ========================================================================
%% 1. PATH
%% ========================================================================

data_dir = ...
    cfg.dataDir;

cache_file = fullfile( ...
    data_dir, ...
    'runoff_salinity_distance_bins_1979_2024.mat');

out_png = fullfile( ...
    data_dir, ...
    'FigS01_LaggedRunoffSalinity.png');


%% ========================================================================
%% 2. CHECK INPUT
%% ========================================================================

if ~exist('discharge','var')
    error(['Variable "discharge" is not available in workspace. ', ...
           'Please load the 6 x 636 discharge matrix first.']);
end

if size(discharge,1) ~= 6
    error('discharge must contain six rows.');
end

if ~exist(cache_file,'file')
    error('Processed salinity file not found:\n%s',cache_file);
end


%% ========================================================================
%% 3. LOAD PROCESSED SALINITY
%% ========================================================================

C = load(cache_file);

SAL_DIST      = C.SAL_DIST;
N_VALID_CELLS = C.N_VALID_CELLS;
TIME_COMMON   = C.TIME_COMMON;
dist_edges    = C.dist_edges;
dist_centers  = C.dist_centers;

fprintf('\nProcessed salinity loaded.\n');

fprintf('Period: %s to %s\n', ...
    datestr(TIME_COMMON(1),'yyyy-mm'), ...
    datestr(TIME_COMMON(end),'yyyy-mm'));


%% ========================================================================
%% 4. BASIC INFORMATION
%% ========================================================================

basin_names = {
    'Mackenzie–Beaufort'
    'Yukon–Chukchi'
    'Kolyma–East Siberian'
    'Lena–Laptev'
    'Yenisei–Kara'
    'Ob–Kara'
    };

panel_letters = {
    '(a)','(b)','(c)', ...
    '(d)','(e)','(f)'
    };

radii = [150 300 450 600];

lags = 0:3;

nBasin  = 6;
nRadius = length(radii);
nLag    = length(lags);

year_start = 1979;
year_end   = 2024;

MIN_N = 100;

MIN_TOTAL_CELLS = 20;

ALPHA = 0.05;


%% ========================================================================
%% 5. COLORBAR SETTINGS
%% ========================================================================

% Current correlations are mostly within about +/-0.30.
% A narrower range gives much better contrast than [-1,1].

COLOR_LIMIT = 0.30;

COLOR_TICKS = ...
    -0.30:0.10:0.30;


%% ========================================================================
%% 6. DISCHARGE TIME AXIS
%% ========================================================================

nQ = size(discharge,2);

q_dates = ...
    datetime(1972,1,15) + ...
    calmonths(0:nQ-1);

q_dates = q_dates(:);


%% ========================================================================
%% 7. ALIGN DISCHARGE WITH SALINITY PERIOD
%% ========================================================================

nTime = length(TIME_COMMON);

Q_COMMON = ...
    nan(nBasin,nTime);

for t = 1:nTime

    yy = year(TIME_COMMON(t));
    mm = month(TIME_COMMON(t));

    iq = ...
        year(q_dates)==yy & ...
        month(q_dates)==mm;

    if sum(iq) ~= 1
        continue
    end

    Q_COMMON(:,t) = ...
        discharge(:,iq);

end


%% ========================================================================
%% 8. SALINITY WITHIN EACH SAMPLING RADIUS
%% ========================================================================

S_RADIUS = ...
    nan(nBasin,nTime,nRadius);

N_RADIUS = ...
    nan(nBasin,nTime,nRadius);


for b = 1:nBasin

    for ir = 1:nRadius

        R0 = radii(ir);

        % Include all 50-km annuli located inside radius R0
        use_bins = ...
            dist_edges(2:end) <= R0 + 1e-10;

        Sring = ...
            squeeze(SAL_DIST(b,:,use_bins));

        Wring = ...
            squeeze(N_VALID_CELLS(b,:,use_bins));

        if size(Sring,1) ~= nTime
            Sring = Sring';
        end

        if size(Wring,1) ~= nTime
            Wring = Wring';
        end

        valid = ...
            isfinite(Sring) & ...
            isfinite(Wring) & ...
            Wring > 0;

        W = Wring;
        W(~valid) = 0;

        S0 = Sring;
        S0(~valid) = 0;

        numerator = ...
            sum(S0 .* W,2);

        denominator = ...
            sum(W,2);

        Smean = ...
            numerator ./ denominator;

        Smean(denominator < MIN_TOTAL_CELLS) = NaN;

        S_RADIUS(b,:,ir) = ...
            Smean;

        N_RADIUS(b,:,ir) = ...
            denominator;

    end

end


%% ========================================================================
%% 9. MONTHLY ANOMALIES + DETRENDING
%% ========================================================================

Q_ANOM = ...
    nan(size(Q_COMMON));

S_ANOM = ...
    nan(size(S_RADIUS));


for b = 1:nBasin

    q = ...
        squeeze(Q_COMMON(b,:))';

    Q_ANOM(b,:) = ...
        monthly_anomaly_detrend( ...
        q, ...
        TIME_COMMON);

    for ir = 1:nRadius

        s = ...
            squeeze(S_RADIUS(b,:,ir))';

        S_ANOM(b,:,ir) = ...
            monthly_anomaly_detrend( ...
            s, ...
            TIME_COMMON);

    end

end


%% ========================================================================
%% 10. LAGGED CORRELATIONS
%% ========================================================================

R_LAG = ...
    nan(nBasin,nRadius,nLag);

P_LAG = ...
    nan(nBasin,nRadius,nLag);

NEFF = ...
    nan(nBasin,nRadius,nLag);

NPAIR = ...
    nan(nBasin,nRadius,nLag);


for b = 1:nBasin

    q = ...
        squeeze(Q_ANOM(b,:))';

    for ir = 1:nRadius

        s = ...
            squeeze(S_ANOM(b,:,ir))';

        for ilag = 1:nLag

            L = lags(ilag);

            if L == 0

                x = q;
                y = s;

            else

                x = q(1:end-L);
                y = s(1+L:end);

            end

            good = ...
                isfinite(x) & ...
                isfinite(y);

            x = x(good);
            y = y(good);

            n = length(x);

            NPAIR(b,ir,ilag) = n;

            if n < MIN_N
                continue
            end

            %% Pearson correlation

            CC = ...
                corrcoef(x,y);

            r = ...
                CC(1,2);

            R_LAG(b,ir,ilag) = r;
            %% Two-sided Pearson significance
            [~,Ptmp] = corrcoef(x,y);
            P_LAG(b,ir,ilag) = Ptmp(1,2);

        end

    end

end


%% ========================================================================
%% 11. CREATE FIGURE: 3 ROWS x 2 COLUMNS
%% ========================================================================

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 19.5 22.0]);


tl = tiledlayout( ...
    fig, ...
    3,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


colormap(fig,blue_white_red_soft(256));


for b = 1:nBasin

    ax = ...
        nexttile(tl,b);

    Z = ...
        squeeze(R_LAG(b,:,:));

    P = ...
        squeeze(P_LAG(b,:,:));


    %% --------------------------------------------------------------------
    % Heatmap
    % ---------------------------------------------------------------------

    himg = imagesc( ...
        ax, ...
        1:nLag, ...
        1:nRadius, ...
        Z);

    set(ax,'YDir','normal');

    set(himg, ...
        'AlphaData', ...
        double(isfinite(Z)));

    ax.Color = [1 1 1];

    clim(ax, ...
        [-COLOR_LIMIT COLOR_LIMIT]);

    hold(ax,'on');


    %% --------------------------------------------------------------------
    % White cell boundaries
    % ---------------------------------------------------------------------

    for xx = 1.5:1:(nLag-0.5)

        xline(ax,xx, ...
            '-', ...
            'Color',[1 1 1], ...
            'LineWidth',0.8);

    end

    for yy = 1.5:1:(nRadius-0.5)

        yline(ax,yy, ...
            '-', ...
            'Color',[1 1 1], ...
            'LineWidth',0.8);

    end


    %% --------------------------------------------------------------------
    % Correlation numbers
    % ---------------------------------------------------------------------

    for ir = 1:nRadius

        for ilag = 1:nLag

            r = Z(ir,ilag);
            p = P(ir,ilag);

            if ~isfinite(r)
                continue
            end

            if p < ALPHA

                txt = ...
                    sprintf('%.2f*',r);

            else

                txt = ...
                    sprintf('%.2f',r);

            end


            % White text only in strongly saturated cells
            if abs(r) >= 0.22

                txtColor = [1 1 1];

            else

                txtColor = [0.10 0.10 0.10];

            end


            text( ...
                ax, ...
                ilag, ...
                ir, ...
                txt, ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontName','Times New Roman', ...
                'FontSize',9.4, ...
                'Color',txtColor);

        end

    end


    %% --------------------------------------------------------------------
    % Axes
    % ---------------------------------------------------------------------

    xlim(ax,[0.5 nLag+0.5]);

    ylim(ax,[0.5 nRadius+0.5]);


    xticks(ax,1:nLag);

    xticklabels(ax, ...
        {'Lag 0','Lag 1','Lag 2','Lag 3'});


    yticks(ax,1:nRadius);

    yticklabels(ax, ...
        {'150','300','450','600'});


    %% --------------------------------------------------------------------
    % X label only on bottom row: panels e and f
    % ---------------------------------------------------------------------

    if b >= 5

        xlabel(ax, ...
            'Runoff lag (month)', ...
            'FontName','Times New Roman', ...
            'FontSize',10.2);

    else

        xlabel(ax,'');

    end


    %% --------------------------------------------------------------------
    % Y label only on left column
    % ---------------------------------------------------------------------

    if mod(b,2) == 1

        ylabel(ax, ...
            'Sampling radius (km)', ...
            'FontName','Times New Roman', ...
            'FontSize',10.2);

    else

        ylabel(ax,'');

    end


    %% --------------------------------------------------------------------
    % Panel title
    % ---------------------------------------------------------------------

    title(ax, ...
        [panel_letters{b} ' ' basin_names{b}], ...
        'FontName','Times New Roman', ...
        'FontSize',12.2, ...
        'FontWeight','Bold');


    %% --------------------------------------------------------------------
    % Style
    % ---------------------------------------------------------------------

    set(ax, ...
        'FontName','Times New Roman', ...
        'FontSize',9.5, ...
        'LineWidth',0.75, ...
        'TickDir','out', ...
        'TickLength',[0.012 0.012], ...
        'Box','on', ...
        'Layer','top');

    ax.XColor = ...
        [0.15 0.15 0.15];

    ax.YColor = ...
        [0.15 0.15 0.15];

end


%% ========================================================================
%% 12. SHARED COLORBAR
%% ========================================================================

cb = colorbar;

cb.Layout.Tile = 'east';


cb.FontName = ...
    'Times New Roman';

cb.FontSize = ...
    9.5;

cb.LineWidth = ...
    0.7;

cb.TickDirection = ...
    'out';


cb.Ticks = ...
    COLOR_TICKS;


ylabel(cb, ...
    'Correlation coefficient, R', ...
    'FontName','Times New Roman', ...
    'FontSize',10.2);


%% ========================================================================
%% 13. EXPORT PNG ONLY
%% ========================================================================

set(gcf,'InvertHardcopy','off');
exportgraphics(gcf,out_png,'Resolution',600);



fprintf('\n============================================================\n');
fprintf('Figure completed successfully.\n');
fprintf('============================================================\n');

fprintf('Saved PNG:\n%s\n',out_png);



%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================


function y = monthly_anomaly_detrend(x,dates)

x = double(x(:));
dates = dates(:);

y = nan(size(x));


%% Remove monthly climatology

for m = 1:12

    ii = ...
        month(dates)==m & ...
        isfinite(x);

    if ~any(ii)
        continue
    end

    clim = ...
        mean(x(ii),'omitnan');

    y(ii) = ...
        x(ii)-clim;

end


%% Remove linear trend

good = ...
    isfinite(y);

if sum(good) < 3
    return
end

t = ...
    (1:length(y))';

p = ...
    polyfit( ...
    t(good), ...
    y(good), ...
    1);

trend = ...
    polyval(p,t);

y(good) = ...
    y(good)-trend(good);

end



function r = lag1_corr(x)

x = double(x(:));

if length(x)<3

    r = NaN;
    return

end

x1 = x(1:end-1);
x2 = x(2:end);

good = ...
    isfinite(x1) & ...
    isfinite(x2);

if sum(good)<3

    r = NaN;
    return

end

C = ...
    corrcoef( ...
    x1(good), ...
    x2(good));

r = ...
    C(1,2);

end



function cmap = blue_white_red_soft(n)

if nargin < 1
    n = 256;
end


blue_dark = ...
    [0.08 0.28 0.68];

blue_light = ...
    [0.70 0.82 0.94];

white = ...
    [1.00 1.00 1.00];

red_light = ...
    [0.96 0.76 0.72];

red_dark = ...
    [0.72 0.06 0.10];


n_half = ...
    floor(n/2);

n1 = ...
    round(n_half*0.55);

n2 = ...
    n_half-n1;

n3 = ...
    n2;

n4 = ...
    n-n1-n2-n3;


c1 = ...
    interp_colors( ...
    blue_dark, ...
    blue_light, ...
    n1);

c2 = ...
    interp_colors( ...
    blue_light, ...
    white, ...
    n2);

c3 = ...
    interp_colors( ...
    white, ...
    red_light, ...
    n3);

c4 = ...
    interp_colors( ...
    red_light, ...
    red_dark, ...
    n4);


cmap = [
    c1
    c2
    c3
    c4
    ];

end



function C = interp_colors(c1,c2,n)

if n <= 1

    C = c1;
    return

end


C = ...
    zeros(n,3);


for k = 1:3

    C(:,k) = ...
        linspace( ...
        c1(k), ...
        c2(k), ...
        n)';

end

end
