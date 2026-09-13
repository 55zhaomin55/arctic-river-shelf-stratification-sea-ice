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


%% ========================================================================
%% 2. CHECK INPUT
%% ========================================================================

if ~exist('discharge','var')
    error(['Variable "discharge" is not available in the workspace. ', ...
           'Please load the 6 x 636 monthly discharge matrix first.']);
end

if size(discharge,1) ~= 6
    error('discharge must contain six rows.');
end

if ~exist(cache_file,'file')
    error('Processed salinity file not found:\n%s',cache_file);
end


%% ========================================================================
%% 3. LOAD PROCESSED ORAS5 SALINITY
%% ========================================================================

C = load(cache_file);

SAL_DIST      = C.SAL_DIST;
N_VALID_CELLS = C.N_VALID_CELLS;
TIME_COMMON   = C.TIME_COMMON;
dist_edges    = C.dist_edges;
dist_centers  = C.dist_centers;

fprintf('\nProcessed ORAS5 salinity loaded.\n');

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

months_plot = 5:10;

month_labels = {
    'May','Jun','Jul','Aug','Sep','Oct'
    };

nBasin = 6;
nMonth = numel(months_plot);
nDist  = numel(dist_centers);

year_start = 1979;
year_end   = 2024;


%% ========================================================================
%% 5. ANALYSIS SETTINGS
%% ========================================================================

MIN_VALID_YEARS = 25;
MIN_VALID_OCEAN_CELLS = 5;

DO_DETREND = true;

% [] = automatically determine robust symmetric color range
COLOR_LIMIT = [];


%% ========================================================================
%% 6. DISCHARGE TIME AXIS
%% ========================================================================

nQ = size(discharge,2);

q_dates = ...
    datetime(1972,1,15) + ...
    calmonths(0:nQ-1);

q_dates = q_dates(:);


%% ========================================================================
%% 7. INITIALIZE RESULTS
%% ========================================================================

BETA  = nan(nBasin,nMonth,nDist);
RVAL  = nan(nBasin,nMonth,nDist);
PVAL  = nan(nBasin,nMonth,nDist);
NPAIR = nan(nBasin,nMonth,nDist);


%% ========================================================================
%% 8. SAME-MONTH RUNOFF -> SAME-MONTH SALINITY
%% ========================================================================

for b = 1:nBasin

    fprintf('\n===============================================\n');
    fprintf('%s\n',basin_names{b});
    fprintf('===============================================\n');

    sal_all  = squeeze(SAL_DIST(b,:,:));
    cell_all = squeeze(N_VALID_CELLS(b,:,:));

    for im = 1:nMonth

        mon = months_plot(im);

        %% Same-calendar-month salinity
        idx_sal = ...
            month(TIME_COMMON) == mon & ...
            year(TIME_COMMON) >= year_start & ...
            year(TIME_COMMON) <= year_end;

        sal_month  = sal_all(idx_sal,:);
        cell_month = cell_all(idx_sal,:);
        sal_years  = year(TIME_COMMON(idx_sal));


        %% Same-calendar-month runoff
        Q = nan(length(sal_years),1);

        for iy = 1:length(sal_years)

            yy = sal_years(iy);

            idx_q = ...
                year(q_dates) == yy & ...
                month(q_dates) == mon;

            qtmp = discharge(b,idx_q);

            if numel(qtmp) == 1 && isfinite(qtmp)
                Q(iy) = qtmp;
            end

        end


        %% Detrend and standardize runoff
        if DO_DETREND

            Qanom = detrend_nan(Q);

        else

            Qanom = ...
                Q - mean(Q,'omitnan');

        end

        sigmaQ = std(Qanom,'omitnan');

        if ~isfinite(sigmaQ) || sigmaQ == 0
            continue
        end

        Qz = Qanom ./ sigmaQ;


        %% Distance bins
        for d = 1:nDist

            median_cells = ...
                median(cell_month(:,d),'omitnan');

            % Exclude unreliable coastal bins
            if ...
                    ~isfinite(median_cells) || ...
                    median_cells < MIN_VALID_OCEAN_CELLS

                continue
            end


            %% Salinity anomaly
            S = sal_month(:,d);

            if DO_DETREND

                Sanom = detrend_nan(S);

            else

                Sanom = ...
                    S - mean(S,'omitnan');

            end


            %% Valid pairs
            good = ...
                isfinite(Qz) & ...
                isfinite(Sanom);

            NPAIR(b,im,d) = sum(good);

            if sum(good) < MIN_VALID_YEARS
                continue
            end

            x = Qz(good);
            y = Sanom(good);


            %% Linear regression
            X = [ones(length(x),1), x];

            bb = X \ y;

            BETA(b,im,d) = bb(2);


            %% Pearson correlation
            [R,P] = corrcoef(x,y);

            RVAL(b,im,d) = R(1,2);
            PVAL(b,im,d) = P(1,2);

        end

    end
end


%% ========================================================================
%% 9. DETERMINE COMMON COLOR RANGE
%% ========================================================================

valid_beta = BETA(isfinite(BETA));

if isempty(valid_beta)
    error('No valid regression results were calculated.');
end


if isempty(COLOR_LIMIT)

    abs_beta = sort(abs(valid_beta));

    p = 0.98;

    ii = ceil(p*numel(abs_beta));

    ii = max(1,min(ii,numel(abs_beta)));

    cmax = abs_beta(ii);


    if cmax <= 0.30

        cmax = ceil(cmax*20)/20;

    elseif cmax <= 0.80

        cmax = ceil(cmax*10)/10;

    else

        cmax = ceil(cmax*5)/5;

    end

else

    cmax = COLOR_LIMIT;

end


fprintf('\nFinal colorbar range: %.2f to %.2f psu\n', ...
    -cmax,cmax);


%% ========================================================================
%% 10. CREATE FINAL FIGURE
%% ========================================================================

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 20.5 23.0]);


tl = tiledlayout( ...
    fig, ...
    3,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


colormap(fig,blue_white_red_soft(256));


for b = 1:nBasin

    ax = nexttile(tl,b);

    Z = squeeze(BETA(b,:,:));


    %% Heatmap
    himg = imagesc( ...
        ax, ...
        dist_centers, ...
        1:nMonth, ...
        Z);

    set(ax,'YDir','normal');


    % NaN / unreliable coastal bins shown as white
    set(himg, ...
        'AlphaData', ...
        double(isfinite(Z)));

    ax.Color = [1 1 1];


    %% Limits and ticks
    xlim(ax,[0 800]);

    ylim(ax,[0.5 nMonth+0.5]);


    xticks(ax,0:100:800);

    yticks(ax,1:nMonth);

    yticklabels(ax,month_labels);


    clim(ax,[-cmax cmax]);


    %% X label only on bottom row
    if b >= 5

        xlabel(ax, ...
            'Distance from river mouth (km)', ...
            'FontName','Times New Roman', ...
            'FontSize',10.5);

    else

        xlabel(ax,'');

    end


    %% Y label only on left column
    if mod(b,2) == 1

        ylabel(ax, ...
            'Month', ...
            'FontName','Times New Roman', ...
            'FontSize',10.5);

    else

        ylabel(ax,'');

    end


    %% Panel titles
    title(ax, ...
        [panel_letters{b} ' ' basin_names{b}], ...
        'FontName','Times New Roman', ...
        'FontSize',12.2, ...
        'FontWeight','bold');


    %% Axis style
    set(ax, ...
        'FontName','Times New Roman', ...
        'FontSize',9.8, ...
        'LineWidth',0.75, ...
        'TickDir','out', ...
        'TickLength',[0.012 0.012], ...
        'Box','on', ...
        'Layer','top');

    ax.XColor = [0.15 0.15 0.15];
    ax.YColor = [0.15 0.15 0.15];

end


%% ========================================================================
%% 11. SHARED COLORBAR
%% ========================================================================

cb = colorbar;

cb.Layout.Tile = 'east';

cb.FontName = 'Times New Roman';
cb.FontSize = 9.8;
cb.LineWidth = 0.70;
cb.TickDirection = 'out';


ylabel(cb, ...
    'Runoff-related salinity anomaly (psu)', ...
    'FontName','Times New Roman', ...
    'FontSize',10.5);


%% ========================================================================
%% 12. EXPORT PNG ONLY
%% ========================================================================

set(gcf,'InvertHardcopy','off');

exportgraphics( ...
    gcf, ...
    fullfile(cfg.outputDir,'Fig05_RunoffRelatedSalinity.png'), ...
    'Resolution',600);


fprintf('\nFigure 5 completed successfully.\n');


%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function y = detrend_nan(x)

x = double(x(:));

y = nan(size(x));

good = isfinite(x);

if sum(good) < 3
    return
end

t = (1:length(x))';

p = polyfit( ...
    t(good), ...
    x(good), ...
    1);

trend = polyval(p,t);

y(good) = ...
    x(good) - trend(good);

end



function cmap = blue_white_red_soft(n)

if nargin < 1
    n = 256;
end

blue_dark  = [0.12 0.30 0.72];
blue_light = [0.72 0.80 0.94];

white = [1.00 1.00 1.00];

red_light = [0.95 0.76 0.73];
red_dark  = [0.76 0.08 0.10];

n_half = floor(n/2);

n1 = round(n_half*0.55);
n2 = n_half-n1;

n3 = n2;
n4 = n-n1-n2-n3;

c1 = interp_colors(blue_dark,blue_light,n1);
c2 = interp_colors(blue_light,white,n2);
c3 = interp_colors(white,red_light,n3);
c4 = interp_colors(red_light,red_dark,n4);

cmap = [c1; c2; c3; c4];

end



function C = interp_colors(c1,c2,n)

if n <= 1
    C = c1;
    return
end

C = zeros(n,3);

for k = 1:3

    C(:,k) = ...
        linspace(c1(k),c2(k),n)';

end

end
