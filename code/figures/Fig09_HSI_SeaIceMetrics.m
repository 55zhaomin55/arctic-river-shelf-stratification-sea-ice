%% Repository-relative setup
thisDir = fileparts(mfilename('fullpath'));
codeDir = fileparts(thisDir);
addpath(codeDir,fullfile(codeDir,'functions'));
cfg = project_config();
sys = river_shelf_systems();
close all;

%% ========================================================================
% Figure 9
% Relationships of JAS HSI anomalies with:
%
% (a) Open-water fraction (OWF)
% (b) Southern sea-ice-edge latitude (SIC_edge)
%
% Period: 1979-2024
%
% FINAL VERSION
%
% HSI:
%   directly loaded from annual_river_shelf_metrics_1979_2024.mat
%
% OWF:
%   fraction of sampling domain with SIC < 0.15
%
% SIC_edge:
%   southernmost latitude of the 15% SIC contour
%   within each 450-km sampling domain
%
% Correlation:
%   Pearson correlation between detrended annual anomalies
%   ordinary two-sided p value
%% ========================================================================




%% ========================================================================
% 1. PATHS
%% ========================================================================

HSI_file = ...
    fullfile(cfg.dataDir,'annual_river_shelf_metrics_1979_2024.mat');

sicFile = ...
    fullfile(cfg.dataDir,'sea_ice_monthly_1979_2024.mat');

sicLonLatFile = ...
    fullfile(cfg.dataDir,'sea_ice_grid_coordinates.mat');

outDir = cfg.outputDir;

if ~exist(outDir,'dir')
    mkdir(outDir);
end


%% ========================================================================
% 2. LOAD FINAL HSI DATA
%% ========================================================================

D = load(HSI_file);

if ~isfield(D,'HSI_JAS_dt')
    error('Annual metrics file does not contain HSI_JAS_dt.');
end

if ~isfield(D,'years')
    error('Annual metrics file does not contain years.');
end


HSI_dt = ...
    double(D.HSI_JAS_dt);

years = ...
    double(D.years(:));


% ------------------------------------------------------------
% Convert to 46 × 6
% ------------------------------------------------------------

if size(HSI_dt,1) == 6 && ...
   size(HSI_dt,2) == 46

    HSI_dt = HSI_dt';

end


assert( ...
    isequal(size(HSI_dt),[46 6]), ...
    'HSI_JAS_dt must be 46 × 6.');


nYear = ...
    length(years);

nBasin = ...
    6;


%% ------------------------------------------------------------------------
% Saved final regional SIC from annual_river_shelf_metrics_1979_2024.mat
% Used only for consistency checking
%% ------------------------------------------------------------------------

SIC_saved = [];

if isfield(D,'SIC_JAS')

    SIC_saved = ...
        double(D.SIC_JAS);

    if size(SIC_saved,1) == 6
        SIC_saved = SIC_saved';
    end

end


%% ========================================================================
% 3. SYSTEM NAMES AND FINAL 450-km SAMPLING DOMAINS
%% ========================================================================

Rname = { ...
    'Mackenzie-Beaufort', ...
    'Yukon-Chukchi', ...
    'Kolyma-East Siberian', ...
    'Lena-Laptev', ...
    'Yenisei-Kara', ...
    'Ob-Kara'};


% ========================================================================
% IMPORTANT
%
% Sampling centers are centralized in river_shelf_systems.m and are
% shared with the final domain-sensitivity and sampling-domain analyses.
% ========================================================================

mouth_lon = sys.mouthLon;
mouth_lat = sys.mouthLat;


radius_km = ...
    450;


Re_km = ...
    6371;


% Unified SIC threshold
SIC_THRESHOLD = ...
    0.15;


%% ========================================================================
% 4. LOAD SIC
%% ========================================================================

S = load(sicFile);
C = load(sicLonLatFile);


sic = ...
    double(S.sic);

lon = ...
    double(C.lon);

lat = ...
    double(C.lat);


%% ------------------------------------------------------------------------
% Keep the same SIC decoding used in your previous scripts
%% ------------------------------------------------------------------------

sic = ...
    sic / 0.004;


sic(sic >= 1.004) = ...
    NaN;

sic(sic < 0 | sic > 1) = ...
    NaN;


[nx,ny,nt] = ...
    size(sic);


if nt ~= nYear*12

    error( ...
        'SIC time dimension must be 46 × 12 = 552.');

end


%% ========================================================================
% 5. CALCULATE ANNUAL JAS SIC FIELDS
%% ========================================================================

mon_target = ...
    [7 8 9];


SIC_JAS_grid = ...
    nan(nx,ny,nYear,'single');


for iy = 1:nYear

    idxJAS = ...
        (iy-1)*12 + ...
        mon_target;


    SIC_JAS_grid(:,:,iy) = ...
        single( ...
        mean( ...
        sic(:,:,idxJAS), ...
        3, ...
        'omitnan'));

end


clear sic S C


%% ========================================================================
% 6. AREA WEIGHTS
%% ========================================================================

W = ...
    cosd(lat);


%% ========================================================================
% 7. CREATE FINAL 450-km DOMAIN MASKS
%% ========================================================================

MASK = ...
    false(nx,ny,nBasin);


for ib = 1:nBasin

    Dkm = ...
        gc_distance_km( ...
        lon, ...
        lat, ...
        mouth_lon(ib), ...
        mouth_lat(ib), ...
        Re_km);


    MASK(:,:,ib) = ...
        isfinite(Dkm) & ...
        Dkm <= radius_km;

end


%% ========================================================================
% 8. CALCULATE OWF, SIC_EDGE AND DOMAIN-MEAN SIC
%% ========================================================================

OWF_raw = ...
    nan(nYear,nBasin);


SICedge_raw = ...
    nan(nYear,nBasin);


SIC_domain = ...
    nan(nYear,nBasin);


for ib = 1:nBasin

    M = ...
        MASK(:,:,ib);


    for iy = 1:nYear

        F = ...
            double(SIC_JAS_grid(:,:,iy));


        %% ------------------------------------------------------------
        % Domain-mean SIC
        % only for checking consistency with annual_river_shelf_metrics_1979_2024.mat
        %% ------------------------------------------------------------

        validMean = ...
            M & ...
            isfinite(F) & ...
            isfinite(W);


        if nnz(validMean) >= 3

            ww = ...
                W(validMean);

            ff = ...
                F(validMean);


            SIC_domain(iy,ib) = ...
                sum(ff.*ww) / ...
                sum(ww);

        end


        %% ------------------------------------------------------------
        % Open-water fraction
        %
        % Area-weighted fraction where SIC < 0.15
        %% ------------------------------------------------------------

        valid = ...
            M & ...
            isfinite(F) & ...
            isfinite(W);


        if nnz(valid) >= 3

            ww = ...
                W(valid);

            ff = ...
                F(valid);


            OWF_raw(iy,ib) = ...
                sum( ...
                ww .* ...
                double(ff < SIC_THRESHOLD)) / ...
                sum(ww);

        end


        %% ------------------------------------------------------------
        % Southern sea-ice-edge latitude
        %
        % Southernmost latitude of 15% SIC contour
        %% ------------------------------------------------------------

        SICedge_raw(iy,ib) = ...
            local_southern_ice_edge_lat( ...
            F, ...
            lon, ...
            lat, ...
            M, ...
            SIC_THRESHOLD, ...
            mouth_lon(ib));

    end

end


%% ========================================================================
% 9. CHECK DOMAIN CONSISTENCY AGAINST FINAL annual_river_shelf_metrics_1979_2024.mat
%
% This is important because Figure 9 should use the same 450-km domains
% as Figure 8.
%% ========================================================================

if ~isempty(SIC_saved) && ...
   isequal(size(SIC_saved),size(SIC_domain))

    fprintf('\n============================================================\n');
    fprintf('DOMAIN CONSISTENCY CHECK\n');
    fprintf('Computed 450-km SIC vs SIC_JAS saved in annual metrics file\n');
    fprintf('============================================================\n');


    for ib = 1:nBasin

        x = ...
            SIC_domain(:,ib);

        y = ...
            SIC_saved(:,ib);


        ok = ...
            isfinite(x) & ...
            isfinite(y);


        if nnz(ok) >= 3

            rr = ...
                corr( ...
                x(ok), ...
                y(ok));


            rmse = ...
                sqrt( ...
                mean( ...
                (x(ok)-y(ok)).^2));


            fprintf( ...
                '%-24s  R = %.4f   RMSE = %.5f\n', ...
                Rname{ib}, ...
                rr, ...
                rmse);

        end

    end


    fprintf('\n');


    % Optional warning
    for ib = 1:nBasin

        x = SIC_domain(:,ib);
        y = SIC_saved(:,ib);

        ok = isfinite(x) & isfinite(y);

        if nnz(ok) >= 3

            rr = corr(x(ok),y(ok));

            if rr < 0.95

                warning( ...
                    ['%s: regional SIC does not closely reproduce ' ...
                     'SIC_JAS in annual metrics file. Check river-mouth coordinates.'], ...
                    Rname{ib});

            end

        end

    end

end


%% ========================================================================
% 10. DETREND OWF AND SIC_EDGE
%
% HSI is already detrended in HSI_JAS_dt.
%% ========================================================================

OWF_dt = ...
    nan(size(OWF_raw));


SICedge_dt = ...
    nan(size(SICedge_raw));


for ib = 1:nBasin

    OWF_dt(:,ib) = ...
        local_detrend_anom( ...
        OWF_raw(:,ib), ...
        years);


    SICedge_dt(:,ib) = ...
        local_detrend_anom( ...
        SICedge_raw(:,ib), ...
        years);

end


%% ========================================================================
% 11. CALCULATE FINAL CORRELATIONS
%% ========================================================================

R_OWF = ...
    nan(nBasin,1);

P_OWF = ...
    nan(nBasin,1);

N_OWF = ...
    nan(nBasin,1);


R_EDGE = ...
    nan(nBasin,1);

P_EDGE = ...
    nan(nBasin,1);

N_EDGE = ...
    nan(nBasin,1);


for ib = 1:nBasin

    %% --------------------------------------------------------
    % HSI - OWF
    %% --------------------------------------------------------

    x = ...
        HSI_dt(:,ib);

    y = ...
        OWF_dt(:,ib);


    [R_OWF(ib), ...
     P_OWF(ib), ...
     N_OWF(ib)] = ...
        local_corr( ...
        x,y);


    %% --------------------------------------------------------
    % HSI - SIC_edge
    %% --------------------------------------------------------

    x = ...
        HSI_dt(:,ib);

    y = ...
        SICedge_dt(:,ib);


    [R_EDGE(ib), ...
     P_EDGE(ib), ...
     N_EDGE(ib)] = ...
        local_corr( ...
        x,y);

end


%% ========================================================================
% 12. DISPLAY RESULTS
%% ========================================================================

T = table( ...
    Rname', ...
    R_OWF, ...
    P_OWF, ...
    N_OWF, ...
    R_EDGE, ...
    P_EDGE, ...
    N_EDGE, ...
    ...
    'VariableNames',{ ...
    'System', ...
    'R_HSI_OWF', ...
    'p_HSI_OWF', ...
    'N_HSI_OWF', ...
    'R_HSI_SICedge', ...
    'p_HSI_SICedge', ...
    'N_HSI_SICedge'});


fprintf('\n============================================================\n');
fprintf('FINAL FIGURE 9 RESULTS\n');
fprintf('============================================================\n');

disp(T);



%% ========================================================================
% 14. PLOTTING RANGE
%% ========================================================================

% ------------------------------------------------------------------------
% HSI is now ~10^-3 s^-2.
%
% Multiply by 1000 ONLY for plotting.
% Correlations above use original HSI units.
% ------------------------------------------------------------------------

HSI_SCALE = ...
    1e3;


HSI_plot = ...
    HSI_dt * ...
    HSI_SCALE;


x_all = ...
    HSI_plot(:);


x_all = ...
    x_all(isfinite(x_all));


xmax = ...
    max(abs(x_all));


xmax = ...
    ceil(xmax/0.1)*0.1;


xlim_common = ...
    [-xmax xmax];


%% ------------------------------------------------------------------------
% OWF
%% ------------------------------------------------------------------------

y1 = ...
    OWF_dt(:);


y1 = ...
    y1(isfinite(y1));


if isempty(y1)

    ylim_left = ...
        [-0.1 0.1];

else

    ymax1 = ...
        max(abs(y1));

    ymax1 = ...
        max(ymax1,0.05);

    ymax1 = ...
        ceil(ymax1/0.05)*0.05;

    ylim_left = ...
        [-ymax1 ymax1];

end


%% ------------------------------------------------------------------------
% SIC edge
%% ------------------------------------------------------------------------

y2 = ...
    SICedge_dt(:);


y2 = ...
    y2(isfinite(y2));


if isempty(y2)

    ylim_right = ...
        [-2 2];

else

    ymax2 = ...
        max(abs(y2));

    ymax2 = ...
        max(ymax2,1);

    ymax2 = ...
        ceil(ymax2/0.5)*0.5;

    ylim_right = ...
        [-ymax2 ymax2];

end


%% ========================================================================
% 15. FIGURE
%% ========================================================================

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 17.5 20.5]);


tiledlayout( ...
    6,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


fontAxis  = 7.2;
fontLabel = 8.0;
fontTitle = 8.8;
fontStat  = 6.7;


markerSize = 13;
fitLW = 0.9;


cPoint = ...
    [0.78 0.05 0.03];


cFit = ...
    [0.25 0.25 0.25];


cZero = ...
    [0.70 0.70 0.70];


cGrid = ...
    [0.78 0.78 0.78];


%% ========================================================================
% 16. PLOT ALL SIX SYSTEMS
%% ========================================================================

for ib = 1:nBasin


    %% ================================================================
    % LEFT: HSI vs OWF
    %% ================================================================

    ax1 = ...
        nexttile( ...
        (ib-1)*2+1);


    hold(ax1,'on');


    x = ...
        HSI_plot(:,ib);


    y = ...
        OWF_dt(:,ib);


    ok = ...
        isfinite(x) & ...
        isfinite(y);


    xp = ...
        x(ok);

    yp = ...
        y(ok);


    [slope,intercept] = ...
        local_linear_fit( ...
        xp,yp);


    %% zero lines

    xline( ...
        ax1, ...
        0, ...
        '-', ...
        'Color',cZero, ...
        'LineWidth',0.55);


    yline( ...
        ax1, ...
        0, ...
        '-', ...
        'Color',cZero, ...
        'LineWidth',0.55);


    %% regression line

    if isfinite(slope)

        xx = ...
            linspace( ...
            xlim_common(1), ...
            xlim_common(2), ...
            100);


        yy = ...
            intercept + ...
            slope*xx;


        plot( ...
            ax1, ...
            xx,yy, ...
            '--', ...
            'Color',cFit, ...
            'LineWidth',fitLW);

    end


    %% scatter

    hs = ...
        scatter( ...
        ax1, ...
        xp,yp, ...
        markerSize, ...
        'MarkerFaceColor',cPoint, ...
        'MarkerEdgeColor','k', ...
        'LineWidth',0.35);


    try

        hs.MarkerFaceAlpha = ...
            0.68;

        hs.MarkerEdgeAlpha = ...
            0.85;

    catch

    end


    %% statistics

    statText = ...
        local_stat_text( ...
        R_OWF(ib), ...
        P_OWF(ib), ...
        N_OWF(ib));


    text( ...
        ax1, ...
        0.035, ...
        0.10, ...
        statText, ...
        'Units','normalized', ...
        'FontSize',fontStat, ...
        'FontName','Times New Roman', ...
        'FontWeight','bold', ...
        'VerticalAlignment','bottom');


    %% system name

    text( ...
        ax1, ...
        0.97, ...
        0.86, ...
        Rname{ib}, ...
        'Units','normalized', ...
        'FontSize',fontTitle, ...
        'FontName','Times New Roman', ...
        'FontWeight','bold', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','top');


    if ib == 1

        title( ...
            ax1, ...
            '(a) Open-water fraction', ...
            'FontName','Times New Roman', ...
            'FontSize',fontTitle, ...
            'FontWeight','bold');

    end


    if ib == 6

        xlabel( ...
            ax1, ...
            'Detrended HSI anomaly (10^{-3} s^{-2})', ...
            'FontSize',fontLabel, ...
            'FontName','Times New Roman', ...
            'FontWeight','bold', ...
            'Interpreter','tex');

    else

        ax1.XTickLabel = [];

    end


    ylabel( ...
        ax1, ...
        'OWF anomaly', ...
        'FontSize',fontLabel, ...
        'FontName','Times New Roman', ...
        'FontWeight','bold');


    xlim(ax1,xlim_common);
    ylim(ax1,ylim_left);


    box(ax1,'on');
    grid(ax1,'on');


    set( ...
        ax1, ...
        'FontSize',fontAxis, ...
        'FontName','Times New Roman', ...
        'FontWeight','bold', ...
        'LineWidth',0.75, ...
        'TickDir','out', ...
        'Layer','top', ...
        'GridColor',cGrid, ...
        'GridAlpha',0.35);


    %% ================================================================
    % RIGHT: HSI vs southern SIC edge
    %% ================================================================

    ax2 = ...
        nexttile( ...
        (ib-1)*2+2);


    hold(ax2,'on');


    x = ...
        HSI_plot(:,ib);


    y = ...
        SICedge_dt(:,ib);


    ok = ...
        isfinite(x) & ...
        isfinite(y);


    xp = ...
        x(ok);

    yp = ...
        y(ok);


    [slope,intercept] = ...
        local_linear_fit( ...
        xp,yp);


    %% zero lines

    xline( ...
        ax2, ...
        0, ...
        '-', ...
        'Color',cZero, ...
        'LineWidth',0.55);


    yline( ...
        ax2, ...
        0, ...
        '-', ...
        'Color',cZero, ...
        'LineWidth',0.55);


    %% regression line

    if isfinite(slope)

        xx = ...
            linspace( ...
            xlim_common(1), ...
            xlim_common(2), ...
            100);


        yy = ...
            intercept + ...
            slope*xx;


        plot( ...
            ax2, ...
            xx,yy, ...
            '--', ...
            'Color',cFit, ...
            'LineWidth',fitLW);

    end


    %% scatter

    hs = ...
        scatter( ...
        ax2, ...
        xp,yp, ...
        markerSize, ...
        'MarkerFaceColor',cPoint, ...
        'MarkerEdgeColor','k', ...
        'LineWidth',0.35);


    try

        hs.MarkerFaceAlpha = ...
            0.68;

        hs.MarkerEdgeAlpha = ...
            0.85;

    catch

    end


    %% statistics

    statText = ...
        local_stat_text( ...
        R_EDGE(ib), ...
        P_EDGE(ib), ...
        N_EDGE(ib));


    text( ...
        ax2, ...
        0.035, ...
        0.10, ...
        statText, ...
        'Units','normalized', ...
        'FontSize',fontStat, ...
        'FontName','Times New Roman', ...
        'FontWeight','bold', ...
        'VerticalAlignment','bottom');


    if ib == 1

        title( ...
            ax2, ...
            '(b) Southern sea-ice-edge latitude', ...
            'FontName','Times New Roman', ...
            'FontSize',fontTitle, ...
            'FontWeight','bold');

    end


    if ib == 6

        xlabel( ...
            ax2, ...
            'Detrended HSI anomaly (10^{-3} s^{-2})', ...
            'FontSize',fontLabel, ...
            'FontName','Times New Roman', ...
            'FontWeight','bold', ...
            'Interpreter','tex');

    else

        ax2.XTickLabel = [];

    end


    ylabel( ...
        ax2, ...
        'SIC_{edge} anomaly (^{\circ}N)', ...
        'FontSize',fontLabel, ...
        'FontName','Times New Roman', ...
        'FontWeight','bold', ...
        'Interpreter','tex');


    xlim(ax2,xlim_common);
    ylim(ax2,ylim_right);


    box(ax2,'on');
    grid(ax2,'on');


    set( ...
        ax2, ...
        'FontSize',fontAxis, ...
        'FontName','Times New Roman', ...
        'FontWeight','bold', ...
        'LineWidth',0.75, ...
        'TickDir','out', ...
        'Layer','top', ...
        'GridColor',cGrid, ...
        'GridAlpha',0.35);

end


%% ========================================================================
% 17. EXPORT
%% ========================================================================
drawnow;

set(fig,'InvertHardcopy','off');

exportgraphics( ...
    gcf, ...
    fullfile(cfg.outputDir,'Fig09_HSI_SeaIceMetrics.png'), ...
    'Resolution',600);




%% ========================================================================
% LOCAL FUNCTIONS
%% ========================================================================


function y = local_detrend_anom(x,years)

    x = ...
        double(x(:));

    years = ...
        double(years(:));


    y = ...
        nan(size(x));


    ok = ...
        isfinite(x) & ...
        isfinite(years);


    if nnz(ok) < 3
        return
    end


    p = ...
        polyfit( ...
        years(ok), ...
        x(ok), ...
        1);


    trend = ...
        polyval( ...
        p, ...
        years(ok));


    y(ok) = ...
        x(ok)-trend;

end


function [r,p,N] = local_corr(x,y)

    x = ...
        double(x(:));

    y = ...
        double(y(:));


    ok = ...
        isfinite(x) & ...
        isfinite(y);


    x = x(ok);
    y = y(ok);


    N = ...
        numel(x);


    r = NaN;
    p = NaN;


    if N < 3 || ...
       std(x) == 0 || ...
       std(y) == 0

        return

    end


    [r,p] = ...
        corr( ...
        x,y, ...
        'Rows','complete');

end


function [slope,intercept] = local_linear_fit(x,y)

    ok = ...
        isfinite(x) & ...
        isfinite(y);


    x = x(ok);
    y = y(ok);


    slope = NaN;
    intercept = NaN;


    if numel(x) < 3 || ...
       std(x) == 0

        return

    end


    pp = ...
        polyfit(x,y,1);


    slope = ...
        pp(1);


    intercept = ...
        pp(2);

end


function txt = local_stat_text(r,p,N)

    if ~isfinite(r)

        txt = ...
            sprintf( ...
            'R = --, p = --, N = %d', ...
            N);

        return

    end


    if ~isfinite(p)

        pstr = ...
            '--';

    elseif p < 0.001

        pstr = ...
            '<0.001';

    else

        pstr = ...
            sprintf( ...
            '%.2f', ...
            p);

    end


    txt = ...
        sprintf( ...
        'R = %.2f, p = %s, N = %d', ...
        r, ...
        pstr, ...
        N);

end


function d = gc_distance_km( ...
    LON,LAT,lon0,lat0,Re)

    phi = ...
        deg2rad(LAT);


    lam = ...
        deg2rad(LON);


    phi0 = ...
        deg2rad(lat0);


    lam0 = ...
        deg2rad(lon0);


    cosang = ...
        sin(phi).*sin(phi0) + ...
        cos(phi).*cos(phi0).* ...
        cos(lam-lam0);


    cosang = ...
        max(-1,min(1,cosang));


    d = ...
        acos(cosang)*Re;

end


function edgeLat = local_southern_ice_edge_lat( ...
    F,lon,lat,M,threshold,lon0)

% ========================================================================
% Southernmost latitude of the threshold SIC contour
% within the selected circular sampling domain.
%
% Unlike the old version:
%   - does NOT select the longest contour and average its latitude
%   - uses the SOUTHERNMOST latitude of all valid 15% contour segments
%
% Longitude is converted relative to the basin center to avoid
% dateline problems in Yukon/Chukchi and Kolyma/East Siberian.
%% ========================================================================

    edgeLat = ...
        NaN;


    valid = ...
        M & ...
        isfinite(F) & ...
        isfinite(lon) & ...
        isfinite(lat);


    if nnz(valid) < 30
        return
    end


    sicv = ...
        F(valid);


    % Need values on both sides of the threshold
    if min(sicv,[],'omitnan') >= threshold || ...
       max(sicv,[],'omitnan') < threshold

        return

    end


    lonv = ...
        lon(valid);


    latv = ...
        lat(valid);


    %% --------------------------------------------------------
    % Local longitude relative to river-mouth longitude
    % Range approximately [-180,180]
    %% --------------------------------------------------------

    lonRel = ...
        mod( ...
        lonv-lon0+180, ...
        360)-180;


    %% --------------------------------------------------------
    % Interpolation grid
    %% --------------------------------------------------------

    lonMin = ...
        min(lonRel);

    lonMax = ...
        max(lonRel);


    latMin = ...
        min(latv);

    latMax = ...
        max(latv);


    if lonMax <= lonMin || ...
       latMax <= latMin

        return

    end


    [LONi,LATi] = ...
        meshgrid( ...
        linspace(lonMin,lonMax,220), ...
        linspace(latMin,latMax,220));


    Fi = ...
        griddata( ...
        lonRel, ...
        latv, ...
        sicv, ...
        LONi, ...
        LATi, ...
        'natural');


    %% --------------------------------------------------------
    % Find threshold contour
    %% --------------------------------------------------------

    C = ...
        contourc( ...
        LONi(1,:), ...
        LATi(:,1), ...
        Fi, ...
        [threshold threshold]);


    if isempty(C)
        return
    end


    %% --------------------------------------------------------
    % Read all contour segments
    %% --------------------------------------------------------

    allLat = [];


    kk = 1;


    while kk < size(C,2)

        level = ...
            C(1,kk); %#ok<NASGU>


        npt = ...
            C(2,kk);


        if npt <= 0 || ...
           kk+npt > size(C,2)

            break

        end


        xy = ...
            C(:,kk+1:kk+npt);


        segLat = ...
            xy(2,:);


        segLat = ...
            segLat(isfinite(segLat));


        if ~isempty(segLat)

            allLat = ...
                [allLat segLat]; %#ok<AGROW>

        end


        kk = ...
            kk+npt+1;

    end


    if isempty(allLat)
        return
    end


    %% --------------------------------------------------------
    % SOUTHERNMOST latitude
    %% --------------------------------------------------------

    edgeLat = ...
        min(allLat);

end
