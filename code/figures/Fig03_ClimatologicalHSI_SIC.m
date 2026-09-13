%% Repository-relative setup
thisDir = fileparts(mfilename('fullpath'));
codeDir = fileparts(thisDir);
addpath(codeDir,fullfile(codeDir,'functions'));
cfg = project_config();
sys = river_shelf_systems();
close all;

%% ========================================================================
% Figure 3
% Climatological JAS distributions of:
%
% (a) Near-surface Haline Stratification Index (HSI), upper 15 m
% (b) Sea-ice concentration (SIC)
%
% Period: 1979-2024
%
% HSI:
%   loaded from pre-calculated intermediate file
%
% SIC:
%   calculated from original monthly SIC data
%% ========================================================================




%% ========================================================================
% 1. Paths
%% ========================================================================

% ---------------- HSI ----------------
HSI_file = ...
    fullfile(cfg.dataDir,'arctic_hsi_climatology_1979_2024.mat');

% ---------------- SIC ----------------
SIC_file = ...
    fullfile(cfg.dataDir,'sea_ice_monthly_1979_2024.mat');

SIC_coord_file = ...
    fullfile(cfg.dataDir,'sea_ice_grid_coordinates.mat');

% ---------------- Marginal-sea boundaries ----------------
sea_shp = ...
    fullfile(cfg.geoDir,'arctic_marginal_seas.shp');

% ---------------- Output ----------------
outDir = ...
    cfg.outputDir;

if ~exist(outDir,'dir')
    mkdir(outDir);
end


%% ========================================================================
% 2. Load HSI climatology
%% ========================================================================

H = load(HSI_file);

lon_HSI = double(H.lonGrid);
lat_HSI = double(H.latGrid);

HSI_plot = double(H.HSI_JAS_clim);


% ------------------------------------------------------------
% Check consistency
% ------------------------------------------------------------

if ~isequal(size(lon_HSI),size(lat_HSI),size(HSI_plot))
    error('HSI longitude, latitude, and HSI dimensions do not match.');
end


% Only retain Arctic region
HSI_plot(lat_HSI < 60) = NaN;


%% ========================================================================
% 3. HSI diagnostics and plotting range
%% ========================================================================

vals_HSI = HSI_plot( ...
    isfinite(HSI_plot) & ...
    lat_HSI >= 60 & ...
    HSI_plot >= 0);


fprintf('\n============================================\n');
fprintf('HSI statistics north of 60N\n');
fprintf('============================================\n');

fprintf('Minimum : %.6f s^-2\n', ...
    min(vals_HSI,[],'omitnan'));

fprintf('Mean    : %.6f s^-2\n', ...
    mean(vals_HSI,'omitnan'));

fprintf('P90     : %.6f s^-2\n', ...
    prctile(vals_HSI,90));

fprintf('P95     : %.6f s^-2\n', ...
    prctile(vals_HSI,95));

fprintf('P97.5   : %.6f s^-2\n', ...
    prctile(vals_HSI,97.5));

fprintf('P99     : %.6f s^-2\n', ...
    prctile(vals_HSI,99));

fprintf('Maximum : %.6f s^-2\n', ...
    max(vals_HSI,[],'omitnan'));


%% ------------------------------------------------------------------------
% HSI color limits
%
% Use P95 to increase visible spatial contrast.
%
% If later you feel too much area is saturated:
% change 95 -> 97.5 or 99.
%% ------------------------------------------------------------------------

HSI_cmin = 0;

HSI_cmax = ...
    prctile(vals_HSI,95);


% Round upward to a visually clean value
if HSI_cmax < 0.01

    HSI_cmax = ...
        ceil(HSI_cmax*1000)/1000;

else

    HSI_cmax = ...
        ceil(HSI_cmax*100)/100;

end


fprintf('\nHSI color range used:\n');
fprintf('0 - %.4f s^-2\n',HSI_cmax);


%% ========================================================================
% 4. Load SIC
%% ========================================================================

D = load(SIC_file);
C = load(SIC_coord_file);

sic = double(D.sic);

lon = double(C.lon);
lat = double(C.lat);


%% ------------------------------------------------------------------------
% Keep original SIC scaling
%% ------------------------------------------------------------------------

sic = sic / 0.004;


%% ========================================================================
% 5. SIC time processing
%% ========================================================================

yrs_all = 1979:2024;

n_years = numel(yrs_all);

mon_target = [7 8 9];    % JAS


[nx,ny,nt] = size(sic);


if nt ~= n_years*12

    error('SIC time dimension must be 46 × 12 = 552.');

end


%% ------------------------------------------------------------------------
% Convert SIC into 0-1 if needed
%% ------------------------------------------------------------------------

sic_valid = sic(isfinite(sic));


if max(sic_valid) > 1.5

    sic = sic ./ 100;

end


sic(sic < 0) = NaN;
sic(sic > 1) = NaN;


%% ------------------------------------------------------------------------
% Reshape:
%
% [nx × ny × 552]
% ->
% [nx × ny × 12 × 46]
%% ------------------------------------------------------------------------

sic_4d = ...
    reshape( ...
    sic, ...
    nx,ny,12,n_years);


%% ------------------------------------------------------------------------
% Annual JAS mean
%% ------------------------------------------------------------------------

SIC_JAS_Y = ...
    squeeze( ...
    mean( ...
    sic_4d(:,:,mon_target,:), ...
    3, ...
    'omitnan'));


%% ------------------------------------------------------------------------
% 1979-2024 climatological JAS SIC
%% ------------------------------------------------------------------------

SIC_JAS_clim = ...
    mean( ...
    SIC_JAS_Y, ...
    3, ...
    'omitnan');


SIC_plot = SIC_JAS_clim;


% Arctic only
SIC_plot(lat < 60) = NaN;


%% ========================================================================
% 6. General plotting settings
%% ========================================================================

fontAxis  = 6;
fontLabel = 7;
fontPanel = 10;


edgeColor = ...
    [0.30 0.50 0.70];


%% ------------------------------------------------------------------------
% Marginal-sea labels
%% ------------------------------------------------------------------------

lon_txt = [ ...
    -169, ...
    -138.6, ...
     165.7, ...
      80.0, ...
     121.3];


lat_txt = [ ...
    70, ...
    77, ...
    75.8, ...
    79.3, ...
    78.7];


name_txt = { ...
    'Chukchi', ...
    'Beaufort', ...
    'East Siberian', ...
    'Kara', ...
    'Laptev'};


%% ========================================================================
% 7. Load marginal-sea boundaries
%% ========================================================================

Ssea = [];


try

    Ssea = shaperead(sea_shp);

catch

    warning('Marginal-sea shapefile could not be loaded.');

end


%% ========================================================================
% 8. Colormaps
%% ========================================================================

nColor = 256;


%% ------------------------------------------------------------------------
% HSI colormap
%
% Light pink -> medium red -> deep red
%
% Important:
% lowest color is NOT pure white anymore.
%% ------------------------------------------------------------------------

anchorX = [ ...
    0.00 ...
    0.20 ...
    0.45 ...
    0.70 ...
    1.00];


anchorRGB = [ ...
    1.00 0.93 0.93; ...
    1.00 0.78 0.78; ...
    1.00 0.52 0.52; ...
    0.96 0.25 0.25; ...
    0.78 0.00 0.00];


xi = ...
    linspace(0,1,nColor);


cmap_HSI = ...
    interp1( ...
    anchorX, ...
    anchorRGB, ...
    xi, ...
    'linear');


%% ------------------------------------------------------------------------
% SIC colormap
% Keep original white-to-blue style
%% ------------------------------------------------------------------------

cmap_SIC = [ ...
    linspace(1.00,0.10,nColor)', ...
    linspace(1.00,0.30,nColor)', ...
    linspace(1.00,0.80,nColor)'];


%% ========================================================================
% 9. Create figure
%% ========================================================================

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 18.5 8.5]);


tl = tiledlayout( ...
    1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% ========================================================================
% 10. Panel (a): HSI
%% ========================================================================

ax1 = nexttile(tl,1);

axes(ax1);

hold(ax1,'on');


m_proj( ...
    'stereographic', ...
    'lat',90, ...
    'long',0, ...
    'radius',30);


h1 = ...
    m_pcolor( ...
    lon_HSI, ...
    lat_HSI, ...
    HSI_plot);


set(h1, ...
    'EdgeColor','none');


colormap(ax1,cmap_HSI);


caxis(ax1, ...
    [HSI_cmin HSI_cmax]);


%% ------------------------------------------------------------------------
% Bathymetry
%% ------------------------------------------------------------------------

safe_draw_isobath200( ...
    [0.25 0.25 0.25], ...
    0.8, ...
    '--');


%% ------------------------------------------------------------------------
% Coast + grid
%% ------------------------------------------------------------------------

safe_coast_and_grid(fontAxis);


%% ------------------------------------------------------------------------
% Deep contour
%% ------------------------------------------------------------------------

safe_deep_contour();


%% ------------------------------------------------------------------------
% Marginal-sea boundaries
%% ------------------------------------------------------------------------

safe_plot_sea_boundaries( ...
    Ssea, ...
    edgeColor);


%% ------------------------------------------------------------------------
% Sea names
%% ------------------------------------------------------------------------

safe_plot_names( ...
    lon_txt, ...
    lat_txt, ...
    name_txt);


%% ------------------------------------------------------------------------
% HSI colorbar
%% ------------------------------------------------------------------------

cb1 = colorbar(ax1);


cb1.FontName = ...
    'Times New Roman';


cb1.FontSize = ...
    fontAxis;


cb1.Label.String = ...
    'Climatological JAS HSI (s^{-2})';


cb1.Label.FontName = ...
    'Times New Roman';


cb1.Label.FontSize = ...
    fontLabel;


cb1.Label.FontWeight = ...
    'bold';


%% ------------------------------------------------------------------------
% HSI colorbar ticks
%% ------------------------------------------------------------------------

nTicksHSI = 6;


cb1.Ticks = ...
    linspace( ...
    HSI_cmin, ...
    HSI_cmax, ...
    nTicksHSI);


%% ------------------------------------------------------------------------
% Panel label
%% ------------------------------------------------------------------------

text( ...
    ax1, ...
    0.03, ...
    0.98, ...
    '(a)', ...
    'Units','normalized', ...
    'FontName','Times New Roman', ...
    'FontSize',fontPanel, ...
    'FontWeight','bold', ...
    'VerticalAlignment','top');


%% ========================================================================
% 11. Panel (b): SIC
%% ========================================================================

ax2 = nexttile(tl,2);

axes(ax2);

hold(ax2,'on');


m_proj( ...
    'stereographic', ...
    'lat',90, ...
    'long',0, ...
    'radius',30);


h2 = ...
    m_pcolor( ...
    lon, ...
    lat, ...
    SIC_plot);


set(h2, ...
    'EdgeColor','none');


colormap(ax2,cmap_SIC);


caxis(ax2,[0 1]);


%% ------------------------------------------------------------------------
% Bathymetry
%% ------------------------------------------------------------------------

safe_draw_isobath200( ...
    [0.25 0.25 0.25], ...
    0.8, ...
    '--');


%% ------------------------------------------------------------------------
% Coast + grid
%% ------------------------------------------------------------------------

safe_coast_and_grid(fontAxis);


%% ------------------------------------------------------------------------
% Deep contour
%% ------------------------------------------------------------------------

safe_deep_contour();


%% ------------------------------------------------------------------------
% Marginal-sea boundaries
%% ------------------------------------------------------------------------

safe_plot_sea_boundaries( ...
    Ssea, ...
    edgeColor);


%% ------------------------------------------------------------------------
% Sea names
%% ------------------------------------------------------------------------

safe_plot_names( ...
    lon_txt, ...
    lat_txt, ...
    name_txt);


%% ------------------------------------------------------------------------
% SIC colorbar
%% ------------------------------------------------------------------------

cb2 = colorbar(ax2);


cb2.FontName = ...
    'Times New Roman';


cb2.FontSize = ...
    fontAxis;


cb2.Label.String = ...
    'Climatological JAS SIC';


cb2.Label.FontName = ...
    'Times New Roman';


cb2.Label.FontSize = ...
    fontLabel;


cb2.Label.FontWeight = ...
    'bold';


cb2.Ticks = ...
    0:0.1:1;


%% ------------------------------------------------------------------------
% Panel label
%% ------------------------------------------------------------------------

text( ...
    ax2, ...
    0.03, ...
    0.98, ...
    '(b)', ...
    'Units','normalized', ...
    'FontName','Times New Roman', ...
    'FontSize',fontPanel, ...
    'FontWeight','bold', ...
    'VerticalAlignment','top');


%% ========================================================================
% 12. Final appearance
%% ========================================================================

set(fig, ...
    'InvertHardcopy','off');


exportgraphics( ...
    gcf, ...
    fullfile(cfg.outputDir,'Fig03_ClimatologicalHSI_SIC.png'), ...
    'Resolution',600);


%% ========================================================================
% LOCAL FUNCTIONS
%% ========================================================================


function safe_draw_isobath200(c,lw,ls)

    try

        m_etopo2( ...
            'contour', ...
            [-200 -200], ...
            'color',c, ...
            'linewidth',lw, ...
            'linestyle',ls);

    catch

    end

end


function safe_coast_and_grid(font)

    m_coast( ...
        'patch', ...
        [229 203 139]/256, ...
        'edgecolor','k');


    m_grid( ...
        'box','fancy', ...
        'tickdir','in', ...
        'gridlines','no', ...
        'fontsize',font, ...
        'xticklabel',[], ...
        'yticklabel',[]);

end


function safe_deep_contour()

    try

        m_etopo2( ...
            'contour', ...
            [-1000 -1000], ...
            'color',[0.30 0.30 0.30], ...
            'LineWidth',1.0);

    catch

    end

end


function safe_plot_sea_boundaries( ...
    Ssea, ...
    edgeColor)

    if isempty(Ssea)
        return
    end


    for ii = 1:length(Ssea)

        try

            m_line( ...
                Ssea(ii).X, ...
                Ssea(ii).Y, ...
                'Color',edgeColor, ...
                'LineWidth',0.8);

        catch

        end

    end

end


function safe_plot_names( ...
    lon, ...
    lat, ...
    name)

    for i = 1:length(name)

        m_text( ...
            lon(i), ...
            lat(i), ...
            name{i}, ...
            'FontName','Times New Roman', ...
            'FontSize',5, ...
            'FontWeight','bold', ...
            'HorizontalAlignment','center');

    end

end
