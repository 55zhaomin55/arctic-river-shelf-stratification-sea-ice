%% Repository-relative setup
thisDir = fileparts(mfilename('fullpath'));
codeDir = fileparts(thisDir);
addpath(codeDir,fullfile(codeDir,'functions'));
cfg = project_config();
sys = river_shelf_systems();
close all;

%% ========================================================================
%% 1. PATHS
%% ========================================================================

data_dir = cfg.dataDir;

cache_file = fullfile( ...
    data_dir, ...
    'marginal_sea_salinity_profiles_1979_2024.mat');

out_png = ...
    fullfile(cfg.outputDir,'Fig02_VerticalSalinityProfiles.png');

if ~exist(cache_file,'file')
    error('Cannot find intermediate file:\n%s',cache_file);
end


%% ========================================================================
%% 2. LOAD PROCESSED DATA
%% ========================================================================

C = load(cache_file);

SAL_MONTHLY = C.SAL_MONTHLY;
TIME_COMMON = C.TIME_COMMON;
depth       = C.depth;
sea_names   = C.sea_names;

% Prefer area-weighted numerator / denominator
if isfield(C,'NUM_MONTHLY') && isfield(C,'DEN_MONTHLY')

    NUM_MONTHLY = C.NUM_MONTHLY;
    DEN_MONTHLY = C.DEN_MONTHLY;

    USE_WEIGHTED_CACHE = true;

else

    USE_WEIGHTED_CACHE = false;

end

fprintf('\nIntermediate ORAS5 data loaded successfully.\n');

fprintf('Period: %s to %s\n', ...
    datestr(TIME_COMMON(1),'yyyy-mm'), ...
    datestr(TIME_COMMON(end),'yyyy-mm'));

fprintf('SAL_MONTHLY size: ');
disp(size(SAL_MONTHLY));


%% ========================================================================
%% 3. SETTINGS
%% ========================================================================

JAS_MONTHS = [7 8 9];

MAX_PLOT_DEPTH = 60;

panel_names = {
    'Beaufort Sea'
    'Chukchi Sea'
    'East Siberian Sea'
    'Laptev Sea'
    'Kara Sea'
    'Five-sea composite'
    };

panel_letters = {
    '(a)','(b)','(c)', ...
    '(d)','(e)','(f)'
    };

nSea = 5;


%% ---------------- Figure appearance -------------------------------------

FONT_NAME = 'Times New Roman';

FONT_AXIS  = 8.2;     % tick labels
FONT_LABEL = 8.8;     % axis labels
FONT_TITLE = 11.5;    % panel titles
FONT_TEXT  = 8.5;     % annotations

LINE_COLOR = [0.055 0.30 0.62];
LINE_WIDTH = 2.00;

SHADE_COLOR = [0.91 0.95 0.98];
GRID_COLOR  = [0.82 0.82 0.82];
REF_COLOR   = [0.52 0.52 0.52];


%% ========================================================================
%% 4. CALCULATE JAS CLIMATOLOGICAL PROFILES
%% ========================================================================

jas_idx = ismember(month(TIME_COMMON),JAS_MONTHS);

nDepth = length(depth);

PROFILE_SEA = nan(nSea,nDepth);


%% ------------------------------------------------------------------------
% Area-weighted JAS climatology
%% ------------------------------------------------------------------------

if USE_WEIGHTED_CACHE

    NUM_JAS = nan(nSea,nDepth);
    DEN_JAS = nan(nSea,nDepth);

    for s = 1:nSea

        num = squeeze(NUM_MONTHLY(s,jas_idx,:));
        den = squeeze(DEN_MONTHLY(s,jas_idx,:));

        NUM_JAS(s,:) = ...
            sum(num,1,'omitnan');

        DEN_JAS(s,:) = ...
            sum(den,1,'omitnan');

        PROFILE_SEA(s,:) = ...
            NUM_JAS(s,:) ./ DEN_JAS(s,:);

    end


    % Five-sea area-weighted composite
    NUM_COMPOSITE = ...
        sum(NUM_JAS,1,'omitnan');

    DEN_COMPOSITE = ...
        sum(DEN_JAS,1,'omitnan');

    PROFILE_COMPOSITE = ...
        NUM_COMPOSITE ./ DEN_COMPOSITE;


else

    %% --------------------------------------------------------------------
    % Fallback if numerator / denominator were not saved
    %% --------------------------------------------------------------------

    for s = 1:nSea

        tmp = ...
            squeeze(SAL_MONTHLY(s,jas_idx,:));

        PROFILE_SEA(s,:) = ...
            mean(tmp,1,'omitnan');

    end

    PROFILE_COMPOSITE = ...
        mean(PROFILE_SEA,1,'omitnan');

end


PROFILE_ALL = [
    PROFILE_SEA
    PROFILE_COMPOSITE
    ];


%% ========================================================================
%% 5. SELECT DEPTH RANGE
%% ========================================================================

plot_idx = ...
    depth <= MAX_PLOT_DEPTH;

dep_plot = ...
    depth(plot_idx);

PLOT_PROFILE = ...
    PROFILE_ALL(:,plot_idx);


%% ========================================================================
%% 6. CALCULATE NEAR-SURFACE SALINITY DIFFERENCE
%% ========================================================================
%
% Delta S = S(15 m) - S(shallowest available ORAS5 level)
%
% ========================================================================

dS15 = nan(6,1);

for p = 1:6

    prof = ...
        squeeze(PLOT_PROFILE(p,:))';

    good = ...
        isfinite(dep_plot) & ...
        isfinite(prof);

    z  = dep_plot(good);
    ss = prof(good);

    if length(z) < 2
        continue
    end

    % Sort by depth
    [z,sortID] = sort(z);
    ss = ss(sortID);

    % Shallowest model level
    S_surface = ss(1);

    % Interpolate salinity at 15 m
    if min(z) <= 15 && max(z) >= 15

        S15 = ...
            interp1( ...
            z, ...
            ss, ...
            15, ...
            'linear');

        dS15(p) = ...
            S15 - S_surface;

    end
end


fprintf('\nNear-surface salinity differences:\n');

for p = 1:6

    fprintf('%-24s DeltaS = %5.2f psu\n', ...
        panel_names{p}, ...
        dS15(p));

end


%% ========================================================================
%% 7. COMMON AXIS RANGE
%% ========================================================================

all_sal = ...
    PLOT_PROFILE(isfinite(PLOT_PROFILE));

xmin = ...
    floor((min(all_sal)-0.30)*2)/2;

xmax = ...
    ceil((max(all_sal)+0.30)*2)/2;

xmin = floor(xmin);
xmax = ceil(xmax);

fprintf('\nCommon salinity axis: %.1f to %.1f psu\n', ...
    xmin,xmax);


%% ========================================================================
%% 8. CREATE FIGURE
%% ========================================================================

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 15.0 15]);

tl = tiledlayout( ...
    fig, ...
    3,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% ========================================================================
%% 9. DRAW SIX PANELS
%% ========================================================================

for p = 1:6

    ax = nexttile(tl,p);

    hold(ax,'on');


    %% --------------------------------------------------------------------
    % 0–15 m shaded layer
    %% --------------------------------------------------------------------

    patch( ...
        ax, ...
        [xmin xmax xmax xmin], ...
        [0 0 15 15], ...
        SHADE_COLOR, ...
        'EdgeColor','none', ...
        'FaceAlpha',1);


    %% --------------------------------------------------------------------
    % 15 m reference line
    %% --------------------------------------------------------------------

    yline( ...
        ax, ...
        15, ...
        '--', ...
        'Color',REF_COLOR, ...
        'LineWidth',0.75);


    %% --------------------------------------------------------------------
    % Vertical salinity profile
    %% --------------------------------------------------------------------

    prof = ...
        squeeze(PLOT_PROFILE(p,:));

    plot( ...
        ax, ...
        prof, ...
        dep_plot, ...
        '-', ...
        'Color',LINE_COLOR, ...
        'LineWidth',LINE_WIDTH);


    %% --------------------------------------------------------------------
    % Axis range
    %% --------------------------------------------------------------------

    set(ax,'YDir','reverse');

    xlim(ax,[xmin xmax]);
    ylim(ax,[0 MAX_PLOT_DEPTH]);

    yticks(ax,0:10:MAX_PLOT_DEPTH);


    %% --------------------------------------------------------------------
    % X axis
    %
    % Upper four panels:
    %   ticks retained but labels hidden
    %
    % Bottom two panels:
    %   labels shown every 2 psu, horizontal
    %% --------------------------------------------------------------------

    if p <= 4

        xticks(ax,xmin:2:xmax);

        ax.XTickLabel = [];

        xlabel(ax,'');

    else

        xticks(ax,xmin:2:xmax);

        ax.XTickLabelRotation = 0;

        xlabel( ...
            ax, ...
            'Salinity (psu)', ...
            'FontName',FONT_NAME, ...
            'FontSize',FONT_LABEL, ...
            'FontWeight','normal');

    end


    %% --------------------------------------------------------------------
    % Y axis labels only on left column
    %% --------------------------------------------------------------------

    if mod(p,2)==1

        ylabel( ...
            ax, ...
            'Depth (m)', ...
            'FontName',FONT_NAME, ...
            'FontSize',FONT_LABEL, ...
            'FontWeight','normal');

    else

        ylabel(ax,'');

    end


    %% --------------------------------------------------------------------
    % Grid
    %% --------------------------------------------------------------------

    grid(ax,'on');

    ax.GridColor = GRID_COLOR;
    ax.GridAlpha = 0.38;
    ax.MinorGridAlpha = 0;


    %% --------------------------------------------------------------------
    % Panel title
    %% --------------------------------------------------------------------

    title( ...
        ax, ...
        [panel_letters{p} ' ' panel_names{p}], ...
        'FontName',FONT_NAME, ...
        'FontSize',FONT_TITLE, ...
        'FontWeight','bold');


    %% --------------------------------------------------------------------
    % 0–15 m annotation
    %% --------------------------------------------------------------------

    text( ...
        ax, ...
        xmin + 0.045*(xmax-xmin), ...
        5.5, ...
        '0–15 m', ...
        'FontName',FONT_NAME, ...
        'FontSize',FONT_TEXT, ...
        'Color',[0.12 0.12 0.12]);


    %% --------------------------------------------------------------------
    % Delta S annotation
    %% --------------------------------------------------------------------

    if isfinite(dS15(p))

        text( ...
            ax, ...
            xmin + 0.25*(xmax-xmin), ...
            51, ...
            sprintf('\\DeltaS = %.2f psu',dS15(p)), ...
            'FontName',FONT_NAME, ...
            'FontSize',FONT_TEXT, ...
            'Color',[0.15 0.15 0.15]);

    end


    %% --------------------------------------------------------------------
    % Final axis style
    %% --------------------------------------------------------------------

    set(ax, ...
        'FontName',FONT_NAME, ...
        'FontSize',FONT_AXIS, ...
        'LineWidth',0.75, ...
        'TickDir','out', ...
        'TickLength',[0.012 0.012], ...
        'Box','on', ...
        'Layer','top');

    ax.XColor = [0.12 0.12 0.12];
    ax.YColor = [0.12 0.12 0.12];

    % Prevent MATLAB from automatically rotating labels
    ax.XTickLabelRotation = 0;

end


%% ========================================================================
%% 10. EXPORT PNG
%% ========================================================================

drawnow;

set(fig,'InvertHardcopy','off');

exportgraphics( ...
    gcf, ...
    fullfile(cfg.outputDir,'Fig02_VerticalSalinityProfiles.png'), ...
    'Resolution',600);


fprintf('\n==============================================\n');
fprintf('Figure completed successfully.\n');
fprintf('Saved PNG:\n%s\n',out_png);
fprintf('==============================================\n');
