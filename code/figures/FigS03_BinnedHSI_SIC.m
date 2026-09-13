%% Repository-relative setup
thisDir = fileparts(mfilename('fullpath'));
codeDir = fileparts(thisDir);
addpath(codeDir,fullfile(codeDir,'functions'));
cfg = project_config();
sys = river_shelf_systems();
close all;

%% ========================================================================
% Figure S3
% Binned relationships between standardized JAS HSI
% and JAS sea-ice concentration
%
% Period: 1979-2024
%
% FINAL VERSION:
%   - HSI and SIC both taken from annual_river_shelf_metrics_1979_2024.mat
%   - x axis: detrended + standardized JAS HSI
%   - y axis: JAS SIC (%)
%   - Bin width: 0.5 standard deviation
%   - Error bars: 95% CI of bin-mean SIC
%% ========================================================================




%% ========================================================================
% 1. LOAD FINAL DATA
%% ========================================================================

dataFile = fullfile(cfg.dataDir,'annual_river_shelf_metrics_1979_2024.mat');

D = load(dataFile);


%% ------------------------------------------------------------------------
% Years
%% ------------------------------------------------------------------------

if isfield(D,'years')
    years = double(D.years(:));
else
    years = (1979:2024)';
end


%% ------------------------------------------------------------------------
% Standardized HSI
%
% Prefer HSI_JAS_z directly.
% If not available, standardize HSI_JAS_dt.
%% ------------------------------------------------------------------------

if isfield(D,'HSI_JAS_z')

    HSI_z = double(D.HSI_JAS_z);

elseif isfield(D,'HSI_JAS_dt')

    HSI_dt = double(D.HSI_JAS_dt);

    if size(HSI_dt,1) == 6
        HSI_dt = HSI_dt';
    end

    HSI_z = nan(size(HSI_dt));

    for ib = 1:6
        HSI_z(:,ib) = local_zscore(HSI_dt(:,ib));
    end

else

    error('Annual metrics file does not contain HSI_JAS_z or HSI_JAS_dt.');

end


%% ------------------------------------------------------------------------
% Raw JAS SIC
%% ------------------------------------------------------------------------

if ~isfield(D,'SIC_JAS')
    error('Annual metrics file does not contain SIC_JAS.');
end

SIC = double(D.SIC_JAS);


%% ------------------------------------------------------------------------
% Convert orientation to 46 × 6
%% ------------------------------------------------------------------------

if size(HSI_z,1) == 6 && size(HSI_z,2) == 46
    HSI_z = HSI_z';
end

if size(SIC,1) == 6 && size(SIC,2) == 46
    SIC = SIC';
end


assert(isequal(size(HSI_z),[46 6]), ...
    'HSI_JAS_z must be 46 × 6.');

assert(isequal(size(SIC),[46 6]), ...
    'SIC_JAS must be 46 × 6.');


%% ========================================================================
% 2. CONVERT SIC TO %
%% ========================================================================

sicValid = SIC(isfinite(SIC));

if isempty(sicValid)
    error('SIC_JAS is empty.');
end


sicMax = max(sicValid);


if sicMax <= 1.5

    % SIC stored as 0-1
    SIC_pct = SIC * 100;

elseif sicMax <= 100.5

    % Already percentage
    SIC_pct = SIC;

else

    error('Unexpected SIC range. Check SIC_JAS units.');

end


fprintf('\nSIC range used for Figure S3:\n');
fprintf('%.2f - %.2f %%\n', ...
    min(SIC_pct(:),[],'omitnan'), ...
    max(SIC_pct(:),[],'omitnan'));


%% ========================================================================
% 3. SYSTEM NAMES
%% ========================================================================

riverNames = { ...
    'Mackenzie-Beaufort', ...
    'Yukon-Chukchi', ...
    'Kolyma-East Siberian', ...
    'Lena-Laptev', ...
    'Yenisei-Kara', ...
    'Ob-Kara'};


nBasin = 6;


%% ========================================================================
% 4. DEFINE COMMON HSI BINS
%
% Bin width = 0.5 standard deviation
%
% Automatically extend the range so that ALL observations are retained.
%% ========================================================================

binWidth = 0.5;


allX = ...
    HSI_z(isfinite(HSI_z));


xMin = ...
    floor(min(allX)/binWidth)*binWidth;


xMax = ...
    ceil(max(allX)/binWidth)*binWidth;


% At least retain the familiar -2 to +2 range
xMin = min(xMin,-2);
xMax = max(xMax, 2);


edges = ...
    xMin:binWidth:xMax;


% Safety in case floating point misses upper endpoint
if edges(end) < xMax
    edges(end+1) = xMax;
end


nbins = ...
    numel(edges)-1;


binMid = ...
    (edges(1:end-1)+edges(2:end))/2;


fprintf('\nHSI bin range:\n');
fprintf('%.1f to %.1f, bin width = %.1f\n', ...
    xMin,xMax,binWidth);


%% ========================================================================
% 5. PREALLOCATE BIN STATISTICS
%% ========================================================================

binMean = ...
    nan(nbins,nBasin);


binCILow = ...
    nan(nbins,nBasin);


binCIHigh = ...
    nan(nbins,nBasin);


binN = ...
    zeros(nbins,nBasin);


%% ========================================================================
% 6. CALCULATE BIN STATISTICS
%% ========================================================================

MIN_BIN_N = 3;


for ib = 1:nBasin

    x = HSI_z(:,ib);
    y = SIC_pct(:,ib);


    valid = ...
        isfinite(x) & ...
        isfinite(y);


    x = x(valid);
    y = y(valid);


    % ------------------------------------------------------------
    % Discretize standardized HSI
    % ------------------------------------------------------------

    groups = ...
        discretize( ...
        x, ...
        edges);


    % ------------------------------------------------------------
    % Bin statistics
    % ------------------------------------------------------------

    for k = 1:nbins

        yy = ...
            y(groups == k);


        yy = ...
            yy(isfinite(yy));


        binN(k,ib) = ...
            numel(yy);


        if binN(k,ib) >= MIN_BIN_N

            binMean(k,ib) = ...
                mean(yy);


            % ----------------------------------------------------
            % 95% confidence interval of the bin mean
            % ----------------------------------------------------

            se = ...
                std(yy,0) / ...
                sqrt(binN(k,ib));


            tcrit = ...
                tinv( ...
                0.975, ...
                binN(k,ib)-1);


            binCILow(k,ib) = ...
                binMean(k,ib) - ...
                tcrit*se;


            binCIHigh(k,ib) = ...
                binMean(k,ib) + ...
                tcrit*se;

        end

    end

end


%% ========================================================================
% 7. PRINT BIN RESULTS
%% ========================================================================

fprintf('\n============================================================\n');
fprintf('FIGURE S3 BIN STATISTICS\n');
fprintf('============================================================\n');


for ib = 1:nBasin

    fprintf('\n%s\n',riverNames{ib});

    fprintf('Bin midpoint   N     Mean SIC\n');

    for k = 1:nbins

        if binN(k,ib) > 0

            fprintf('%8.2f    %3d    %8.2f\n', ...
                binMid(k), ...
                binN(k,ib), ...
                binMean(k,ib));

        end

    end

end


%% ========================================================================
% 8. SAVE RESULTS
%% ========================================================================

outDir = cfg.outputDir;


if ~exist(outDir,'dir')
    mkdir(outDir);
end


save( ...
    fullfile( ...
    outDir, ...
    'FigureS3_binned_HSI_SIC_results.mat'), ...
    ...
    'HSI_z', ...
    'SIC_pct', ...
    'years', ...
    'riverNames', ...
    'edges', ...
    'binMid', ...
    'binWidth', ...
    'binMean', ...
    'binCILow', ...
    'binCIHigh', ...
    'binN', ...
    'MIN_BIN_N');


%% ========================================================================
% 9. COMMON Y-AXIS RANGE
%% ========================================================================

allY = ...
    SIC_pct(isfinite(SIC_pct));


yMin = ...
    floor(min(allY)/5)*5 - 5;


yMax = ...
    ceil(max(allY)/5)*5 + 5;


% Keep a sensible range
yMin = max(0,yMin);
yMax = min(100,yMax);


fprintf('\nPlotting SIC range: %.0f - %.0f %%\n', ...
    yMin,yMax);


%% ========================================================================
% 10. FIGURE SETTINGS
%% ========================================================================

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 17.5 17.5]);


tl = tiledlayout( ...
    3,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


fontAxis  = 7.5;
fontLabel = 8.5;
fontTitle = 9.2;


markerSize = 16;


cPointEdge = ...
    [0.12 0.12 0.12];


cPointFace = ...
    [0.72 0.72 0.72];


cBin = ...
    [0.00 0.45 0.85];


cGrid = ...
    [0.75 0.75 0.75];


lineWidthBin = ...
    1.6;


%% ========================================================================
% 11. PLOT SIX SYSTEMS
%% ========================================================================

for ib = 1:nBasin

    ax = ...
        nexttile(tl,ib);


    hold(ax,'on');


    %% --------------------------------------------------------
    % Annual observations
    %% --------------------------------------------------------

    x = ...
        HSI_z(:,ib);


    y = ...
        SIC_pct(:,ib);


    valid = ...
        isfinite(x) & ...
        isfinite(y);


    xObs = ...
        x(valid);


    yObs = ...
        y(valid);


    hObs = ...
        scatter( ...
        ax, ...
        xObs, ...
        yObs, ...
        markerSize, ...
        'o', ...
        'MarkerEdgeColor',cPointEdge, ...
        'MarkerFaceColor',cPointFace, ...
        'LineWidth',0.55);


    try

        hObs.MarkerFaceAlpha = ...
            0.42;

        hObs.MarkerEdgeAlpha = ...
            0.90;

    catch

    end


    %% --------------------------------------------------------
    % Bin-mean SIC and 95% CI
    %% --------------------------------------------------------

    hasBin = ...
        isfinite(binMean(:,ib)) & ...
        isfinite(binCILow(:,ib)) & ...
        isfinite(binCIHigh(:,ib));


    xb = ...
        binMid(hasBin);


    yb = ...
        binMean(hasBin,ib);


    errLow = ...
        yb - ...
        binCILow(hasBin,ib);


    errHigh = ...
        binCIHigh(hasBin,ib) - ...
        yb;


    hBin = ...
        errorbar( ...
        ax, ...
        xb, ...
        yb, ...
        errLow, ...
        errHigh, ...
        '-o', ...
        'Color',cBin, ...
        'LineWidth',lineWidthBin, ...
        'MarkerFaceColor',cBin, ...
        'MarkerEdgeColor',cBin, ...
        'MarkerSize',4.5, ...
        'CapSize',5);


    %% --------------------------------------------------------
    % Axes
    %% --------------------------------------------------------

    xlim( ...
        ax, ...
        [xMin xMax]);


    ylim( ...
        ax, ...
        [yMin yMax]);


    xticks( ...
        ax, ...
        ceil(xMin):1:floor(xMax));


    % use 10% SIC ticks
    yTickStart = ...
        ceil(yMin/10)*10;


    yTickEnd = ...
        floor(yMax/10)*10;


    yticks( ...
        ax, ...
        yTickStart:10:yTickEnd);


    grid(ax,'on');


    ax.GridColor = ...
        cGrid;


    ax.GridAlpha = ...
        0.30;


    set( ...
        ax, ...
        'FontName','Times New Roman', ...
        'FontSize',fontAxis, ...
        'FontWeight','bold', ...
        'LineWidth',0.8, ...
        'TickDir','out', ...
        'Box','off', ...
        'Layer','top');


    %% --------------------------------------------------------
    % X labels only bottom row
    %% --------------------------------------------------------

    if ib == 5 || ib == 6

        xlabel( ...
            ax, ...
            'Standardized JAS HSI', ...
            'FontName','Times New Roman', ...
            'FontSize',fontLabel, ...
            'FontWeight','bold');

    else

        ax.XTickLabel = [];

    end


    %% --------------------------------------------------------
    % Y labels only left column
    %% --------------------------------------------------------

    if ismember(ib,[1 3 5])

        ylabel( ...
            ax, ...
            'JAS SIC (%)', ...
            'FontName','Times New Roman', ...
            'FontSize',fontLabel, ...
            'FontWeight','bold');

    end


    %% --------------------------------------------------------
    % System title
    %% --------------------------------------------------------

    text( ...
        ax, ...
        0.03, ...
        0.94, ...
        riverNames{ib}, ...
        'Units','normalized', ...
        'FontName','Times New Roman', ...
        'FontSize',fontTitle, ...
        'FontWeight','bold', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top');


    %% --------------------------------------------------------
    % Legend
    %% --------------------------------------------------------

    if ib == 1

        lgd = legend( ...
            ax, ...
            [hBin hObs], ...
            { ...
            'Bin-mean SIC', ...
            'Annual observations'}, ...
            'Location','southeast', ...
            'Box','off', ...
            'FontName','Times New Roman', ...
            'FontSize',fontAxis, ...
            'FontWeight','bold');


        lgd.ItemTokenSize = ...
            [10 5];

    end

end


%% ========================================================================
% 12. EXPORT
%% ========================================================================

drawnow;

set(fig,'InvertHardcopy','off');

exportgraphics( ...
    gcf, ...
    fullfile(cfg.outputDir,'FigS03_BinnedHSI_SIC.png'), ...
    'Resolution',600);



%% ========================================================================
% LOCAL FUNCTION
%% ========================================================================

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
        std(x(valid),0);


    if sigma == 0
        return
    end


    z(valid) = ...
        (x(valid)-mu) / ...
        sigma;

end
