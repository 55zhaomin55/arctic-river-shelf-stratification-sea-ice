%% Repository-relative setup
thisDir = fileparts(mfilename('fullpath'));
codeDir = fileparts(thisDir);
addpath(codeDir,fullfile(codeDir,'functions'));
cfg = project_config();
sys = river_shelf_systems();
close all;

%% ================================================================
% Figure 8
% Correlations of JAS HSI with May-September cumulative runoff
% and JAS sea-ice concentration
% 1979-2024
%% ================================================================




%% ================================================================
% 1. Load saved results
%% ================================================================

dataFile = fullfile(cfg.dataDir,'annual_river_shelf_metrics_1979_2024.mat');

load(dataFile, ...
    'systems', ...
    'R_HSI_Q', ...
    'P_HSI_Q', ...
    'R_HSI_SIC', ...
    'P_HSI_SIC');


%% ================================================================
% 2. Basic check
%% ================================================================

R_HSI_Q   = R_HSI_Q(:);
P_HSI_Q   = P_HSI_Q(:);

R_HSI_SIC = R_HSI_SIC(:);
P_HSI_SIC = P_HSI_SIC(:);

nSystem = length(systems);

if nSystem ~= 6
    error('Expected six river-shelf systems.');
end

if length(R_HSI_Q) ~= 6 || length(R_HSI_SIC) ~= 6
    error('Correlation results must contain six values.');
end


%% ================================================================
% 3. Display correlation results
%% ================================================================

ResultTable = table( ...
    systems(:), ...
    R_HSI_Q, ...
    P_HSI_Q, ...
    R_HSI_SIC, ...
    P_HSI_SIC, ...
    'VariableNames', { ...
    'System', ...
    'R_HSI_Q', ...
    'P_HSI_Q', ...
    'R_HSI_SIC', ...
    'P_HSI_SIC'});

disp(ResultTable);


%% ================================================================
% 4. Prepare Figure 8 data
%% ================================================================

Rplot = [R_HSI_Q, R_HSI_SIC];
Pplot = [P_HSI_Q, P_HSI_SIC];

x = 1:6;


%% ================================================================
% 5. Figure settings
%% ================================================================

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[3 3 18 10]);

ax = axes(fig);

hold(ax,'on');


%% ================================================================
% 6. Grouped bar plot
%% ================================================================

bh = bar( ...
    ax, ...
    x, ...
    Rplot, ...
    'grouped', ...
    'BarWidth',0.68);


% Colors
blue   = [0.18 0.47 0.82];
orange = [0.94 0.56 0.20];


bh(1).FaceColor = blue;
bh(1).EdgeColor = 'none';

bh(2).FaceColor = orange;
bh(2).EdgeColor = 'none';


%% ================================================================
% 7. Zero line
%% ================================================================

yline( ...
    ax, ...
    0, ...
    '-', ...
    'Color',[0.30 0.30 0.30], ...
    'LineWidth',0.8);

%% ================================================================
% 8. Significance stars
%
% * = p < 0.05
%% ================================================================

drawnow;


for j = 1:2

    xpos = bh(j).XEndPoints;

    if j == 1
        starColor = blue;
    else
        starColor = orange;
    end


    for i = 1:nSystem

        if Pplot(i,j) < 0.05

            yvalue = Rplot(i,j);


            if yvalue >= 0

                ystar = yvalue + 0.035;

                valign = 'bottom';

            else

                ystar = yvalue - 0.035;

                valign = 'top';

            end


            text( ...
                ax, ...
                xpos(i), ...
                ystar, ...
                '*', ...
                ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment',valign, ...
                ...
                'FontName','Times New Roman', ...
                'FontSize',13, ...
                'FontWeight','bold', ...
                'Color',starColor);

        end

    end

end


%% ================================================================
% 9. Axis labels
%% ================================================================

xlim(ax,[0.45 6.55]);

ylim(ax,[-0.25 1.00]);

yticks(ax,-0.2:0.2:1.0);


xticks(ax,1:6);

xticklabels(ax,{ ...
    'Mackenzie-Beaufort', ...
    'Yukon-Chukchi', ...
    'Kolyma-East Siberian', ...
    'Lena-Laptev', ...
    'Yenisei-Kara', ...
    'Ob-Kara'});

ylabel( ...
    ax, ...
    'Correlation R', ...
    'FontName','Times New Roman', ...
    'FontSize',11, ...
    'FontWeight','normal');


%% ================================================================
% 10. Legend
%% ================================================================

lgd = legend( ...
    ax, ...
    {'HSI-Q','HSI-SIC'}, ...
    'Location','southeast', ...
    'Box','off');

set(lgd, ...
    'FontName','Times New Roman', ...
    'FontSize',9);


%% ================================================================
% 11. General formatting
%% ================================================================

set(ax, ...
    'FontName','Times New Roman', ...
    'FontSize',9, ...
    ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    ...
    'Box','on', ...
    ...
    'XGrid','off', ...
    'YGrid','on', ...
    'GridAlpha',0.18, ...
    ...
    'TickLabelInterpreter','none', ...
    'Layer','top');


% Make x labels slightly smaller
ax.XAxis.FontSize = 8;


%% ================================================================
% 12. Adjust axes position
%% ================================================================

ax.Position = [0.09 0.14 0.88 0.82];


ax.XAxis.FontWeight = 'bold';
ax.XAxis.FontSize   = 7;
%% ================================================================
% 13. Export
%% ================================================================
drawnow;

set(fig,'InvertHardcopy','off');

exportgraphics( ...
    gcf, ...
    fullfile(cfg.outputDir,'Fig08_Runoff_HSI_SIC.png'), ...
    'Resolution',600);


