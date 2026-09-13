%% Repository-relative setup
thisDir = fileparts(mfilename('fullpath'));
codeDir = fileparts(thisDir);
addpath(codeDir,fullfile(codeDir,'functions'));
cfg = project_config();
sys = river_shelf_systems();
close all;

mouth_lon = sys.mouthLon;
mouth_lat = sys.mouthLat;

% ---------------- Projection & Figure Setup ----------------
figure('Units','centimeters','Position',[2 2 12 12],'Color','w');

m_proj('stereographic','lat',90,'long',0,'radius',40);

hold on
[CS, h] = m_tbase('contour',[-6000:800:0]);
set(h,'LineWidth',0.5,'LineColor',[0.7 0.8 0.9]);

m_coast('patch',[0.92 0.9 0.8],'edgecolor','k','LineWidth',0.8);
m_gshhs('fc','color',[.5 .5 .5]);
m_grid('box','fancy','tickdir','in','gridlines','no','fontsize',12, ...
        'yticklabel',[],'xticklabel',[]);

[CS,h1] = m_etopo2('contour',[-1000 -1000],'color',[0.3 0.3 0.3],'LineWidth',1.5);


% ---------------- River & Basin Outlines ----------------
color = [51 0 0]/256; size1 = 1.5;

B = load_basin_boundary(fullfile(cfg.geoDir,'basin_mackenzie.mat'));
    m_plot(B(1,:),B(2,:),'color',color,'LineWidth',size1);

B = load_basin_boundary(fullfile(cfg.geoDir,'basin_yukon.mat'));
    m_plot(B(1,:),B(2,:),'color',color,'LineWidth',size1);

B = load_basin_boundary(fullfile(cfg.geoDir,'basin_kolyma.mat'));
    m_plot(B(1,:),B(2,:),'color',color,'LineWidth',size1);

B = load_basin_boundary(fullfile(cfg.geoDir,'basin_lena.mat'));
    m_plot(B(1,:),B(2,:),'color',color,'LineWidth',size1);

B = load_basin_boundary(fullfile(cfg.geoDir,'basin_yenisei.mat'));
    m_plot(B(1,:),B(2,:),'color',color,'LineWidth',size1);

B = load_basin_boundary(fullfile(cfg.geoDir,'basin_ob.mat'));
    m_plot(B(1,:),B(2,:),'color',color,'LineWidth',size1);


S = shaperead(fullfile(cfg.geoDir,'arctic_marginal_seas.shp'));
% for i = 1:length(S)
%     if isfield(S(i),'X') && ~isempty(S(i).X)
%         m_line(S(i).X, S(i).Y, 'color',edgeColor, 'LineWidth', 1.1);
%     end
% end

shelfHaloColor = [0.95 0.97 1.00];   % very light blue-white halo
shelfEdgeColor = [0.00 0.25 0.50];   % dark blue boundary

for i = 1:length(S)
    if isfield(S(i), 'X') && ~isempty(S(i).X)
        x = S(i).X;
        y = S(i).Y;

        % halo line
        m_line(x, y, ...
            'Color', shelfHaloColor, ...
            'LineWidth', 2.4);

        % main boundary line
        m_line(x, y, ...
            'Color', shelfEdgeColor, ...
            'LineWidth', 1.15);
    end
end

lon=[-169	-138.652971157731	165.759154582323	80.0845283034114	121.387763108681];
lat=[70	 77	75.8060985968788	79.3823483107168	78.7374330759235];
name={'Chukchi Sea','Beaufort Sea','East Siberian Sea','Kara Sea','Laptev Sea'};
hold on
for i = 1:numel(name)

    m_text(lon(i), lat(i), name{i}, ...
        'FontSize',7,'FontName','Times New Roman', ...
        'FontWeight','bold','color',[0.25 0.25 0.25], ...
        'HorizontalAlignment','center','FontAngle','italic');
end


labels = {'1','2','3','4','5','6'};
badgeSize = 9; badgeBlue = [50 90 160]/255; badgeEdge = [30 60 110]/255;
numFont = 8;
pos = [-120.5 60.8; -147.5 66.2; 157.0 65.5; 124.0 63.0; 94.0 61.0; 76.0 58.5];
for r = 1:numel(mouth_lon)
    hold on
    lon0 = mouth_lon(r); lat0 = mouth_lat(r);

    m_line(lon0, lat0, 'marker','^','markersize',8, ...
           'MarkerFaceColor',[204 0 0]/256,'color',[204 0 0]/256,'LineWidth',0.8);
    m_text(lon0, lat0, labels{r}, 'color','w','FontSize',numFont-1, ...
           'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
end
for k = 1:6
    m_line(pos(k,1), pos(k,2), 'marker','o','markersize',badgeSize, ...
           'MarkerFaceColor',badgeBlue,'color',badgeEdge,'LineWidth',0.8);
    m_text(pos(k,1), pos(k,2), labels{k}, 'color','w','FontSize',numFont, ...
           'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
end

%% ============================================================

set(gcf,'InvertHardcopy','off');
exportgraphics(gcf,fullfile(cfg.outputDir,'Fig01_StudyArea.png'),'Resolution',600);
exportgraphics(gcf,fullfile(cfg.outputDir,'Fig01_StudyArea.tif'),'Resolution',600);
