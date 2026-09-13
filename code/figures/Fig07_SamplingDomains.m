%% Repository-relative setup
thisDir = fileparts(mfilename('fullpath'));
codeDir = fileparts(thisDir);
addpath(codeDir,fullfile(codeDir,'functions'));
cfg = project_config();
sys = river_shelf_systems();
close all;

%% ==============================================================
%  Figure 7. Sampling domains and sensitivity radii
%  Projection: polar stereographic (90°N)
% ==============================================================

Re_km = 6371;
% ---------------- River mouth coordinates ----------------
mouth_lon = sys.mouthLon;
mouth_lat = sys.mouthLat;

% ---------------- Figure & Projection ----------------
figure('Units','centimeters','Position',[2 2 13 13],'Color','w');
subplot('Position',[0.01 0.02 0.98 0.95])
m_proj('stereographic','lat',90,'long',0,'radius',32);
hold on;

% ---------------- Background (light & clean) ----------------
[~, h] = m_tbase('contour',[-6000:1000:0]);
set(h,'LineWidth',0.3,'LineColor',[0.85 0.88 0.9]);

% m_coast('patch',[0.92 0.9 0.85],'edgecolor',[0.4 0.4 0.4],'LineWidth',0.5);
m_gshhs('fc','color',[.5 .5 .5]);
m_grid('box','fancy','tickdir','in','gridlines','no', ...
       'fontsize',9,'yticklabel',[],'xticklabel',[]);
%% ---------------- Sampling circles & robustness setup ----------------
az = linspace(0,2*pi,300);
mainColor = [0.85 0.1 0.1];
sensColor = [0.3 0.3 0.3];
arrowColor = [0 0 1];
% rad_main1=[150;150;250;250;400;150];

for r = 1:numel(mouth_lon)
    lon0 = mouth_lon(r); lat0 = mouth_lat(r);
    lat0r = deg2rad(lat0); lon0r = deg2rad(lon0);
    az_r = az;

    rad_main = 450;
    ang = rad_main / Re_km;
    latc = asin(sin(lat0r).*cos(ang) + cos(lat0r).*sin(ang).*cos(az_r));
    lonc = lon0r + atan2(sin(az_r).*sin(ang).*cos(lat0r), cos(ang)-sin(lat0r).*sin(latc));
    latc = rad2deg(latc); lonc = rad2deg(lonc);
    lonc = mod(lonc+180,360)-180;
    m_plot(lonc,latc,'Color',mainColor,'LineWidth',2);

    for rr = [150 300 450 600]
        ang = rr / Re_km;
        latc = asin(sin(lat0r).*cos(ang) + cos(lat0r).*sin(ang).*cos(az_r));
        lonc = lon0r + atan2(sin(az_r).*sin(ang).*cos(lat0r), cos(ang)-sin(lat0r).*sin(latc));
        latc = rad2deg(latc); lonc = rad2deg(lonc);
        lonc = mod(lonc+180,360)-180;
        m_plot(lonc,latc,'LineStyle',':','Color',sensColor,'LineWidth',0.8);
    end

    shift_deg = 1.5;
    lat_shift = lat0 + shift_deg;
    lon_shift = lon0;
    [x0,y0] = m_ll2xy(lon0,lat0);
    [x1,y1] = m_ll2xy(lon_shift,lat_shift);
    quiver(x0,y0,(x1-x0)*0.7,(y1-y0)*0.7,0,'Color',arrowColor,'LineWidth',0.7,'MaxHeadSize',1);

end
m_coast('patch',[0.92 0.9 0.85],'edgecolor',[0.4 0.4 0.4],'LineWidth',0.5);


color = [51 0 0]/256; size1 = 1.5;
% ---------------- River basin outlines (light tone) ----------------
colorBasin = [0.5 0.2 0.2];
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

% ---------------- Marginal Seas boundaries (faint blue) ----------------
S = shaperead(fullfile(cfg.geoDir,'arctic_marginal_seas.shp'));
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
lat=[70	 75	75.8060985968788	79.3823483107168	76.7374330759235];
name={'Chukchi Sea','Beaufort Sea','East Siberian Sea','Kara Sea','Laptev Sea'};
for i = 1:numel(name)

    m_text(lon(i), lat(i), name{i}, ...
        'FontSize',7,'FontName','Times New Roman', ...
        'FontWeight','bold','color',[0.25 0.25 0.25], ...
        'HorizontalAlignment','center','FontAngle','italic');
end
badgeSize = 9; badgeBlue = [50 90 160]/255; badgeEdge = [30 60 110]/255;
numFont = 8;
pos = [-120.5 60.8; -147.5 66.2; 157.0 65.5; 124.0 63.0; 94.0 61.0; 76.0 62];
labels = {'1','2','3','4','5','6'};
for k = 1:6
    m_line(pos(k,1), pos(k,2), 'marker','o','markersize',badgeSize, ...
           'MarkerFaceColor',badgeBlue,'color',badgeEdge,'LineWidth',0.8);
    m_text(pos(k,1), pos(k,2), labels{k}, 'color','w','FontSize',numFont, ...
           'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
end
labels = {'1','2','3','4','5','6'};numFont = 8;
for r = 1:6
    lon0 = mouth_lon(r); lat0 = mouth_lat(r);
    m_line(lon0, lat0, 'marker','^','markersize',8, ...
           'MarkerFaceColor',[204 0 0]/256,'color',[204 0 0]/256,'LineWidth',0.8);
    m_text(lon0, lat0, labels{r}, 'color','w','FontSize',numFont-1, ...
           'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
end

%%

set(gcf,'InvertHardcopy','off');
exportgraphics(gcf,fullfile(cfg.outputDir,'Fig07_SamplingDomains.png'),'Resolution',600);
% exportgraphics(gcf,fullfile(cfg.outputDir,'Fig07_SamplingDomains.tif'),'Resolution',600);
