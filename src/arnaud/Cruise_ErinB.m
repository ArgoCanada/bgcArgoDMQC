clearvars

% ID Lat    Lon    Depth
Stations = ...
[0	47.5667	-52.7	0
1	52.9726	-51.4684	1972
2	53.2278	-51.1042	2837
3	56.5	-52.3	3522
4	58.8143	-54.3218	3363
5	58.7829	-51.2256	3501
6	58.7932	-47.6507	3077
7	57.7082	-46.3546	2776
8	57.6965	-49.2829	3529
9	57.7082	-52.364	3500
10	57.7397	-55.1331	3231
11	56.5095	-55.0504	3019
12	56.5	-52.3	3522
13	56.5248	-50.1841	3605
14	56.5897	-48.2988	3641
15	56.5573	-46.06	3504
16	55.3878	-50.3756	3466
17	55.3702	-52.4229	3196
18	56.5	-52.3	3522
19	47.5667	-52.7		0];

Station_names = {'St. John''s','OSNAP','OSNAP','SC-1','','','','','','','','','SC-2','','','','','','SC-3','St. John''s'};

disp('get data')
lakes = shaperead('/misc/10/output/laurenta/data/NaturalEarth/ne_10m_lakes/ne_10m_lakes.shp', 'UseGeoCoords', true);
rivers = shaperead('/misc/10/output/laurenta/data/NaturalEarth/ne_10m_rivers_lake_centerlines/ne_10m_rivers_lake_centerlines.shp', 'UseGeoCoords', true);
ice = shaperead('/misc/10/output/laurenta/data/NaturalEarth/ne_50m_glaciated_areas/ne_50m_glaciated_areas.shp', 'UseGeoCoords', true);

% ETOPO1 data available here: https://maps.ngdc.noaa.gov/viewers/wcs-client/
ncload('/misc/10/output/laurenta/data/ETOPO/etopo1_m120_10_30_80.nc')
addpath('/misc/10/output/laurenta/data/cptcmap')

Z = Band1;
R = georasterref('RasterSize', size(Z), 'Latlim', [min(lat) max(lat)], 'Lonlim', [min(lon) max(lon)]);

disp('make figure')
pos = [1 1500 1500 1500];
axpos = [0.05 0.05 0.9 0.9];
fs = 20;

figure('Position',pos,'Visible','off')

ax1 = axesm('MapProjection','lambert','MapLatLimit',[40 70],'MapLonLimit',[-70 -30]);
set(ax1,'Position',axpos,'XColor','none','YColor','none')
framem;
gridm('on');
geoshow(Z,R, 'DisplayType', 'texturemap');
cptcmap('GMT_globe', 'mapping', 'direct');
caxis([-10000 10000])

hold on
g2 = geoshow(lakes,'FaceColor',[0.9451    0.9882    1.0000],'EdgeColor','None','LineWidth',0.5);
g3 = geoshow(rivers,'Color',[0.9451    0.9882    1.0000],'LineWidth',0.5);
g4 = geoshow(ice,'FaceColor',[1 1 1],'EdgeColor','None');

plotm(Stations(1:end-1,2),Stations(1:end-1,3),'k-','LineWidth',1)
plotm(Stations([5:12,14:18],2),Stations([5:12,14:18],3),'ko','MarkerFaceColor','k','MarkerSize',11)
plotm(Stations(2:3,2),Stations(2:3,3),'ks','MarkerFaceColor','k','MarkerSize',13)
%textm(Stations([2,3,5:12,14:18],2),Stations([2,3,5:12,14:18],3),num2cell(Stations([2,3,5:12,14:18],1)),'HorizontalAlignment','Left','VerticalAlignment','Top')
plotm(Stations(4,2),Stations(4,3),'rs','MarkerFaceColor','r','MarkerSize',25)
plotm(Stations(1,2),Stations(1,3),'ro','MarkerFaceColor','r','MarkerSize',11)

tightmap
setm(ax1,'glinestyle','-','glinewidth',0.25,'gcolor',[0.1500 0.1500 0.1500],'MLineLocation',10,'PLineLocation',10,...
         'ParallelLabel','on','MeridianLabel','on','MLabelLocation',-70:10:-30,'PLabelLocation',40:10:70,...
         'LabelRotation','on','FontSize',fs,'FontWeight','bold','FontName','helvetica')

disp('save figure')
plotname = 'Cruise_map_ErinB';
set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
print(gcf,'-dpng','-r100',plotname)
close

