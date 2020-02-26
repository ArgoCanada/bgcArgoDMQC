clearvars
zones = {'NAtl','AC_LB_BB','SO_50','Global','NWAtl'};% Atlantic Canada, LS+BB; Southern Ocean (>50S); Global case
zones_names = {'North Atlantic','Atlantic Canada','Southern Ocean','Global','Northwest North Atlantic',};
date_start_argo = datenum(1997,1,1);
date_start = datenum(1950,1,1);
date_stop = datenum(2019,9,1);

vtypes = {'DOXY','NITRATE','PH_IN_SITU_TOTAL','CHLA','BBP700','DOWNWELLING_PAR','CDOM'};
vnames = {'OXYGEN','NITRATE','PH','CHLA','BBP700','PAR','CDOM'};
vnames_low = {'do','no3','ph','chl','bbp','par','cdom'};
lnames = {'Oxygen','Nitrate','pH','Chlorophyll','Backscatter','Light','CDOM'};

wod_names.osd = {'Oxygen','Nitrate','pH','Chlorophyll','','',''};
wod_names.ctd = {'Oxygen','','pH','Chlorophyll','','',''};
wod_type = {'osd','ctd'};

disp('Set region')

NA_ll = [32  32      40        41 50  80   80   32;
         -100 -4.6474 -4.6474 -0.5 10  10 -100 -100];

% transect_xy = [
%
%    217    84
%    217   188
%    209   188
%    209   334
%    196   337
%    193   341
%    167   344
%    160   344
%    160   342
%    157   331
%    154   327
%    155   323
%    154   319
%    140   319
%    140   264
%    124   239
%    133   224
%    134   224
%    134   221
%    135   219
%    136   219
%    136   217
%    135   209
%    135   164
%     86   136
%     86    84
%     95    84];
%

WOD = load('WODstats_data','dataset');

ll = [42.00  -40.00
      60.00  -40.00
      60.00  -43.7563
      76.9374  -43.7563
      76.9374  -63.0622
      76.9374  -68.8143
      77.0445  -71.2021
      75.2687  -79.5064
      74.7177  -81.1278
      74.5924  -80.5044
      73.6255  -77.9679
      73.1188  -77.5887
      72.8949  -76.3045
      72.5177  -75.5397
      71.5341  -78.7835
      66.7053  -68.8382
      63.3478  -69.9071
      62.0954  -66.3152
      62.1360  -66.0641
      61.7930  -65.8167
      61.6015  -65.4046
      61.6404  -65.1520
      61.4073  -64.9955
      60.4241  -64.6588
      54.6531  -62.2532
      49.7625  -72.7457
      41.7871  -72.7457
      41.7871  -40.00];
latp = ll(:,1);
lonp = ll(:,2);

ll2 = [30.00  -35.00
      70.00  -35.00
      70  -75.5
      66.7053  -68.8382
      63.3478  -69.9071
      62.0954  -66.3152
      62.1360  -66.0641
      61.7930  -65.8167
      61.6015  -65.4046
      61.6404  -65.1520
      61.4073  -64.9955
      60.4241  -64.6588
      54.6531  -62.2532
      49.7625  -72.7457
      32  -83
      30  -83
      30  -35.00];
latp_large = ll2(:,1);
lonp_large = ll2(:,2);

if 0
  ncload('/misc/10/output/laurenta/data/ETOPO/etopo1_m100_10_20_80.nc')
  [long,latg] = meshgrid(lon,lat);
  idx = long>min(lonp_large)-1&long<max(lonp_large)+1&latg>min(latp_large)-1&latg<max(latp)+1;
  bathy = Band1(idx);
  longy = long(idx);
  latgy = latg(idx);

  for itype = 1:numel(wod_type)
      %wod_names_osd = {'Oxygen','Nitrate','pH','Chlorophyll','',''};
      %wod_names_ctd = {'Oxygen','','pH','','',''};
      disp(wod_type{itype})
      for ivar = 1:numel(wod_names.(wod_type{itype}))
          if ~isempty(wod_names.(wod_type{itype}){ivar})
              disp(wod_names.(wod_type{itype}){ivar})
              % Get indices of all files with user-specified polygon and timeframe:
              lona = WOD.dataset.(wod_type{itype}).(wod_names.(wod_type{itype}){ivar})(:,2);
              lata = WOD.dataset.(wod_type{itype}).(wod_names.(wod_type{itype}){ivar})(:,3);
              iloc = find(inpolygon(lona,lata,lonp,latp));
              idx = zeros(numel(lona),1);
              counter = floor(((numel(iloc)/10):numel(iloc)/10:numel(iloc)));
              for ista = 1:numel(iloc)
                  if ismember(ista,counter)
                    disp([num2str(ista/numel(iloc)*100,'%1.0f'),'%'])
                  end
                  sqrt_dist = sqrt((longy-lona(iloc(ista))).^2+(latgy-lata(iloc(ista))).^2);
                  idx(iloc(ista),1) = double(bathy(find(sqrt_dist==min(sqrt_dist(:))))<-200);
              end
              [in1,on1] = inpolygon(lona,lata,NA_ll(2,:),NA_ll(1,:));
              [in2,on2] = inpolygon(lona.*idx,lata.*idx,lonp,latp);
              [in3,on3] = inpolygon(lona.*idx,lata.*idx,lonp_large,latp_large);
              WOD_inpoly{1}.(wod_type{itype}).(wod_names.(wod_type{itype}){ivar}) = in1+on1;
              WOD_inpoly{2}.(wod_type{itype}).(wod_names.(wod_type{itype}){ivar}) = in2+on2;
              WOD_inpoly{3}.(wod_type{itype}).(wod_names.(wod_type{itype}){ivar}) = lata<-50;
              WOD_inpoly{4}.(wod_type{itype}).(wod_names.(wod_type{itype}){ivar}) = lona*0+1;
              WOD_inpoly{5}.(wod_type{itype}).(wod_names.(wod_type{itype}){ivar}) = in3+on3;
          end
      end
  end
  save WOD_inpoly_CFI WOD_inpoly
else
  load WOD_inpoly_CFI
end
fs = 30;
pos = [1 800 1800 800];

if 1
    for izone = 5 % 5: NW North Atlantic only, 2: Atlantic Canada Only
        fname = ['ArgoCount_',zones{izone},'_',datestr(date_start_argo,'yyyymmdd'),'_',datestr(date_stop,'yyyymmdd')];
        disp(['Load ',fname])
        load(fname)

        tvec = datevec(t);
        tvec_start = datevec(date_start);
        tvec_stop = datevec(date_stop);
        years = tvec_start(1):tvec_stop(1);

        dataset = [];
        itime = 1;
        for iyear = 1:numel(years)-1
            dataset.(zones{izone}).glob.all(itime,1) = sum(count.glob.all(tvec(:,1)==years(iyear)));
            dataset.(zones{izone}).glob.time(itime,1) = datenum(years(iyear),7,1);
            for ivar = 1:numel(vtypes)-2
                dataset.(zones{izone}).bio.(vnames_low{ivar})(itime,1) = sum(count.bio.(vnames_low{ivar})(tvec(:,1)==years(iyear)));
                if ivar==1
                    dataset.(zones{izone}).bio.time(itime,1) = datenum(years(iyear),7,1);
                end
                for itype = 1:numel(wod_type)
                    if ~isempty(wod_names.(wod_type{itype}){ivar})
                        tvec_tmp = datevec(WOD.dataset.(wod_type{itype}).(wod_names.(wod_type{itype}){ivar})(:,1));
                        idx = tvec_tmp(:,1)==years(iyear)&WOD_inpoly{izone}.(wod_type{itype}).(wod_names.(wod_type{itype}){ivar});
                        dataset.(zones{izone}).wod.(wod_type{itype}).(vnames_low{ivar})(itime,1) = sum(idx);
                    end
                end
            end
            itime = itime + 1;
        end

        if 0
            fname = '/misc/10/output/laurenta/data/ETOPO/ETOPO1_Ice_g_gmt4.grd';
            disp(['Load ',fname])
            ncload(fname)
            [X,Y] = meshgrid(x,y);
            R = georasterref('RasterSize', size(z), 'Latlim', [min(y) max(y)], 'Lonlim', [min(x) max(x)]);
            [in1,on1] = inpolygon(X,Y,lonp_large,latp_large);
            NWAg = (in1+on1)&z<-2000;
            sareakm2 = areamat(NWAg, R, wgs84Ellipsoid('kilometers'));
            disp(['NWA area: ',num2str(sareakm2,'%06.0f'),' km2'])
            return
        end
        if 1 % New figure for CFI2020 proposal
            % Numer of stations
            fs1 = 55;
            ymax = [3750,3000,3000,3500,3000,3000,3000];
            time = dataset.(zones{izone}).bio.time;
            timef = [datenum(2021:2024,7,1),datenum(2025:2026,1,1)];

            tvec = datevec(time);
            time = datenum(tvec(:,1),1,1);
            tvecf = datevec(timef);
            idx_0014 = tvec(:,1)>1999&tvec(:,1)<2015;
            nfloats = [0,460,1184,1914,2190,2190];

            float_ndep = cumsum([0,0,4,0,3,0,3,0,3,0,4,0,3,0,3,0,4,0,3,0,0]);
            float_tdep = datenum([2022,1,1;2022,4,1;2022,4,1;2022,5,1;2022,5,1;2022,9,1;2022,9,1;2023,4,1;2023,4,1;...
                                  2023,5,1;2023,5,1;2023,9,1;2023,9,1;2024,4,1;2024,4,1;2024,5,1;2024,5,1;2024,9,1;2024,9,1;2025,1,1;2026,1,1]);
            tvec_dep = datevec(float_tdep);
            
            for ivar = 1:numel(vtypes)-2
                disp(vnames_low{ivar})
                figure('Position',pos,'Visible','off')

                if ~isempty(wod_names.ctd{ivar})
                    ctd = dataset.(zones{izone}).wod.ctd.(vnames_low{ivar});
                else
                    ctd = dataset.(zones{izone}).wod.ctd.do*0;
                end
                if ~isempty(wod_names.osd{ivar})
                    osd = dataset.(zones{izone}).wod.osd.(vnames_low{ivar});
                else
                    osd = dataset.(zones{izone}).wod.osd.do*0;
                end
                ctd_avg0014 = mean(ctd(idx_0014));
                osd_avg0014 = mean(osd(idx_0014));
                % WARNING
                ctd(tvec(:,1)>2014) = ctd_avg0014;
                osd(tvec(:,1)>2014) = osd_avg0014;
                %
                ax1 = axes('Position',[0.115 0.1 0.65 0.85]); %[0.1 0.1 0.7 0.85],[0.1 0.1 0.625 0.85]
                h1 = area(time,[dataset.(zones{izone}).bio.(vnames_low{ivar}),ctd,osd]);
                hold on
                set(h1,'LineStyle','none')
                h3 = plot(time,dataset.(zones{izone}).bio.(vnames_low{ivar}),'k',...
                          time,dataset.(zones{izone}).bio.(vnames_low{ivar})+ctd,'k',...
                          time,dataset.(zones{izone}).bio.(vnames_low{ivar})+ctd+osd,'k');
                set(h3,'LineWidth',2)
% %                title([zones_names{izone},' - ',lnames{ivar}],'Interpreter','none');
%                 lg = legend(h1,'Argo','CTD','Bottle');
%                 set(lg,'Location','NorthWest','Box','Off')
                set(ax1,'XLim',datenum([1965,tvec(end,1)],1,1),'YLim',[0 ymax(ivar)],'XTick',datenum(1965:5:tvec(end,1),1,1),'YTick',0:500:3500,...
                        'XTickLabel',datestr(datenum(1965:5:tvec(end,1),1,1),'yy'),'FontWeight','bold','FontSize',fs1,'FontName','Arial','LineWidth',3,...
                        'Layer','Top')
                ylabel('# of stations')
%                 t1 = text(ax1.XLim(1)+diff(ax1.XLim)/2,ax1.YLim(1)+diff(ax1.YLim)/1.15,lnames{ivar});
                t1 = text(ax1.XLim(1)+diff(ax1.XLim)/1.25,ax1.YLim(1)+diff(ax1.YLim)/1.075,lnames{ivar});
                set(t1,'FontWeight','bold','FontSize',fs1,'FontName','Arial','HorizontalAlignment','center')
                ax2 = axes('Position',[0.785 0.1 0.2 0.85]); %[0.825 0.1 0.125 0.85],[0.75 0.1 0.2 0.85]
                h2 = area(timef,[nfloats',repmat(dataset.(zones{izone}).bio.(vnames_low{ivar})(end),[numel(timef),1]),repmat(ctd_avg0014,[numel(timef),1]),repmat(osd_avg0014,[numel(timef),1])]);
                % h2 = area(float_tdep,[float_ndep'*365/5,repmat(dataset.(zones{izone}).bio.(vnames_low{ivar})(end),[numel(float_tdep),1]),repmat(ctd_avg0014,[numel(float_tdep),1]),repmat(osd_avg0014,[numel(float_tdep),1])]);
                               
                h2(1).FaceColor = [0 0.75 1];
                for itype = 2:4
                    h2(itype).FaceColor = h1(itype-1).FaceColor;
                end
                hold on
                h4 = plot(timef,nfloats,'k--',...
                          timef,nfloats'+repmat(dataset.(zones{izone}).bio.(vnames_low{ivar})(end),[numel(timef),1]),'k--',...
                          timef,nfloats'+repmat(dataset.(zones{izone}).bio.(vnames_low{ivar})(end),[numel(timef),1])+repmat(ctd_avg0014,[numel(timef),1]),'k--',...
                          timef,nfloats'+repmat(dataset.(zones{izone}).bio.(vnames_low{ivar})(end),[numel(timef),1])+repmat(ctd_avg0014,[numel(timef),1])+...
                          repmat(osd_avg0014,[numel(timef),1]),'k--');
                set(h2,'LineStyle','none')
                set(h4,'LineWidth',2)
                set(ax2,'XLim',float_tdep([1 end]),'YLim',[0 ymax(ivar)],'XTick',datenum(unique(tvec_dep(1:end-1,1)),1,1),'YTick',0:500:3500,...
                        'XTickLabel',datestr(datenum(unique(tvec_dep(1:end-1,1)),1,1),'yy'),'FontWeight','bold','FontSize',fs1,'FontName','Arial','LineWidth',3,...
                        'Layer','Top','YtickLabel','')
                    
                    
                lg = legend(fliplr(h2),'Bottle','CTD','Argo','NWNA BGC Argo (this proposal)');
                set(lg,'Position',[0.2,0.75,0.3,0.1],'Box','Off')
                
%                 h4 = plot(float_tdep,float_ndep'*365/5,'k--',...
%                           float_tdep,float_ndep'*365/5'+repmat(dataset.(zones{izone}).bio.(vnames_low{ivar})(end),[numel(float_tdep),1]),'k--',...
%                           float_tdep,float_ndep'*365/5'+repmat(dataset.(zones{izone}).bio.(vnames_low{ivar})(end),[numel(float_tdep),1])+repmat(ctd_avg0014,[numel(float_tdep),1]),'k--',...
%                           float_tdep,float_ndep'*365/5'+repmat(dataset.(zones{izone}).bio.(vnames_low{ivar})(end),[numel(float_tdep),1])+repmat(ctd_avg0014,[numel(float_tdep),1])+...
%                           repmat(osd_avg0014,[numel(float_tdep),1]),'k--');
%                 set(h2,'LineStyle','none')
%                 set(h4,'LineWidth',2)
%                 set(ax2,'XLim',float_tdep([1 end]),'YLim',[0 ymax(ivar)],'XTick',datenum(unique(tvec_dep(1:end-1,1)),1,1),'YTick',0:500:3500,...
%                         'XTickLabel',datestr(datenum(unique(tvec_dep(1:end-1,1)),1,1),'yy'),'FontWeight','bold','FontSize',fs1,'FontName','Arial','LineWidth',3,...
%                         'Layer','Top','YtickLabel','')

                set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
                plot_name = ['CFI_ArgoWODCount_',zones{izone},'_1965_2025_year_area_',vnames_low{ivar}];
                print(gcf,'-depsc2','-painters','-r300',plot_name);
                print(gcf,'-dpng','-r150',plot_name);
                close
            end
        end
    end
end

if 0
    disp('Map with model grids coordinates')
	disp('get data')
	lakes = shaperead('/misc/10/output/laurenta/data/NaturalEarth/ne_10m_lakes/ne_10m_lakes.shp', 'UseGeoCoords', true);
	rivers = shaperead('/misc/10/output/laurenta/data/NaturalEarth/ne_10m_rivers_lake_centerlines/ne_10m_rivers_lake_centerlines.shp', 'UseGeoCoords', true);
    ice = shaperead('/misc/10/output/laurenta/data/NaturalEarth/ne_50m_glaciated_areas/ne_50m_glaciated_areas.shp', 'UseGeoCoords', true);

%     ncload('/misc/10/output/laurenta/data/ETOPO/etopo1_m120_10_30_80.nc')
    ncload('/misc/10/output/laurenta/data/ETOPO/etopo1_m100_10_20_80.nc')
    addpath('/misc/10/output/laurenta/data/cptcmap')

    [lattrk,lontrk] = track([latp;latp(1)],[lonp;lonp(1)]);

    Z = Band1;
    R = georasterref('RasterSize', size(Z), 'Latlim', [min(lat) max(lat)], 'Lonlim', [min(lon) max(lon)]);

 	disp('make figure')

    pos = [1 1250 1250 1250];
    axpos = [0.01 0.01 0.98 0.98];

    if 0
        figure('Position',pos,'Visible','off')

        ax1 = axesm('MapProjection','lambert','MapLatLimit',[32 80],'MapLonLimit',[-90 -30],'frame','on','grid','on');
        set(ax1,'Position',axpos,'XColor','none','YColor','none')
        geoshow(Z,R, 'DisplayType', 'texturemap');
        cptcmap('GMT_globe', 'mapping', 'direct');
        caxis([-10000 10000])

        hold on
        g2 = geoshow(lakes,'FaceColor',[0.9451    0.9882    1.0000],'EdgeColor','None','LineWidth',0.5);
        g3 = geoshow(rivers,'Color',[0.9451    0.9882    1.0000],'LineWidth',0.5);
        g4 = geoshow(ice,'FaceColor',[1 1 1],'EdgeColor','None');

        plotm(lattrk,lontrk,'k-','LineWidth',3)

        tightmap
        setm(ax1,'glinestyle','-','glinewidth',0.25,'gcolor',[0.1500 0.1500 0.1500],'MLineLocation',10,'PLineLocation',10,...
                 'ParallelLabel','off','MeridianLabel','off','MLabelLocation',-90:10:-30,'PLabelLocation',40:10:80,...
                 'LabelRotation','off','FontSize',fs,'FontWeight','bold','FontName','helvetica')

        disp('save figure')
        plotname = 'Map_CFI_AtlanticCanada_maplimits_narrow';
        set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
        print(gcf,'-dpng','-r100',plotname)
        close
    end

    if 0
        figure('Position',pos,'Visible','off')

        ax1 = axesm('MapProjection','lambert','MapLatLimit',[32 80],'MapLonLimit',[-90 -30],'frame','on','grid','on');
        set(ax1,'Position',axpos,'XColor','none','YColor','none','Color','none')

        plotm(lattrk,lontrk,'k-','LineWidth',2)

        tightmap
        setm(ax1,'glinestyle','-','glinewidth',0.25,'gcolor',[0.1500 0.1500 0.1500],'MLineLocation',10,'PLineLocation',10,...
                 'ParallelLabel','off','MeridianLabel','off','MLabelLocation',-90:10:-30,'PLabelLocation',40:10:80,...
                 'LabelRotation','off','FontSize',fs,'FontWeight','bold','FontName','helvetica')

        disp('save figure')
        plotname = 'Map_CFI_AtlanticCanada_limitsonly_narrow';
        set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
        print(gcf,'-depsc2','-painters',plotname);
        close
    end

    pos = [1 1500 2000 1500];

    if 0 % Map with AC region limits
        figure('Position',pos,'Visible','off')

        ax1 = axesm('MapProjection','lambert','MapLatLimit',[32 80],'MapLonLimit',[-100 10],'frame','on','grid','on');
        set(ax1,'Position',axpos,'XColor','none','YColor','none')
        geoshow(Z,R, 'DisplayType', 'texturemap');
        cptcmap('GMT_globe', 'mapping', 'direct');
        caxis([-10000 10000])

        hold on
        g2 = geoshow(lakes,'FaceColor',[0.9451    0.9882    1.0000],'EdgeColor','None','LineWidth',0.5);
        g3 = geoshow(rivers,'Color',[0.9451    0.9882    1.0000],'LineWidth',0.5);
        g4 = geoshow(ice,'FaceColor',[1 1 1],'EdgeColor','None');

        plotm(lattrk,lontrk,'k-','LineWidth',3)

        tightmap
        setm(ax1,'glinestyle','-','glinewidth',0.25,'gcolor',[0.1500 0.1500 0.1500],'MLineLocation',10,'PLineLocation',10,...
                 'ParallelLabel','off','MeridianLabel','off','MLabelLocation',-100:10:10,'PLabelLocation',40:10:80,...
                 'LabelRotation','off','FontSize',fs,'FontWeight','bold','FontName','helvetica')

        disp('save figure')
        plotname = 'Map_CFI_AtlanticCanada_maplimits_wide';
        set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
        print(gcf,'-dpng','-r100',plotname)
        close
    end

    if 0 % AC region limits only
        figure('Position',pos,'Visible','off')

        ax1 = axesm('MapProjection','lambert','MapLatLimit',[32 80],'MapLonLimit',[-100 10],'frame','on','grid','on');
        set(ax1,'Position',axpos,'XColor','none','YColor','none','Color','none')

        plotm(lattrk,lontrk,'k-','LineWidth',2)

        tightmap
        setm(ax1,'glinestyle','-','glinewidth',0.25,'gcolor',[0.1500 0.1500 0.1500],'MLineLocation',10,'PLineLocation',10,...
                 'ParallelLabel','off','MeridianLabel','off','MLabelLocation',-100:10:10,'PLabelLocation',40:10:80,...
                 'LabelRotation','off','FontSize',fs,'FontWeight','bold','FontName','helvetica')

        disp('save figure')
        plotname = 'Map_CFI_AtlanticCanada_limitsonly_wide';
        set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
        print(gcf,'-depsc2','-painters',plotname);
        close
    end

    if 0 % Map with Erin cruise + AZOMP + AR7W transect

        fname = '/misc/10/output/laurenta/data/AZMP/AR7W_stations.xlsx';
        AR7W = readtable(fname,'Sheet','AR7W_stations'); % AR07W is sampled once a year in May

        load /misc/1/home/laurenta/ROMS/workspace/matlab/acm/AZOMP_coord.mat

        % ID Lat    Lon    Depth
        Stations_EB = ...
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

        Station_EB_names = {'St. John''s','OSNAP','OSNAP','SC-1','','','','','','','','','SC-2','','','','','','SC-3','St. John''s'};

        disp('make figure')
%         pos = [1 935 1500 935]; % For long 30-80
        pos = [1 945 1500 945]; % For long 20-80
        axpos = [0.0 0.0 1 1];
        fs = 20;

        figure('Position',pos,'Visible','off')

        ax1 = axesm('MapProjection','lambert','MapLatLimit',[20 80],'MapLonLimit',[-100 10],'frame','on','grid','on');
        set(ax1,'Position',axpos,'XColor','none','YColor','none')
        geoshow(Z,R, 'DisplayType', 'texturemap');
        cptcmap('GMT_globe', 'mapping', 'direct');
        caxis([-10000 10000])

        hold on
        g2 = geoshow(lakes,'FaceColor',[0.9451    0.9882    1.0000],'EdgeColor','None','LineWidth',0.5);
        g3 = geoshow(rivers,'Color',[0.9451    0.9882    1.0000],'LineWidth',0.5);
        g4 = geoshow(ice,'FaceColor',[1 1 1],'EdgeColor','None');

        plotm(Stations_EB(2:end-1,2),Stations_EB(2:end-1,3),'m-','LineWidth',1)
        plotm(Stations_EB(4,2),Stations_EB(4,3),'ms','MarkerFaceColor','m','MarkerSize',15)
        plotm(AR7W.Latitude_N,-AR7W.Longitude_W,'ko','MarkerFaceColor','k','MarkerSize',6)
        plotm(Stations_EB([5:12,14:18],2),Stations_EB([5:12,14:18],3),'mo','MarkerFaceColor','m','MarkerSize',6)
        plotm(Stations_EB(2:3,2),Stations_EB(2:3,3),'ms','MarkerFaceColor','m','MarkerSize',8)
        %textm(Stations_EB([2,3,5:12,14:18],2),Stations_EB([2,3,5:12,14:18],3),num2cell(Stations_EB([2,3,5:12,14:18],1)),'HorizontalAlignment','Left','VerticalAlignment','Top')
        %plotm(Stations_EB(4,2),Stations_EB(4,3),'ms','MarkerSize',15)
        %plotm(Stations_EB(1,2),Stations_EB(1,3),'ro','MarkerFaceColor','r','MarkerSize',6)
        plotm(HL(1:8,1),HL(1:8,2),'ko','MarkerFaceColor','k','MarkerSize',6)
        plotm(HL(9:11,1),HL(9:11,2),'ko','MarkerFaceColor','k','MarkerSize',6)
        plotm(HL(12:14,1),HL(12:14,2),'ko','MarkerFaceColor','w','MarkerSize',6)

        tightmap
        setm(ax1,'glinestyle','-','glinewidth',0.25,'gcolor',[0.1500 0.1500 0.1500],'MLineLocation',10,'PLineLocation',10,...
                 'ParallelLabel','off','MeridianLabel','off','MLabelLocation',-100:10:10,'PLabelLocation',20:10:80,...
                 'LabelRotation','off','FontSize',fs,'FontWeight','bold','FontName','helvetica')

        disp('save figure')
        plotname = 'Cruise_map_CFI2';
        set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
        print(gcf,'-dpng','-r100',plotname)
        close
    end
    
    if 1 % Map with Erin cruise + AZOMP + AR7W transect (final version)

        fname = '/misc/10/output/laurenta/data/AZMP/AR7W_stations.xlsx';
        AR7W = readtable(fname,'Sheet','AR7W_stations'); % AR07W is sampled once a year in May

        load /misc/1/home/laurenta/ROMS/workspace/matlab/acm/AZOMP_coord.mat

        [lattrk,lontrk] = track([30,30,66.2],[-81.31,-35,-35]);

        % ID Lat    Lon    Depth
        Stations_EB = ...
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

        AZMP_ll = [-52.967,48.73;-52.750,48.800;-52.650,48.833;-52.400,48.917;-52.067,49.025;... % BON
                   -51.830,49.100;-51.542,49.190;-51.280,49.280;-51.017,49.367;-50.533,49.517;...% BON
                   -50.017,49.683;-49.500,49.850;-49.000,50.000;-48.472,50.177;-47.947,50.332;...% BON
                   -59.850,45.825;-59.702,45.658;-59.517,45.492;-59.175,45.158;-58.850,44.817;...% LOU
                   -58.508,44.475;-58.175,44.133;-57.833,43.783;-57.527,43.473;...               % LOU
             -52.933,46.583;-52.850,46.500;-52.733,46.350;-52.608,46.208;...                % SEG
             -52.500,46.070;-52.267,45.788;-52.000,45.458;-51.700,45.095;...                % SEG
             -51.395,44.725;-51.103,44.363;-50.808,44.000;-50.517,43.633;...                % SEG
             -50.167,43.200;-49.942,42.920;-49.887,42.850;-49.828,42.775;...                % SEG
             -49.683,42.588;-49.518,42.395;-49.270,42.082;-48.950,41.700;-48.667,41.333;... % SEG
             -52.832,47.000;-52.705,47.000;-52.580,47.000;-52.322,47.000;-52.033,47.000;... % FC
             -51.485,47.000;-51.000,47.000;-50.667,47.000;-50.000,47.000;-49.117,47.000;... % FC
             -52.832,47.000;-52.705,47.000;-52.580,47.000;-52.322,47.000;-52.033,47.000;... % FC
             -51.485,47.000;-51.000,47.000;-50.667,47.000;-50.000,47.000;-49.117,47.000;... % FC
             -48.617,47.000;-48.117,47.000;-47.817,47.000;-47.500,47.000;-47.250,47.000;... % FC
             -47.168,47.000;-47.017,47.000;-46.833,47.000;-46.670,47.000;-46.483,47.000;... % FC
             -46.017,47.000;-45.730,47.000;-45.500,47.000;-45.213,47.000;-44.988,47.000;... % FC
             -44.772,47.000;-44.578,47.000;-44.433,47.000;-44.232,47.000;-44.083,47.000;... % FC
             -43.833,47.000;-43.750,47.000;-43.400,47.000;-43.250,47.000;-43.000,47.000;... % FC
             -42.750,47.000;-42.500,47.000;-42.000,47.000];                                 % FC

        Station_EB_names = {'St. John''s','OSNAP','OSNAP','SC-1','','','','','','','','','SC-2','','','','','','SC-3','St. John''s'};

        disp('make figure')
        pos = [1 830 1500 830]; % For long 20-70
        axpos = [0.0 0.0 1 1];
        fs = 20;

        figure('Position',pos,'Visible','off')

        ax1 = axesm('MapProjection','lambert','MapLatLimit',[20 70],'MapLonLimit',[-100 10],'frame','on','grid','on');
        set(ax1,'Position',axpos,'XColor','none','YColor','none')
        geoshow(Z,R, 'DisplayType', 'texturemap');
        cptcmap('GMT_globe', 'mapping', 'direct');
        caxis([-10000 10000])

        hold on
        g2 = geoshow(lakes,'FaceColor',[0.9451    0.9882    1.0000],'EdgeColor','None','LineWidth',0.5);
        g3 = geoshow(rivers,'Color',[0.9451    0.9882    1.0000],'LineWidth',0.5);
        g4 = geoshow(ice,'FaceColor',[1 1 1],'EdgeColor','None');

        plotm(Stations_EB(2:end-1,2),Stations_EB(2:end-1,3),'m-','LineWidth',1)
        plotm(Stations_EB(4,2),Stations_EB(4,3),'ms','MarkerFaceColor','m','MarkerSize',15)
        plotm(AR7W.Latitude_N,-AR7W.Longitude_W,'ko','MarkerFaceColor','k','MarkerSize',6)
        plotm(Stations_EB([5:12,14:18],2),Stations_EB([5:12,14:18],3),'mo','MarkerFaceColor','m','MarkerSize',6)
        plotm(Stations_EB(2:3,2),Stations_EB(2:3,3),'ms','MarkerFaceColor','m','MarkerSize',8)
        plotm(HL(1:8,1),HL(1:8,2),'ko','MarkerFaceColor','k','MarkerSize',6)
        plotm(HL(9:11,1),HL(9:11,2),'ko','MarkerFaceColor','k','MarkerSize',6)
        plotm(HL(12:14,1),HL(12:14,2),'ko','MarkerFaceColor','k','MarkerSize',6)
        plotm(AZMP_ll(:,2),AZMP_ll(:,1),'ko','MarkerFaceColor','k','MarkerSize',6)
        plotm(lattrk,lontrk,'k-','LineWidth',3)
%        plotm([62,60.3],[-66.1,-64.8],'k-','LineWidth',3)
%load coastlines
%plotm(coastlat,coastlon,'k-')
        tightmap
        setm(ax1,'glinestyle','-','glinewidth',0.25,'gcolor',[0.1500 0.1500 0.1500],'MLineLocation',10,'PLineLocation',10,...
                 'ParallelLabel','off','MeridianLabel','off','MLabelLocation',-100:10:10,'PLabelLocation',20:10:70,...
                 'LabelRotation','off','FontSize',fs,'FontWeight','bold','FontName','helvetica')

        disp('save figure')
        plotname = 'Cruise_map_CFI_v3';
        set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
        print(gcf,'-dpng','-r100',plotname)
        close
    end
end
