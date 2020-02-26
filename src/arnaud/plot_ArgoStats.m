clearvars
zones = {'NAtl','AC_LB_BB','SO_50','Global'};% Atlantic Canada, LS+BB; Southern Ocean (>50S); Global case
zones_names = {'North Atlantic','Atlantic Canada','Southern Ocean','Global'};
date_start_argo = datenum(1997,1,1);
date_start = datenum(1950,1,1);
date_stop = datenum(2019,9,1);

vtypes = {'DOXY','NITRATE','PH_IN_SITU_TOTAL','CHLA','BBP700','DOWNWELLING_PAR','CDOM'};
vnames = {'OXYGEN','NITRATE','PH','CHLA','BBP700','PAR','CDOM'};
vnames_low = {'do','no3','ph','chl','bbp','par','cdom'};

wod_names.osd = {'Oxygen','Nitrate','pH','Chlorophyll','','',''};
wod_names.ctd = {'Oxygen','','pH','Chlorophyll','','',''};            
wod_type = {'osd','ctd'};

oceans = {'A','I','P'};
floatid_active = [];
for iocean = 1:numel(oceans)
    infoFile = urlread(['http://www.ifremer.fr/co-argoFloats/floatList?detail=false&active=true&ocean=',oceans{iocean}]);
    expr = '&amp;ptfCode=';
    k = strfind(infoFile,expr);
    for ifloat = 1:numel(k)
        str = infoFile(k(ifloat)+numel(expr):k(ifloat)+numel(expr)+20);
        floatid_active(end+1,1) = str2num(str(regexp(str,'[0-9]')));
    end
    clear infoFile
end

sareakm2.Global = 361132000; % km2

%%% Old sensors (missing light)
% sensor_types = [1  1  1  1  1;... % 5, All sensors
%                 1  1  0  1  1;... % 4, no pH
%                 1  1  1  0  0;... % 3, no chl, bbp
%                 1  0  0  1  1;... % 3, no nitrate, pH
%                 1  0  1  0  0;... % 2, O2, pH
%                 0  0  0  1  1;... % 2, chl, bbp
%                 1  0  0  0  0;... % 1, O2
%                 0  0  0  0  0];   % 0, Others, no main sensors

%%% New sensors (including light)

sensor_types = [ 1     1     1     1     1     1;... % 6, All sensors
                 1     1     1     1     1     0;... % 5, No light
                 1     1     0     1     1     1;... % 5, no pH
                 1     1     0     1     1     0;... % 4, no pH, light
                 1     0     0     1     1     1;... % 4, no nitrate, no pH
                 1     1     1     0     0     0;... % 3, O2, nitrate, ph
                 1     0     0     1     1     0;... % 3, O2,chl, bbp
                 0     0     0     1     1     1;... % 3, chl, bbp, light
                 1     0     1     0     0     0;... % 2, O2, pH
                 0     0     0     1     1     0;... % 2, chl, bbp
                 1     0     0     0     0     0];   % 1, O2
   
sensor_types_names = {'6 sensors (do, no3, ph, chl, bbp, light)',...
                 '5 sensors (do, no3, ph, chl, bbp)','5 sensors (do, no3, chl, bbp, light)',...
                 '4 sensors (do, no3, chl, bbp)','4 sensors (do, chl, bbp, light)',...
                 '3 sensors (O2, nitrate, ph)','3 sensors (O2, chl, bbp)','3 sensors (chl, bbp, light)',...
                 '2 sensors (O2, pH)','2 sensors (chl, bbp)','O2 only'};
             
%%% This was calculated as follows:
% instruments = [];k = 1;izone = 4;
% fname = ['ArgoCount_',zones{izone},'_',datestr(date_start_argo,'yyyymmdd'),'_',datestr(date_stop,'yyyymmdd')];
% load(fname)
% nsensors = sum(data_float.bio.types(:,1:6),2);
% disp('add floats')
% for ifloat = 1:numel(data_float.bio.traj)
%     if ismember(data_float.bio.names(ifloat),floatid_active)&nsensors(ifloat)>0
%         instruments(k,:) = data_float.bio.types(ifloat,1:6);
%         k=k+1;
%     end
% end
% flipud(unique(instruments,'rows'))

disp('Set region')

NA_ll = [32  32      40        41 50  80   80   32;
         -100 -4.6474 -4.6474 -0.5 10  10 -100 -100];
    
transect_xy = [

   217    84
   217   188
   209   188
   209   334
   196   337
   193   341
   167   344
   160   344
   160   342
   157   331
   154   327
   155   323
   154   319
   140   319
   140   264
   124   239
   133   224
   134   224
   134   221
   135   219
   136   219
   136   217
   135   209
   135   164
    86   136
    86    84
    95    84];

WOD = load('WODstats_data','dataset');
load /misc/1/home/laurenta/NEMO/workspace/matlab/mat/mask_nemo3oceans
for i0 = 1:numel(transect_xy(:,1))
    latp(i0,1) = lat(transect_xy(i0,1),transect_xy(i0,2));
    lonp(i0,1) = lon(transect_xy(i0,1),transect_xy(i0,2));
end

for itype = 1:numel(wod_type)
    %wod_names_osd = {'Oxygen','Nitrate','pH','Chlorophyll','',''};
    %wod_names_ctd = {'Oxygen','','pH','','',''};
    for ivar = 1:numel(wod_names.(wod_type{itype}))
        if ~isempty(wod_names.(wod_type{itype}){ivar})
            % Get indices of all files with user-specified polygon and timeframe:
            lona = WOD.dataset.(wod_type{itype}).(wod_names.(wod_type{itype}){ivar})(:,2);
            lata = WOD.dataset.(wod_type{itype}).(wod_names.(wod_type{itype}){ivar})(:,3);
            [in1,on1] = inpolygon(lona,lata,NA_ll(2,:),NA_ll(1,:));
            [in2,on2] = inpolygon(lona,lata,lonp,latp);
            WOD_inpoly{1}.(wod_type{itype}).(wod_names.(wod_type{itype}){ivar}) = in1+on1;
            WOD_inpoly{2}.(wod_type{itype}).(wod_names.(wod_type{itype}){ivar}) = in2+on2;
            WOD_inpoly{3}.(wod_type{itype}).(wod_names.(wod_type{itype}){ivar}) = lata<-50;
            WOD_inpoly{4}.(wod_type{itype}).(wod_names.(wod_type{itype}){ivar}) = lona*0+1;
        end
    end
end
    
fname = '/misc/10/output/laurenta/data/ETOPO/ETOPO1_Ice_g_gmt4.grd';
disp(['Load ',fname])
ncload(fname)
[X,Y] = meshgrid(x,y);
R = georasterref('RasterSize', size(z), 'Latlim', [min(y) max(y)], 'Lonlim', [min(x) max(x)]);

izone = 1;
[in1,on1] = inpolygon(X,Y,NA_ll(2,:),NA_ll(1,:));
NAg = (in1+on1)&z<0;
sareakm2.(zones{izone}) = areamat(NAg, R, wgs84Ellipsoid('kilometers'));

izone = 3;
SOg = Y<-50&z<0;
sareakm2.(zones{izone}) = areamat(SOg, R, wgs84Ellipsoid('kilometers'));
clearvars x y z

fs = 30;
pos = [1 800 1600 800];

if 0
    for izone = 1:numel(zones)
        fname = ['ArgoCount_',zones{izone},'_',datestr(date_start_argo,'yyyymmdd'),'_',datestr(date_stop,'yyyymmdd')];
        disp(['Load ',fname])
        load(fname)

        tvec = datevec(t);
        tvec_start = datevec(date_start);
        tvec_stop = datevec(date_stop);
        years = tvec_start(1):tvec_stop(1);

        dataset = [];
        itime = 1;
        for iyear = 1:numel(years)
            for imonth = 1:12
                dataset.(zones{izone}).glob.all(itime,1) = sum(count.glob.all(tvec(:,1)==years(iyear)&tvec(:,2)==imonth));
                dataset.(zones{izone}).glob.time(itime,1) = datenum(years(iyear),imonth,1);
                for itype = 1:numel(vtypes)
                    dataset.(zones{izone}).bio.(vnames_low{itype})(itime,1) = sum(count.bio.(vnames_low{itype})(tvec(:,1)==years(iyear)&tvec(:,2)==imonth));
                    if itype==1
                        dataset.(zones{izone}).bio.time(itime,1) = datenum(years(iyear),imonth,1);
                    end
                end
                itime = itime + 1;
            end
        end
        dataset.(zones{izone}).glob.time(dataset.(zones{izone}).glob.time>datenum(2019,8,31)) = nan;
        dataset.(zones{izone}).bio.time(dataset.(zones{izone}).bio.time>datenum(2019,8,31)) = nan;
        figure('Position',pos,'Visible','off')
        h1 = plot(dataset.(zones{izone}).bio.time,dataset.(zones{izone}).bio.do,'o-',dataset.(zones{izone}).glob.time,dataset.(zones{izone}).glob.all,'o-');
        set(h1,'LineWidth',2)
        set(h1(1),'MarkerFaceColor',h1(1).Color)
        set(h1(2),'MarkerFaceColor',h1(2).Color)
        ax1 = gca;
        title(zones_names{izone},'Interpreter','none');
        lg = legend(h1,'BGC ARGO','ARGO');
        set(lg,'Location','NorthWest','Box','Off')
        set(ax1,'Position',[0.1 0.1 0.85 0.85],'XLim',[date_start,datenum(tvec_stop(1)+1,1,1)],'XTick',datenum([years,tvec_stop(1)+1],1,1),...
                'XTickLabel',datestr(datenum([years,tvec_stop(1)+1],1,1),'yy'),'FontWeight','bold','FontSize',fs,'FontName','Arial','LineWidth',3,...
                'Layer','Top')
        set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
        print(gcf,'-dpng','-r300',fname)
    end
    clear h1
end

if 1
    for izone = 2%[1,3,4]
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
			for ivar = 1:numel(vtypes)
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
            figure('Position',pos,'Visible','off')
            h1 = plot(dataset.(zones{izone}).bio.time,dataset.(zones{izone}).bio.(vnames_low{1}),'o-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).bio.(vnames_low{2}),'o-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).bio.(vnames_low{3}),'o-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).bio.(vnames_low{4}),'o-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).bio.(vnames_low{5}),'o-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).bio.(vnames_low{6}),'o-');
            set(h1,'LineWidth',2)
            set(h1(1),'MarkerFaceColor',h1(1).Color)
            set(h1(2),'MarkerFaceColor',h1(2).Color)
            set(h1(3),'MarkerFaceColor',h1(3).Color)
            set(h1(4),'MarkerFaceColor',h1(4).Color)
            set(h1(5),'MarkerFaceColor',h1(5).Color)
            set(h1(6),'MarkerFaceColor',h1(6).Color)
            ax1 = gca;
            title(zones_names{izone},'Interpreter','none');
            lg = legend(h1,vnames{1:6});
            set(lg,'Location','NorthWest','Box','Off')
            set(ax1,'Position',[0.1 0.1 0.85 0.85],'XLim',[date_start_argo,datenum(tvec_stop(1)+1,1,1)],'XTick',datenum([years,tvec_stop(1)+1],1,1),...
                    'XTickLabel',datestr(datenum([years,tvec_stop(1)+1],1,1),'yy'),'FontWeight','bold','FontSize',fs,'FontName','Arial','LineWidth',3,...
                    'Layer','Top')
            set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
            plot_name = [fname,'_year_6sensors'];
            %print(gcf,'-dpng','-r300',plot_name) 
            print(gcf,'-depsc2','-painters','-r300',plot_name);
        end
        
        if 0
            figure('Position',pos,'Visible','off')
            h1 = plot(dataset.(zones{izone}).bio.time,dataset.(zones{izone}).bio.(vnames_low{1}),'o-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).bio.(vnames_low{2}),'o-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).bio.(vnames_low{3}),'o-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).bio.(vnames_low{4}),'o-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).bio.(vnames_low{5}),'o-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).bio.(vnames_low{6}),'o-');
            hold on

            h2 = plot(dataset.(zones{izone}).bio.time,dataset.(zones{izone}).wod.ctd.(vnames_low{1}),'^-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).wod.ctd.(vnames_low{1})*nan,'^-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).wod.ctd.(vnames_low{3}),'^-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).wod.ctd.(vnames_low{4}),'^-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).wod.ctd.(vnames_low{1})*nan,'^-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).wod.ctd.(vnames_low{1})*nan,'^-');

            h3 = plot(dataset.(zones{izone}).bio.time,dataset.(zones{izone}).wod.osd.(vnames_low{1}),'s-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).wod.osd.(vnames_low{2}),'s-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).wod.osd.(vnames_low{3}),'s-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).wod.osd.(vnames_low{4}),'s-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).wod.osd.(vnames_low{1})*nan,'s-',...
                      dataset.(zones{izone}).bio.time,dataset.(zones{izone}).wod.osd.(vnames_low{1})*nan,'s-');

            set(h1,'LineWidth',2)
            set([h2,h3],'LineWidth',1)
            for ivar = 1:numel(vtypes)
                set(h1(ivar),'MarkerFaceColor',h1(ivar).Color)
    %             set(h2(ivar),'Color',h1(ivar).Color,'MarkerFaceColor',h1(ivar).Color)
    %             set(h3(ivar),'Color',h1(ivar).Color,'MarkerFaceColor',h1(ivar).Color)
                set(h2(ivar),'Color',h1(ivar).Color)
                set(h3(ivar),'Color',h1(ivar).Color)
            end
            ax1 = gca;
            title(zones_names{izone},'Interpreter','none');
            lg = legend(h1,vnames);
            set(lg,'Location','NorthWest','Box','Off')
            set(ax1,'Position',[0.1 0.1 0.85 0.85],'XLim',datenum([1965,tvec_stop(1)],1,1),'XTick',datenum([1965:5:years(end),tvec_stop(1)+1],1,1),...
                    'XTickLabel',datestr(datenum([1965:5:years(end),tvec_stop(1)+1],1,1),'yyyy'),'FontWeight','bold','FontSize',fs,'FontName','Arial','LineWidth',3,...
                    'Layer','Top')
            set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
            plot_name = ['ArgoWODCount_',zones{izone},'_',datestr(datenum(1965,1,1),'yyyymmdd'),'_',datestr(date_stop,'yyyymmdd'),'_year'];
            %print(gcf,'-dpng','-r300',plot_name);
            print(gcf,'-depsc2','-painters','-r300',plot_name);
        end
        
        if 0 % This was used for Katja's OceanObs19 poster
            % Numer of stations
            for ivar = 1:numel(vtypes)
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
                h1 = area(dataset.(zones{izone}).bio.time,[dataset.(zones{izone}).bio.(vnames_low{1}),ctd,osd]);

                set(h1,'LineWidth',1)
                ax1 = gca;
                title(zones_names{izone},'Interpreter','none');
                lg = legend(h1,'ARGO','CTD','Bottle');
                set(lg,'Location','NorthWest','Box','Off')
                set(ax1,'Position',[0.1 0.1 0.85 0.85],'XLim',datenum([1965,tvec_stop(1)],1,1),'XTick',datenum([1965:5:years(end),tvec_stop(1)+1],1,1),...
                        'XTickLabel',datestr(datenum([1965:5:years(end),tvec_stop(1)+1],1,1),'yyyy'),'FontWeight','bold','FontSize',fs,'FontName','Arial','LineWidth',3,...
                        'Layer','Top','YLim',[0 100])
                ylabel('# of stations')
                set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
                plot_name = ['ArgoWODCount_',zones{izone},'_',datestr(datenum(1965,1,1),'yyyymmdd'),'_',datestr(date_stop,'yyyymmdd'),'_year_area_',vnames_low{ivar}];
                %print(gcf,'-dpng','-r300',plot_name);
                print(gcf,'-depsc2','-painters','-r300',plot_name);
            end
        end
        
        if 1 % New figure for CFI2020 proposal
            % Numer of stations
            for ivar = 1:numel(vtypes)
                return
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
                h1 = area(dataset.(zones{izone}).bio.time,[dataset.(zones{izone}).bio.(vnames_low{1}),ctd,osd]);

                set(h1,'LineWidth',1)
                ax1 = gca;
                title(zones_names{izone},'Interpreter','none');
                lg = legend(h1,'ARGO','CTD','Bottle');
                set(lg,'Location','NorthWest','Box','Off')
                set(ax1,'Position',[0.1 0.1 0.85 0.85],'XLim',datenum([1965,tvec_stop(1)],1,1),'XTick',datenum([1965:5:years(end),tvec_stop(1)+1],1,1),...
                        'XTickLabel',datestr(datenum([1965:5:years(end),tvec_stop(1)+1],1,1),'yyyy'),'FontWeight','bold','FontSize',fs,'FontName','Arial','LineWidth',3,...
                        'Layer','Top','YLim',[0 100])
                ylabel('# of stations')
                set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
                plot_name = ['CFI_ArgoWODCount_',zones{izone},'_',datestr(datenum(1965,1,1),'yyyymmdd'),'_',datestr(date_stop,'yyyymmdd'),'_year_area_',vnames_low{ivar}];
                %print(gcf,'-dpng','-r300',plot_name);
                print(gcf,'-depsc2','-painters','-r300',plot_name);
            end
        end
        
        % Numer of stations/area
        if 0%1
            for ivar = 1:numel(vtypes)
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
                h1 = area(dataset.(zones{izone}).bio.time,[dataset.(zones{izone}).bio.(vnames_low{1}),ctd,osd]/(sareakm2.(zones{izone})/1e6));

                set(h1,'LineWidth',1)
                ax1 = gca;
                title(zones_names{izone},'Interpreter','none');
                lg = legend(h1,'ARGO','CTD','Bottle');
                set(lg,'Location','NorthWest','Box','Off')
                set(ax1,'Position',[0.1 0.1 0.85 0.85],'XLim',datenum([1965,tvec_stop(1)],1,1),'XTick',datenum([1965:5:years(end),tvec_stop(1)+1],1,1),...
                        'XTickLabel',datestr(datenum([1965:5:years(end),tvec_stop(1)+1],1,1),'yyyy'),'FontWeight','bold','FontSize',fs,'FontName','Arial','LineWidth',3,...
                        'Layer','Top','YLim',[0 100])
                ylabel('# of stations per million km2')
                set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
                plot_name = ['ArgoWODCount_',zones{izone},'_',datestr(datenum(1965,1,1),'yyyymmdd'),'_',datestr(date_stop,'yyyymmdd'),'_year_areakm2_y100_',vnames_low{ivar}];
                %print(gcf,'-dpng','-r300',plot_name);
                print(gcf,'-depsc2','-painters','-r300',plot_name);
            end
        end
    end
    clear h1 h2
end

if 0
    figure('Position',pos,'Visible','off')
    for izone = 1:numel(zones)
        disp(['Load ',fname])
        load(fname)

        tvec = datevec(t);
        tvec_start = datevec(date_start);
        tvec_stop = datevec(date_stop);
        years = tvec_start(1):tvec_stop(1);

        dataset = [];
        itime = 1;
        for iyear = 1:numel(years)
            for imonth = 1:12
                dataset.(zones{izone}).glob.all(itime,1) = sum(count.glob.all(tvec(:,1)==years(iyear)&tvec(:,2)==imonth));
                dataset.(zones{izone}).glob.time(itime,1) = datenum(years(iyear),imonth,1);
                itime = itime + 1;
            end
        end
        dataset.(zones{izone}).glob.time(dataset.(zones{izone}).glob.time>datenum(2019,8,31)) = nan;

        h1(izone) = plot(dataset.(zones{izone}).glob.time,dataset.(zones{izone}).glob.all,'o-','LineWidth',2);
        hold on
    end

    set(h1(1),'MarkerFaceColor',h1(1).Color)
    set(h1(2),'MarkerFaceColor',h1(2).Color)
    set(h1(3),'MarkerFaceColor',h1(3).Color)
    ax1 = gca;
    title(zones_names{izone},'Interpreter','none');
    lg = legend(h1,zones_names);
    set(lg,'Location','NorthWest','Box','Off')
    set(ax1,'Position',[0.1 0.1 0.85 0.85],'XLim',[date_start,datenum(tvec_stop(1)+1,1,1)],'XTick',datenum([years,tvec_stop(1)+1],1,1),...
            'XTickLabel',datestr(datenum([years,tvec_stop(1)+1],1,1),'yy'),'FontWeight','bold','FontSize',fs,'FontName','Arial','LineWidth',3,...
            'Layer','Top')
    set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
    print(gcf,'-dpng','-r300',['ArgoCount_',datestr(date_start,'yyyymmdd'),'_',datestr(date_stop,'yyyymmdd')])
    clear h1
end

if 0
    figure('Position',pos,'Visible','off')
    for izone = 1:numel(zones)
        disp(['Load ',fname])
        load(fname)

        tvec = datevec(t);
        tvec_start = datevec(date_start);
        tvec_stop = datevec(date_stop);
        years = tvec_start(1):tvec_stop(1);

        dataset = [];
        itime = 1;
        for iyear = 1:numel(years)
            for imonth = 1:12
                for itype = 1:numel(vtypes)
                    dataset.(zones{izone}).bio.(vnames_low{itype})(itime,1) = sum(count.bio.(vnames_low{itype})(tvec(:,1)==years(iyear)&tvec(:,2)==imonth));
                    if itype==1
                        dataset.(zones{izone}).bio.time(itime,1) = datenum(years(iyear),imonth,1);
                    end
                end
                itime = itime + 1;
            end
        end
        dataset.(zones{izone}).bio.time(dataset.(zones{izone}).bio.time>datenum(2019,8,31)) = nan;

        h1(izone) = plot(dataset.(zones{izone}).bio.time,dataset.(zones{izone}).bio.do,'o-','LineWidth',2);
        hold on
    end

    set(h1,'LineWidth',2)
    set(h1(1),'MarkerFaceColor',h1(1).Color)
    set(h1(2),'MarkerFaceColor',h1(2).Color)
    set(h1(3),'MarkerFaceColor',h1(3).Color)
    ax1 = gca;
    title(zones_names{izone},'Interpreter','none');
    lg = legend(h1,zones_names);
    set(lg,'Location','NorthWest','Box','Off')
    set(ax1,'Position',[0.1 0.1 0.85 0.85],'XLim',[date_start,datenum(tvec_stop(1)+1,1,1)],'XTick',datenum([years,tvec_stop(1)+1],1,1),...
            'XTickLabel',datestr(datenum([years,tvec_stop(1)+1],1,1),'yy'),'FontWeight','bold','FontSize',fs,'FontName','Arial','LineWidth',3,...
            'Layer','Top')
    set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
    print(gcf,'-dpng','-r300',['ArgoCount_BGC_',datestr(date_start,'yyyymmdd'),'_',datestr(date_stop,'yyyymmdd')])
end

if 0
%   Map only
	disp('get data')
	lakes = shaperead('/misc/10/output/laurenta/data/NaturalEarth/ne_10m_lakes/ne_10m_lakes.shp', 'UseGeoCoords', true);
	rivers = shaperead('/misc/10/output/laurenta/data/NaturalEarth/ne_10m_rivers_lake_centerlines/ne_10m_rivers_lake_centerlines.shp', 'UseGeoCoords', true);
    ice = shaperead('/misc/10/output/laurenta/data/NaturalEarth/ne_50m_glaciated_areas/ne_50m_glaciated_areas.shp', 'UseGeoCoords', true);

    % ETOPO1 data available here: https://maps.ngdc.noaa.gov/viewers/wcs-client/
    ncload('/misc/10/output/laurenta/data/ETOPO/etopo1_m120_10_30_80.nc')
    addpath('/misc/10/output/laurenta/data/cptcmap')
%  	[X,Y] = meshgrid(lon,lat);
 
    Z = Band1;
    R = georasterref('RasterSize', size(Z), 'Latlim', [min(lat) max(lat)], 'Lonlim', [min(lon) max(lon)]);

 	disp('make figure')
    pos = [1 1500 2000 1500];
    axpos = [0.05 0.09 0.9 0.88];
    figure('Position',pos,'Visible','off')
    
%    ax1 = axesm('MapProjection','ortho','Origin',[60 -65],'FLatLimit',[-Inf 20]); % _20b version
    ax1 = axesm('MapProjection','lambert','MapLatLimit',[32 80],'MapLonLimit',[-100 10]);
    set(ax1,'XColor','none','YColor','none')
%    ax1 = axesm('MapProjection','gstereo','MapLatLimit',[32 80],'MapLonLimit',[-100 10]);
    framem; gridm('on');
%     geoshow(Y,X,Band1, 'DisplayType', 'texturemap');
    geoshow(Z,R, 'DisplayType', 'texturemap');
    cptcmap('GMT_globe', 'mapping', 'direct');
    caxis([-10000 10000])
%     mlabel; 
%     plabel;
    
    hold on
   % shadem([-30, 50])
    g2 = geoshow(lakes,'FaceColor',[0.9451    0.9882    1.0000],'EdgeColor','None','LineWidth',0.5);
    g3 = geoshow(rivers,'Color',[0.9451    0.9882    1.0000],'LineWidth',0.5);
    g4 = geoshow(ice,'FaceColor',[1 1 1],'EdgeColor','None');
    
    tightmap
    setm(ax1,'glinestyle','-','glinewidth',0.25,'gcolor',[0.1500 0.1500 0.1500])
    
    disp('save figure')
    plotname = ['test_map_v1c'];
    set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
    print(gcf,'-dpng','-r100',plotname)
    close
end

if 0
%   Map with model grids coordinates
    disp('Map with model grids coordinates')
    % Parent model: 34N-75N, 40W-80W. Please use black dashed lines.
    % Child model for the Labrador Sea: 52N-63N, 45W-76W. Please use red dashed lines.
    % Child model for the Slope Water: 38N-48N, 50W-66.3W Please use red dashed lines
    glim = [-80,-40,34,75;-76,-45,52,63;-66.3,-50,38,48];
    for igrid = 1:size(glim,1)
        [glat.(['g',num2str(igrid)]), glon.(['g',num2str(igrid)])] = outlinegeoquad(glim(igrid,3:4),glim(igrid,1:2),0.1,0.1);
    end
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
    pos = [1 1250 2000 1250];
    axpos = [0.01 0.01 0.98 0.98];
    figure('Position',pos,'Visible','off')
    
    ax1 = axesm('MapProjection','lambert','MapLatLimit',[32 80],'MapLonLimit',[-100 10],'frame','on','grid','on');
    set(ax1,'Position',axpos,'XColor','none','YColor','none')
%     framem; gridm('on');
    geoshow(Z,R, 'DisplayType', 'texturemap');
    cptcmap('GMT_globe', 'mapping', 'direct');
    caxis([-10000 10000])

    hold on
    g2 = geoshow(lakes,'FaceColor',[0.9451    0.9882    1.0000],'EdgeColor','None','LineWidth',0.5);
    g3 = geoshow(rivers,'Color',[0.9451    0.9882    1.0000],'LineWidth',0.5);
    g4 = geoshow(ice,'FaceColor',[1 1 1],'EdgeColor','None');
    % Quadrangle covering Australia and vicinity
    g5 = geoshow(glat.g1,glon.g1,'DisplayType','polygon','FaceColor','none','EdgeColor','k','LineStyle','--','LineWidth',4);
    g6 = geoshow(glat.g2,glon.g2,'DisplayType','polygon','FaceColor','none','EdgeColor','r','LineStyle','--','LineWidth',4);
    g7 = geoshow(glat.g3,glon.g3,'DisplayType','polygon','FaceColor','none','EdgeColor','r','LineStyle','--','LineWidth',4);

    tightmap
    setm(ax1,'glinestyle','-','glinewidth',0.25,'gcolor',[0.1500 0.1500 0.1500])
    
    disp('save figure')
    plotname = 'NA_map_grids';
    set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
    print(gcf,'-dpng','-r100',plotname)
    close
end

if 0
%   Map with trajectories of active BGC floats (all variables, NA region)
	disp('get data')
	lakes = shaperead('/misc/10/output/laurenta/data/NaturalEarth/ne_10m_lakes/ne_10m_lakes.shp', 'UseGeoCoords', true);
	rivers = shaperead('/misc/10/output/laurenta/data/NaturalEarth/ne_10m_rivers_lake_centerlines/ne_10m_rivers_lake_centerlines.shp', 'UseGeoCoords', true);
    ice = shaperead('/misc/10/output/laurenta/data/NaturalEarth/ne_50m_glaciated_areas/ne_50m_glaciated_areas.shp', 'UseGeoCoords', true);

    ncload('/misc/10/output/laurenta/data/ETOPO/etopo1_m120_10_30_80.nc')
    addpath('/misc/10/output/laurenta/data/cptcmap')
 
    Z = Band1;
    R = georasterref('RasterSize', size(Z), 'Latlim', [min(lat) max(lat)], 'Lonlim', [min(lon) max(lon)]);
    
 	disp('make figure')
    pos = [1 1500 2000 1500];
    axpos = [0.05 0.09 0.9 0.88];
    figure('Position',pos,'Visible','off')
    
    ax1 = axesm('MapProjection','lambert','MapLatLimit',[32 80],'MapLonLimit',[-100 10]);
    set(ax1,'XColor','none','YColor','none')
    framem; gridm('on');
    geoshow(Z,R, 'DisplayType', 'texturemap');
    cptcmap('GMT_globe', 'mapping', 'direct');
    caxis([-10000 10000])
    hold on
    g2 = geoshow(lakes,'FaceColor',[0.9451    0.9882    1.0000],'EdgeColor','None','LineWidth',0.5);
    g3 = geoshow(rivers,'Color',[0.9451    0.9882    1.0000],'LineWidth',0.5);
    g4 = geoshow(ice,'FaceColor',[1 1 1],'EdgeColor','None');

    tightmap
    setm(ax1,'glinestyle','-','glinewidth',0.25,'gcolor',[0.1500 0.1500 0.1500])

    cm = cbrewer('div','Spectral',size(sensor_types,1)+3);
    cm = cm(1:end-3,:);
	izone = 1;
	fname = ['ArgoCount_',zones{izone},'_',datestr(date_start_argo,'yyyymmdd'),'_',datestr(date_stop,'yyyymmdd')];
	load(fname)
    nsensors = sum(data_float.bio.types(:,1:6),2);
    nfloats = 0;
	for ifloat = 1:numel(data_float.bio.traj)
		if ismember(data_float.bio.names(ifloat),floatid_active)&nsensors(ifloat)>0
			p1 = plotm(data_float.bio.traj{ifloat}(:,3),data_float.bio.traj{ifloat}(:,2),'-','Color',[1 1 1]*0.25,'LineWidth',0.5);
            nfloats = nfloats + 1;
		end
    end
    disp(['Total active BGC floats on ',zones{izone},' map: ',num2str(nfloats)])
    idx = [];
	for ifloat = 1:numel(data_float.bio.traj)
		if ismember(data_float.bio.names(ifloat),floatid_active)&nsensors(ifloat)>0
            idx(end+1,1) = strmatch(num2str(data_float.bio.types(ifloat,1:6)),num2str(sensor_types));
            p2 = plotm(data_float.bio.traj{ifloat}(end,3),data_float.bio.traj{ifloat}(end,2),'ko','MarkerFaceColor',cm(idx(end),:),...
                       'MarkerSize',11);
		end
    end
    ax2 = axes('visible','off');
    colormap(ax2,cm)
    caxis([0 1])
    cb = colorbar('Location','South');
    set(cb,'Position',[0.775 0.21 0.125 0.02],'YTick',0.1:0.2:0.9,'YTickLabel',[],'FontWeight','bold','FontSize',25,'FontName','Arial',...
        'YAxisLocation','Top')
    ylabel(cb,'Sensor types')
    disp('save figure')
    plotname = ['map_traj_bgc_all_active_',zones{izone},'_v1e_6sensors'];
    set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
    print(gcf,'-dpng','-r100',plotname)
    close
    
    figure('Position',pos,'Visible','off')
    colormap(flipud(cm(unique(idx),:)))
    cb = colorbar;
    set(cb,'Position',[0.5 0.1, 0.02, 0.8],'YTick',1/(numel(unique(idx)))/2:1/(numel(unique(idx))):(1-1/(numel(unique(idx)))/2),'TickLabels',fliplr(sensor_types_names(unique(idx))),...
           'FontWeight','bold','FontSize',25,'FontName','Arial','LineWidth',1)
    set(gca,'Visible','off')
    plotname = ['cm_',plotname];
    set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
    print(gcf,'-dpng','-r100',plotname)
    close
end

if 0
%   Map with trajectories of active BGC floats (all variables, Global Ocean)
	disp('get shape data')
	lakes = shaperead('/misc/10/output/laurenta/data/NaturalEarth/ne_10m_lakes/ne_10m_lakes.shp', 'UseGeoCoords', true);
	rivers = shaperead('/misc/10/output/laurenta/data/NaturalEarth/ne_10m_rivers_lake_centerlines/ne_10m_rivers_lake_centerlines.shp', 'UseGeoCoords', true);
    ice = shaperead('/misc/10/output/laurenta/data/NaturalEarth/ne_50m_glaciated_areas/ne_50m_glaciated_areas.shp', 'UseGeoCoords', true);

    disp('load etopo')
    ncload('/misc/10/output/laurenta/data/ETOPO/ETOPO1_Ice_g_gmt4.grd')
    addpath('/misc/10/output/laurenta/data/cptcmap')
    disp('georasterref etopo data')
    R = georasterref('RasterSize', size(z), 'Latlim', [min(y) max(y)], 'Lonlim', [min(x) max(x)]);
    
 	disp('make figure')
    res = 1;
    pos = [1 [1500 2000 1500]*res];
    axpos = [0.05 0.09 0.9 0.88];
    figure('Position',pos,'Visible','off')
    
    ax1 = axesm('MapProjection','mollweid');%,'MapLatLimit',[32 80],'MapLonLimit',[-100 10]);
    set(ax1,'XColor','none','YColor','none')
    framem; gridm('on');
    geoshow(z,R, 'DisplayType', 'texturemap');
    cptcmap('GMT_globe', 'mapping', 'direct');
    caxis([-10000 10000])
    hold on
    g2 = geoshow(lakes,'FaceColor',[0.9451    0.9882    1.0000],'EdgeColor','None','LineWidth',0.25);
    g3 = geoshow(rivers,'Color',[0.9451    0.9882    1.0000],'LineWidth',0.25);
    g4 = geoshow(ice,'FaceColor',[1 1 1],'EdgeColor','None');

    tightmap
    setm(ax1,'glinestyle','-','glinewidth',0.25,'gcolor',[0.1500 0.1500 0.1500])

    cm = cbrewer('div','Spectral',size(sensor_types,1));
	izone = 4;
	fname = ['ArgoCount_',zones{izone},'_',datestr(date_start_argo,'yyyymmdd'),'_',datestr(date_stop,'yyyymmdd')];
	load(fname)
    nsensors = sum(data_float.bio.types(:,1:5),2);
    disp('add floats')
    nfloats = 0;
	for ifloat = 1:numel(data_float.bio.traj)
		if ismember(data_float.bio.names(ifloat),floatid_active)&nsensors>0
            idx = strmatch(num2str(data_float.bio.types(ifloat,1:5)),num2str(sensor_types));
            p2 = plotm(data_float.bio.traj{ifloat}(end,3),data_float.bio.traj{ifloat}(end,2),'ko','MarkerFaceColor',cm(idx,:),...
                       'MarkerSize',11);
            nfloats = nfloats + 1;
		end
    end
    disp(['Total active BGC floats on ',zones{izone},' map: ',num2str(nfloats)])
    ax2 = axes('visible','off');
    colormap(ax2,cm)
    caxis([0 1])
    cb = colorbar('Location','South');
    set(cb,'Position',[0.775 0.15 0.15 0.02],'YTick',0.1:0.2:0.9,'YTickLabel','','FontWeight','bold','FontSize',20,'FontName','Arial',...
        'YAxisLocation','Top')
    ylabel(cb,'Sensor types')
    disp('save figure')
    plotname = ['map_traj_bgc_all_active_',zones{izone},'_v1g'];
    set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
    print(gcf,'-dpng','-r100',plotname)
    close
end

if 0
%   Map with trajectories of active BGC floats, continents in grey (all variables, Global Ocean)
    disp('load etopo')
    ncload('/misc/10/output/laurenta/data/ETOPO/ETOPO1_Ice_g_gmt4.grd')
    addpath('/misc/10/output/laurenta/data/cptcmap')
    disp('georasterref etopo data')
    R = georasterref('RasterSize', size(z), 'Latlim', [min(y) max(y)], 'Lonlim', [min(x) max(x)]);
    
 	disp('make figure')
    res = 1;
    pos = [1 [1500 2000 1500]*res];
    axpos = [0.05 0.09 0.9 0.88];
    figure('Position',pos,'Visible','off')

    ax1 = axesm('MapProjection','mollweid');%,'MapLatLimit',[32 80],'MapLonLimit',[-100 10]);
    set(ax1,'XColor','none','YColor','none')
    framem; gridm('on');
    geoshow(z,R, 'DisplayType', 'texturemap');
    cm = cptcmap('GMT_globe', 'mapping', 'direct');
    cm(129:end,1) = 0.75;cm(129:end,2) = 0.75;cm(129:end,3) = 0.75;
    colormap(cm)
    caxis([-10000 10000])
    hold on
%     g2 = geoshow(lakes,'FaceColor',[0.9451    0.9882    1.0000],'EdgeColor','None','LineWidth',0.25);
%     g3 = geoshow(rivers,'Color',[0.9451    0.9882    1.0000],'LineWidth',0.25);
%     g4 = geoshow(ice,'FaceColor',[1 1 1],'EdgeColor','None');

    tightmap
    setm(ax1,'glinestyle','-','glinewidth',0.25,'gcolor',[0.1500 0.1500 0.1500])

    cm = cbrewer('div','Spectral',size(sensor_types,1)+3);
    cm = cm(1:end-3,:);
	izone = 4;
	fname = ['ArgoCount_',zones{izone},'_',datestr(date_start_argo,'yyyymmdd'),'_',datestr(date_stop,'yyyymmdd')];
	load(fname)
    nsensors = sum(data_float.bio.types(:,1:6),2);
    disp('add floats')
    idx = [];
    nfloats = 0;
	for ifloat = 1:numel(data_float.bio.traj)
		if ismember(data_float.bio.names(ifloat),floatid_active)&nsensors(ifloat)>0
            idx(end+1,1) = strmatch(num2str(data_float.bio.types(ifloat,1:6)),num2str(sensor_types));
            p2 = plotm(data_float.bio.traj{ifloat}(end,3),data_float.bio.traj{ifloat}(end,2),'ko','MarkerFaceColor',cm(idx(end),:),...
                       'MarkerSize',11);
            nfloats = nfloats + 1;
		end
    end
    disp(['Total active BGC floats on ',zones{izone},' map: ',num2str(nfloats)])
    ax2 = axes('visible','off');
    colormap(ax2,cm)
    caxis([0 1])
    cb = colorbar('Location','South');
    set(cb,'Position',[0.775 0.15 0.15 0.02],'YTick',0.1:0.2:0.9,'YTickLabel','','FontWeight','bold','FontSize',20,'FontName','Arial',...
        'YAxisLocation','Top')
    ylabel(cb,'Sensor types')
    disp('save figure')
    plotname = ['map_traj_bgc_all_active_',zones{izone},'_v1h_6sensors'];
    set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
    print(gcf,'-dpng','-r100',plotname)
    close
    
    figure('Position',pos,'Visible','off')
    colormap(flipud(cm(unique(idx),:)))
    cb = colorbar;
    set(cb,'Position',[0.5 0.1, 0.02, 0.8],'YTick',1/(numel(unique(idx)))/2:1/(numel(unique(idx))):(1-1/(numel(unique(idx)))/2),'TickLabels',fliplr(sensor_types_names(unique(idx))),...
           'FontWeight','bold','FontSize',25,'FontName','Arial','LineWidth',1)
    set(gca,'Visible','off')
    plotname = ['cm_',plotname];
    set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
    print(gcf,'-dpng','-r100',plotname)
    close
end

if 0
%   Map with trajectories of all BGC floats (active/inactive, all variables, AC region)
	izone = 1;
	fname = ['ArgoCount_',zones{izone},'_',datestr(date_start_argo,'yyyymmdd'),'_',datestr(date_stop,'yyyymmdd')];
	load(fname)
	disp('get data')
	lakes = shaperead('/misc/10/output/laurenta/data/NaturalEarth/ne_10m_lakes/ne_10m_lakes.shp', 'UseGeoCoords', true);
	rivers = shaperead('/misc/10/output/laurenta/data/NaturalEarth/ne_10m_rivers_lake_centerlines/ne_10m_rivers_lake_centerlines.shp', 'UseGeoCoords', true);
    ice = shaperead('/misc/10/output/laurenta/data/NaturalEarth/ne_50m_glaciated_areas/ne_50m_glaciated_areas.shp', 'UseGeoCoords', true);
    ncload('/misc/10/output/laurenta/data/ETOPO/etopo1_m120_10_30_80.nc')
    addpath('/misc/10/output/laurenta/data/cptcmap')
 	[X,Y] = meshgrid(lon,lat);
    
 	disp('make figure')
    pos = [1 1500 2000 1500];
    axpos = [0.05 0.09 0.9 0.88];
    figure('Position',pos,'Visible','off')
    
    ax1 = axesm('MapProjection','gstereo','MapLatLimit',[32 80],'MapLonLimit',[-100 10]);
    set(ax1,'XColor','none','YColor','none')
    framem; gridm('off');
    geoshow(Y,X,Band1, 'DisplayType', 'texturemap');
    cptcmap('GMT_globe', 'mapping', 'direct');
    caxis([-10000 10000])

    hold on
    g2 = geoshow(lakes,'FaceColor',[0.9451    0.9882    1.0000],'EdgeColor','None','LineWidth',0.5);
    g3 = geoshow(rivers,'Color',[0.9451    0.9882    1.0000],'LineWidth',0.5);
    g4 = geoshow(ice,'FaceColor',[1 1 1],'EdgeColor','None');

    cm1 = lines(300);
	for ifloat = 1:numel(data_float.bio.traj)
        p1 = plotm(data_float.bio.traj{ifloat}(:,3),data_float.bio.traj{ifloat}(:,2),'.-','Color',cm1(ifloat,:),'LineWidth',0.25);
    end

    disp('save figure')
    plotname = ['map_traj_bgc_all_',zones{izone},'_v1'];
    set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'points', 'PaperPosition', pos);
    print(gcf,'-dpng','-r100',plotname)
    close
end