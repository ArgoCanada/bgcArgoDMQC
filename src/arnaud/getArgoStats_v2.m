clearvars;%close all

process = 1;
plt = 0;
NA = 0; % North Atlantic
LB = 0; % Labrador Sea
AC = 0; % Atlantic Canada, LS+BB
NWA = 1; % NorthWest North Atlantic, this is the region used for CFI proposal
SO1 = 0;% Southern Ocean (>45S)
SO2 = 0;% Southern Ocean (>50S)
GB = 0; % Global case

date_start = datenum(1997,1,1);
date_stop = datenum(2019,9,1);

update_files = 0;
fglobal = 'ar_index_global_prof.txt';
fbio = 'argo_bio-profile_index.txt';
greylist = 'ar_greylist.txt';

fprintf(['Processing info:\n',datestr(date_start,'yyyymmdd'),'-',datestr(date_stop,'yyyymmdd'),...
      '\nUpdate = ',num2str(update_files),'\nProcessing = ',num2str(process),'\n\n'])
if update_files
    disp('Updating floats info files')
    ftpobj = ftp('usgodae.org');
    cd(ftpobj,'/pub/outgoing/argo/');
    disp(['Downloading ',fglobal])
    mget(ftpobj,fglobal);
    disp(['Downloading ',fbio])
    mget(ftpobj,fbio);
    close(ftpobj)
    disp('Done with updates')
end

vtypes = {'DOXY','NITRATE','PH_IN_SITU_TOTAL','CHLA','BBP700','DOWNWELLING_PAR','CDOM'};
vnames = {'OXYGEN','NITRATE','PH','CHLA','BBP700','DOWNWELLING_PAR','CDOM'};
vnames_low = {'do','no3','ph','chl','bbp','par','cdom'};
% Active floats in Atlantic Basin
disp('Reading info on active floats')
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

% Import the month's data file for BGC Argo:
% file,date,latitude,longitude,ocean,profiler_type,institution,parameters,param_data_mode,date_update
disp('Load data for BGC Argo')
fid = fopen(fbio);
% 1:fname 2:date 3:lat 4:lon 5:basin 6:variables
Dbio = textscan(fid,'%s %s %f %f %s %*f %*s %s %*s %*f','headerlines',9,'delimiter',',', 'whitespace', '');
fclose(fid);
variables_bio = Dbio{6};
floatid_bio = regexp(Dbio{1},'/(\w*)/','tokens');
for ifloat = 1:numel(floatid_bio)
    float_bio(ifloat,1) = str2double(cell2mat(floatid_bio{ifloat}{1}));
end
tmp_bio = Dbio{2};
lata_bio = Dbio{3};
lona_bio = Dbio{4};
basin_bio = Dbio{5};

% Remove empty rows:
ind_bio = cellfun(@isempty,tmp_bio);
float_bio(ind_bio) = [];
variables_bio(ind_bio) = [];
lata_bio(ind_bio) = [];
lona_bio(ind_bio) = [];
tmp_bio(ind_bio) = [];
basin_bio(ind_bio) = [];
% Convert date strings to datenums:
time_bio = datenum(tmp_bio,'yyyymmddHHMMSS');

% Import the month's data file for Global Argo:
% file,date,latitude,longitude,ocean,profiler_type,institution,date_update
disp('Load data for global Argo')
fid = fopen(fglobal);
% 1:fname 2:date 3:lat 4:lon 5:basin
Dglob = textscan(fid,'%s %s %f %f %s %*f %*s %*f','headerlines',9,'delimiter',',', 'whitespace', '');
fclose(fid);
floatid_glob = regexp(Dglob{1},'/(\w*)/','tokens');
for ifloat = 1:numel(floatid_glob)
    float_glob(ifloat,1) = str2double(cell2mat(floatid_glob{ifloat}{1}));
end
tmp_glob = Dglob{2};
lata_glob = Dglob{3};
lona_glob = Dglob{4};
basin_glob = Dglob{5};

% Remove empty rows:
ind_glob = cellfun(@isempty,tmp_glob);
float_glob(ind_glob) = [];
lata_glob(ind_glob) = [];
lona_glob(ind_glob) = [];
tmp_glob(ind_glob) = [];
basin_glob(ind_glob) = [];
% Convert date strings to datenums:
time_glob = datenum(tmp_glob,'yyyymmddHHMMSS');


% The following command was used to get transect_xy
% [gridout,transect_xy] = transect_gridpoints_nemo_v2;
% Lines 216, 217 and 220 of the script were commented

disp('Set region')
if NA
% Limits for North Atlantic lat: -100,10; lon: 32,80 including Labrador Sea and Baffin Bay
zone = 'NAtl';
NA_ll = [32  32      40        41 50  80   80   32;
         -100 -4.6474 -4.6474 -0.5 10  10 -100 -100];

elseif LB
% Limits for Labrador Sea and Baffin Bay
zone = 'LB_BB';
transect_xy = [

   131   157
   157   143
   209   190
   202   215
   204   244
   215   276
   219   307
   200   333
   194   340
   179   351
   173   348
   162   344
   148   340
   132   327
   142   303
   142   267
   128   254
   124   237
   131   227
   136   219
   135   208
   128   188
   114   173];
elseif AC
% Limits for Atlantic region plus Labrador Sea and Baffin Bay
zone = 'AC_LB_BB';
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
elseif NWA
    % Limits for North West North Atlantic including Labrador Sea and southern part of Baffin Bay
    zone = 'NWAtl';
    NWA_ll = [ 30   70   70   66.7053   63.3478   62.0954   62.1360   61.7930   61.6015   61.6404   61.4073   60.4241   54.6531   49.7625   32   30   30;...
              -35  -35  -75.50  -68.8382  -69.9071  -66.3152  -66.0641  -65.8167  -65.4046  -65.1520  -64.9955  -64.6588  -62.2532  -72.7457  -83  -83  -35];
elseif SO1
    zone = 'SO_45';
elseif SO2
    zone = 'SO_50';
else
    zone = 'Global';
end

if process
    disp(['Processing BGC Argo for zone ',zone])
    % BGC ARGO
    if GB
        inpoly = lona_bio*0+1;
    elseif SO1
        inpoly = lata_bio<-45;
    elseif SO2
        inpoly = lata_bio<-50;
    elseif NA
        [inpoly,onpoly] = inpolygon(lona_bio,lata_bio,NA_ll(2,:),NA_ll(1,:));
        inpoly = inpoly+onpoly;
    elseif NWA
        [inpoly,onpoly] = inpolygon(lona_bio,lata_bio,NWA_ll(2,:),NWA_ll(1,:));
        inpoly = inpoly+onpoly;
    else
        load /misc/1/home/laurenta/NEMO/workspace/matlab/mat/mask_nemo3oceans
        for i0 = 1:numel(transect_xy(:,1))
            latp(i0,1) = lat(transect_xy(i0,1),transect_xy(i0,2));
            lonp(i0,1) = lon(transect_xy(i0,1),transect_xy(i0,2));
        end
        % Get indices of all files with user-specified polygon and timeframe:
        [inpoly,onpoly] = inpolygon(lona_bio,lata_bio,lonp,latp);
        inpoly = inpoly+onpoly;
    end

    data_bio = sortrows([time_bio,inpoly,float_bio,lona_bio,lata_bio,(1:numel(time_bio))'],1);
    var_bio = variables_bio(data_bio(:,6));
    var_bio(data_bio(:,2)==0) = [];
    data0_bio = data_bio;
    data_bio(data_bio(:,2)==0,:) = [];

    t = date_start:date_stop;
    count.bio.all = nan(numel(t),1);
    for ivar = 1:numel(vnames_low)
        count.bio.(vnames_low{ivar}) = zeros(numel(t),1);
    end
    floats_bio = [];
    for itime = 1:numel(t)
        disp(datestr(t(itime)))
        ind = find(data_bio(:,1) >= t(itime)& data_bio(:,1) < t(itime)+1);
        flt_bio = unique(data_bio(ind,3));
        vname_bio = {};
        for ifloat = 1:numel(flt_bio)
            idx = find(float_bio==flt_bio(ifloat),1,'first');
            vname_bio{ifloat} = regexp(variables_bio{idx}, ' ', 'split');
            for ivar = 1:numel(vnames_low)
                 count.bio.(vnames_low{ivar})(itime) = count.bio.(vnames_low{ivar})(itime) + ismember(vtypes(ivar),vname_bio{ifloat});
            end
        end
        count.bio.all(itime) = numel(flt_bio);
        floats_bio = [floats_bio;flt_bio];
    end
    data_float.bio.names = unique(floats_bio);
    data_float.bio.types = nan(numel(data_float.bio.names),5);
    data_float.bio.all = [];
    for ifloat = 1:numel(data_float.bio.names)
        idx1 = find(data0_bio(:,3)==data_float.bio.names(ifloat));
        data_float.bio.traj{ifloat} = data0_bio(idx1,[1,4,5]);
        idx2 = find(float_bio==data_float.bio.names(ifloat),1,'first');
        data_float.bio.var{ifloat} = regexp(variables_bio{idx2}, ' ', 'split');
        data_float.bio.traj{ifloat}(abs(data_float.bio.traj{ifloat}(:,2))>360|abs(data_float.bio.traj{ifloat}(:,3))>90,2:3) = nan;
        % Manual QC
        if data_float.bio.names(ifloat)==4901805
            data_float.bio.traj{ifloat}(data_float.bio.traj{ifloat}(:,2)<-64,[2,3]) = nan;
        elseif data_float.bio.names(ifloat)==6902671
            data_float.bio.traj{ifloat}(data_float.bio.traj{ifloat}(:,1)>737014&data_float.bio.traj{ifloat}(:,1)<737283,[2,3]) = nan;
        elseif data_float.bio.names(ifloat)==6902896
            data_float.bio.traj{ifloat}(data_float.bio.traj{ifloat}(:,1)>737362,[2,3]) = nan;
        elseif data_float.bio.names(ifloat)==6902953
            data_float.bio.traj{ifloat}(data_float.bio.traj{ifloat}(:,1)>737365&data_float.bio.traj{ifloat}(:,1)<737612,[2,3]) = nan;
        end
        for itype = 1:numel(vtypes)
            data_float.bio.types(ifloat,itype) = ismember(vtypes(itype),data_float.bio.var{ifloat});
        end
        data_float.bio.all = [data_float.bio.all;data_float.bio.traj{ifloat}];
    end

    % GLOBAL ARGO
    disp(['Processing Global Argo for zone ',zone])
    if GB
        inpoly = lona_glob*0+1;
    elseif SO1
        inpoly = lata_glob<-45;
    elseif SO2
        inpoly = lata_glob<-50;
    elseif NA
        [inpoly,onpoly] = inpolygon(lona_glob,lata_glob,NA_ll(2,:),NA_ll(1,:));
        inpoly = inpoly+onpoly;
    elseif NWA
        [inpoly,onpoly] = inpolygon(lona_glob,lata_glob,NWA_ll(2,:),NWA_ll(1,:));
        inpoly = inpoly+onpoly;
    else
        load /misc/1/home/laurenta/NEMO/workspace/matlab/mat/mask_nemo3oceans
        for i0 = 1:numel(transect_xy(:,1))
            latp(i0,1) = lat(transect_xy(i0,1),transect_xy(i0,2));
            lonp(i0,1) = lon(transect_xy(i0,1),transect_xy(i0,2));
        end
        % Get indices of all files with user-specified polygon and timeframe:
        [inpoly,onpoly] = inpolygon(lona_glob,lata_glob,lonp,latp);
        inpoly = inpoly+onpoly;
    end

    data_glob = sortrows([time_glob,inpoly,float_glob,lona_glob,lata_glob,(1:numel(time_glob))'],1);
    data0_glob = data_glob;
    data_glob(data_glob(:,2)==0,:) = [];

    t = date_start:date_stop;
    count.glob.all = nan(numel(t),1);
    for ivar = 1:numel(vnames_low)
        count.glob.(vnames_low{ivar}) = zeros(numel(t),1);
    end
    floats_glob = [];
    for itime = 1:numel(t)
        disp(datestr(t(itime)))
        ind = find(data_glob(:,1) >= t(itime)& data_glob(:,1) < t(itime)+1);
        flt_glob = unique(data_glob(ind,3));
        count.glob.all(itime) = numel(flt_glob);
        floats_glob = [floats_glob;flt_glob];
    end
    data_float.glob.names = unique(floats_glob);
    data_float.glob.types = nan(numel(data_float.glob.names),5);
    data_float.glob.all = [];
    for ifloat = 1:numel(data_float.glob.names)
        idx1 = find(data0_glob(:,3)==data_float.glob.names(ifloat));
        data_float.glob.traj{ifloat} = data0_glob(idx1,[1,4,5]);
        data_float.glob.traj{ifloat}(abs(data_float.glob.traj{ifloat}(:,2))>360|abs(data_float.glob.traj{ifloat}(:,3))>90,2:3) = nan;
        % Manual QC
        if data_float.glob.names(ifloat)==4901805
            data_float.glob.traj{ifloat}(data_float.glob.traj{ifloat}(:,2)<-64,[2,3]) = nan;
        elseif data_float.glob.names(ifloat)==6902671
            data_float.glob.traj{ifloat}(data_float.glob.traj{ifloat}(:,1)>737014&data_float.glob.traj{ifloat}(:,1)<737283,[2,3]) = nan;
        elseif data_float.glob.names(ifloat)==6902896
            data_float.glob.traj{ifloat}(data_float.glob.traj{ifloat}(:,1)>737362,[2,3]) = nan;
        elseif data_float.glob.names(ifloat)==6902953
            data_float.glob.traj{ifloat}(data_float.glob.traj{ifloat}(:,1)>737365&data_float.glob.traj{ifloat}(:,1)<737612,[2,3]) = nan;
        end
%         for itype = 1:numel(vtypes)
%             data_float.glob.types(ifloat,itype) = ismember(vtypes(itype),data_float.glob.var{ifloat});
%         end
        data_float.glob.all = [data_float.glob.all;data_float.glob.traj{ifloat}];
    end
    fname = ['ArgoCount_',zone,'_',datestr(date_start,'yyyymmdd'),'_',datestr(date_stop,'yyyymmdd')];
    disp(['Save ',fname])
    save(fname,'count','t','data_float')
else
    fname = ['ArgoCount_',zone,'_',datestr(date_start,'yyyymmdd'),'_',datestr(date_stop,'yyyymmdd')];
    disp(['Load ',fname])
    load(fname)
end
if ~plt
    return
end
load coastlines

if 0
    figure;%set(gca,'XLim',[-75 -40],'YLim',[50 80])
    set(gca,'XLim',[-75 0],'YLim',[30 80])
    hold on
    for ifloat = 1:numel(data_float.traj)
        plot(data_float.traj{ifloat}(:,2),data_float.traj{ifloat}(:,3),'.','MarkerSize',2)
    end
    plot(coastlon,coastlat,'k-')
    amerc

    idx = data_float.all(:,1)>datenum(2019,1,1);
    figure;
    set(gca,'XLim',[-75 0],'YLim',[30 80])
    hold on
    scatter(data_float.all(idx,2),data_float.all(idx,3),5,data_float.all(idx,1),'Filled')
    plot(coastlon,coastlat,'k-')
    cb = colorbar;
    datetick(cb)
end

pos = [1 400 2000 400];
fs = 20;
axpos = [0.03 0.1 0.13 0.8;0.19 0.1 0.13 0.8;0.35 0.1 0.13 0.8;0.51 0.1 0.13 0.8;0.67 0.1 0.13 0.8;0.83 0.1 0.13 0.8];

if 0
    % Trajectories of BGC floats
    figure('Position',pos)
    for itype = 1:numel(vtypes)
        ax = axes('Position',axpos(itype,:));
        hold on
        nf = 0;
        for ifloat = 1:numel(data_float.bio.traj)
            if data_float.bio.types(ifloat,itype)
                plot(data_float.bio.traj{ifloat}(:,2),data_float.bio.traj{ifloat}(:,3),'.','MarkerSize',2)
                nf = nf+1;
            end
        end
        set(ax,'XLim',[-75 0],'YLim',[30 80],'FontSize',fs,'FontWeight','bold','FontName','Arial',...
                'LineWidth',1.75,'Layer','Top','Box','on')
        plot(coastlon,coastlat,'k-')
        title([vnames{itype},', #floats (total): ',num2str(nf)])
        amerc
    end
end

pos = [1 400 2000 400];
fs = 20;
axpos = [0.03 0.1 0.13 0.8;0.19 0.1 0.13 0.8;0.35 0.1 0.13 0.8;0.51 0.1 0.13 0.8;0.67 0.1 0.13 0.8;0.83 0.1 0.13 0.8];
if 0
    figure('Position',pos)
    for itype = 1:numel(vtypes)
        ax = axes('Position',axpos(itype,:));
        hold on
        nf = 0;
        for ifloat = 1:numel(data_float.bio.traj)
            if data_float.bio.types(ifloat,itype)
                if ismember(data_float.bio.names(ifloat),floatid_active)
                    p1 = plot(data_float.bio.traj{ifloat}(:,2),data_float.bio.traj{ifloat}(:,3),'-');
                    p2 = plot(data_float.bio.traj{ifloat}(end,2),data_float.bio.traj{ifloat}(end,3),'ko');
                    set(p2,'MarkerFaceColor',p1.Color)
                    nf = nf+1;
                end
            end
        end
        set(ax,'XLim',[-75 0],'YLim',[30 80],'FontSize',fs,'FontWeight','bold','FontName','Arial',...
                'LineWidth',1.75,'Layer','Top','Box','on')
        plot(coastlon,coastlat,'k-')
        title([vnames{itype},', #floats (active): ',num2str(nf)])
        if itype>1
            set(ax,'YTickLabel','')
        end
    end
end

pos = [1 400 2000 400];
fs = 20;
if 0
    axpos = [0.03 0.1 0.13 0.8;0.19 0.1 0.13 0.8;0.35 0.1 0.13 0.8;0.51 0.1 0.13 0.8;0.67 0.1 0.13 0.8;0.83 0.1 0.13 0.8];
    figure('Position',pos)
    for itype = 1:numel(vtypes)
        ax = axes('Position',axpos(itype,:));
        p1 = plot(t,count.bio.(vnames_low{itype}));
        set(ax,'XLim',datenum([2000,2020],1,1),'YLim',[0 10],'FontSize',fs,'FontWeight','bold','FontName','Arial',...
               'LineWidth',1.75,'Layer','Top','Box','on','XTick',datenum(2000:2020,1,1),'XTickLabel',...
               {'2000','','','','','2005','','','','','2010','','','','','2015','','','','','2020'})
        title(vnames{itype})
    end
end
%%
tvec = datevec(t);
tvec_start = datevec(date_start);
tvec_stop = datevec(date_stop);
years = tvec_start(1):tvec_stop(1);
dataset = [];
itime = 1;
for iyear = 1:numel(years)
    for imonth = 1:12
        dataset.glob.all(itime,1) = sum(count.glob.all(tvec(:,1)==years(iyear)&tvec(:,2)==imonth));
        dataset.glob.time(itime,1) = datenum(years(iyear),imonth,1);
        for itype = 1:numel(vtypes)
            dataset.bio.(vnames_low{itype})(itime,1) = sum(count.bio.(vnames_low{itype})(tvec(:,1)==years(iyear)&tvec(:,2)==imonth));
            if itype==1
                dataset.bio.time(itime,1) = datenum(years(iyear),imonth,1);
            end
        end
        itime = itime + 1;
    end
end

dataset.glob.time(dataset.glob.time>datenum(2019,8,31)) = nan;
dataset.bio.time(dataset.bio.time>datenum(2019,8,31)) = nan;
figure
plot(dataset.bio.time,dataset.bio.do,'o-',dataset.glob.time,dataset.glob.all,'o-')
title(zone,'Interpreter','none');
datetick
