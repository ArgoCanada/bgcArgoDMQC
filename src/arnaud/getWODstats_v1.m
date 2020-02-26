clearvars;close all
varnames = {'Temperature','Salinity','Oxygen','Nitrate','Phosphate','Chlorophyll','Silicate','pH','Alkalinity','tCO2','pCO2','DeltaC14','DeltaC13','CFC11','CFC12','Tritium','Helium','DeltaHe3','Neon'};
WODdir = '/misc/10/output/laurenta/data/WOD/';

process = 1;

% The cast count for your request is:
% 3229232	OSD casts
% 1071438	CTD casts
% 4300670	TOTAL casts

if process
    osd_files = dir([WODdir,'ocldb1565204366.26375_OSD*.nc']);
    ctd_files = dir([WODdir,'ocldb1565204366.26375_CTD*.nc']);

    % Check variable names in osd files
    diminfo = {};
    for ifile = 1:numel(osd_files)
        ncid = netcdf.open([osd_files(ifile).folder,'/',osd_files(ifile).name],'NOWRITE');
        finfo = ncinfo([osd_files(ifile).folder,'/',osd_files(ifile).name]); % Get information on the file
        diminfo = cat(2,diminfo,{finfo.Dimensions.Name});
        netcdf.close(ncid);
    end
    diminfo = unique(diminfo);
    Index = find(contains(diminfo,'_obs'));
    varnames_osd = erase(diminfo(Index),'_obs');
% Alkalinity Argon CFC113 CFC11 CFC12 Chlorophyll DeltaC13 DeltaC14 DeltaHe3 Helium Neon Nitrate
% Oxy18 Oxygen Phosphate Pressure Salinity Silicate Temperature Tritium pCO2 pH tCO2 z
    
    
    % Check variable names in ctd files
    diminfo = {};
    for ifile = 1:numel(ctd_files)
        ncid = netcdf.open([ctd_files(ifile).folder,'/',ctd_files(ifile).name],'NOWRITE');
        finfo = ncinfo([ctd_files(ifile).folder,'/',ctd_files(ifile).name]); % Get information on the file
        diminfo = cat(2,diminfo,{finfo.Dimensions.Name});
        netcdf.close(ncid);
    end
    diminfo = unique(diminfo);
    Index = find(contains(diminfo,'_obs'));
    varnames_ctd = erase(diminfo(Index),'_obs');
% Chlorophyll Nitrate Oxygen Pressure Salinity Temperature Transmissiv pH z
    
    % Process water samples
    for ifile = 1:numel(osd_files)
        disp(['Working on file ',osd_files(ifile).folder,'/',osd_files(ifile).name])
        ncid = netcdf.open([osd_files(ifile).folder,'/',osd_files(ifile).name],'NOWRITE'); % Open Netcdf file
        finfo = ncinfo([osd_files(ifile).folder,'/',osd_files(ifile).name]); % Get information on the file
        lon = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'),'double');   % degrees east
        lat = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'),'double');   % degrees north
        time = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'),'double'); % days since 1770-01-01 00:00:00 UTC
        time = time + datenum(1770,1,1);
        % Get data 
        for ivar = 1:numel(varnames_osd)
            if ifile==1
                dataset.osd.(varnames_osd{ivar}) = [];
            end
            if ismember({[varnames_osd{ivar},'_obs']},{finfo.Dimensions.Name})
                disp(['   - ',varnames_osd{ivar}])
                varid = netcdf.inqVarID(ncid,[varnames_osd{ivar},'_row_size']);
                data = netcdf.getVar(ncid,varid,'double');
                idx = find(data>0);
                dataset.osd.(varnames_osd{ivar}) = [dataset.osd.(varnames_osd{ivar});time(idx),lon(idx),lat(idx)];
            end
        end
        netcdf.close(ncid);
    end

    % Process CTD
    for ifile = 1:numel(ctd_files)
        disp(['Working on file ',ctd_files(ifile).folder,'/',ctd_files(ifile).name])
        ncid = netcdf.open([ctd_files(ifile).folder,'/',ctd_files(ifile).name],'NOWRITE'); % Open Netcdf file
        finfo = ncinfo([ctd_files(ifile).folder,'/',ctd_files(ifile).name]); % Get information on the file
        lon = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'),'double');   % degrees east
        lat = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'),'double');   % degrees north
        time = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'),'double'); % days since 1770-01-01 00:00:00 UTC
        time = time + datenum(1770,1,1);
        % Get data 
        for ivar = 1:numel(varnames_ctd)
            if ifile==1
                dataset.ctd.(varnames_ctd{ivar}) = [];
            end
            if ismember({[varnames_ctd{ivar},'_obs']},{finfo.Dimensions.Name})
                disp(['   - ',varnames_ctd{ivar}])
                varid = netcdf.inqVarID(ncid,[varnames_ctd{ivar},'_row_size']);
                data = netcdf.getVar(ncid,varid,'double');
                idx = find(data>0);
                dataset.ctd.(varnames_ctd{ivar}) = [dataset.ctd.(varnames_ctd{ivar});time(idx),lon(idx),lat(idx)];
            end
        end
        netcdf.close(ncid);
    end

    save('WODstats_data','dataset','varnames','WODdir')
else
    load WODstats_data dataset
end

if 0
    figure
    Tvec = datevec(dataset.osd.Temperature(:,1));
    uv = unique(Tvec(:,1));
    n  = histc(Tvec(:,1),uv);
    plot(uv(2:end-1),n(2:end-1),'.-')
    title('Temperature - water samples')
    set(gca,'XLim',[1800,2020])

    figure
    Tvec = datevec(dataset.osd.Oxygen(:,1));
    uv = unique(Tvec(:,1));
    n  = histc(Tvec(:,1),uv);
    plot(uv(2:end-1),n(2:end-1),'.-')
    title('Oxygen - water samples')
    set(gca,'XLim',[1800,2020])

    figure
    Tvec = datevec(dataset.osd.Chlorophyll(:,1));
    uv = unique(Tvec(:,1));
    n  = histc(Tvec(:,1),uv);
    plot(uv(2:end-1),n(2:end-1),'.-')
    title('Chlorophyll - water samples')
    set(gca,'XLim',[1800,2020])

    figure
    Tvec = datevec(dataset.osd.Nitrate(:,1));
    uv = unique(Tvec(:,1));
    n  = histc(Tvec(:,1),uv);
    plot(uv(2:end-1),n(2:end-1),'.-')
    title('Nitrate - water samples')
    set(gca,'XLim',[1800,2020])
end

if 0
    figure
    Tvec = datevec(dataset.osd.Temperature(:,1));
    uv = unique(Tvec(:,1));
    n  = histc(Tvec(:,1),uv);
    plot(uv(2:end-1),n(2:end-1),'.-')
    hold on

    Tvec = datevec(dataset.osd.Oxygen(:,1));
    uv = unique(Tvec(:,1));
    n  = histc(Tvec(:,1),uv);
    plot(uv(2:end-1),n(2:end-1),'.-')

    Tvec = datevec(dataset.osd.Chlorophyll(:,1));
    uv = unique(Tvec(:,1));
    n  = histc(Tvec(:,1),uv);
    plot(uv(2:end-1),n(2:end-1),'.-')

    Tvec = datevec(dataset.osd.Nitrate(:,1));
    uv = unique(Tvec(:,1));
    n  = histc(Tvec(:,1),uv);
    plot(uv(2:end-1),n(2:end-1),'.-')
    set(gca,'XLim',[1850,2020])
    ylabel('# stations')
    lg = legend('Temperature','Oxygen','Chlorophyll','Nitrate');
    set(lg,'Location','NorthWest')
    title('Water samples (#Stations)')
end

if 0
    figure
    Tvec = datevec(dataset.ctd.Temperature(:,1));
    uv = unique(Tvec(:,1));
    n  = histc(Tvec(:,1),uv);
    plot(uv(2:end-1),n(2:end-1),'.-')
    hold on

    Tvec = datevec(dataset.ctd.Oxygen(:,1));
    uv = unique(Tvec(:,1));
    n  = histc(Tvec(:,1),uv);
    plot(uv(2:end-1),n(2:end-1),'.-')

    Tvec = datevec(dataset.ctd.Chlorophyll(:,1));
    uv = unique(Tvec(:,1));
    n  = histc(Tvec(:,1),uv);
    plot(uv(2:end-1),n(2:end-1),'.-')

%     Tvec = datevec(dataset.osd.Nitrate(:,1));
%     uv = unique(Tvec(:,1));
%     n  = histc(Tvec(:,1),uv);
%     plot(uv(2:end-1),n(2:end-1),'.-')
%     set(gca,'XLim',[1800,2020])

    set(gca,'XLim',[1960,2020])
    ylabel('# casts')
    lg = legend('Temperature','Oxygen','Chlorophyll');
    set(lg,'Location','NorthWest')
    title('CTD casts')
end
