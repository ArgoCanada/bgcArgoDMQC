clearvars;

fglobal = 'arnaud/ar_index_global_prof.txt';
fbio = 'arnaud/argo_bio-profile_index.txt';

date_start = datenum(1997,1,1);
date_stop = datenum(2019,9,1);

vtypes = {'DOXY','NITRATE','PH_IN_SITU_TOTAL','CHLA','BBP700','DOWNWELLING_PAR','CDOM'};
vnames = {'OXYGEN','NITRATE','PH','CHLA','BBP700','DOWNWELLING_PAR','CDOM'};
vnames_low = {'do','no3','ph','chl','bbp','par','cdom'};

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

% Limits for North West North Atlantic including Labrador Sea and southern part of Baffin Bay
zone = 'NWAtl';
NWA_ll = [ 30   70   70   66.7053   63.3478   62.0954   62.1360   61.7930   61.6015   61.6404   61.4073   60.4241   54.6531   49.7625   32   30   30;...
          -35  -35  -75.50  -68.8382  -69.9071  -66.3152  -66.0641  -65.8167  -65.4046  -65.1520  -64.9955  -64.6588  -62.2532  -72.7457  -83  -83  -35];

[inpoly,onpoly] = inpolygon(lona_bio,lata_bio,NWA_ll(2,:),NWA_ll(1,:));
inpoly = inpoly+onpoly;

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
