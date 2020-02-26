clearvars;close all

% The following command was used to get transect_xy
% [gridout,transect_xy] = transect_gridpoints_nemo_v2;
% Lines 216, 217 and 220 of the script were commented

if 0
% Limits for North Atlantic (3Oceans lower limit) including Labrador Sea and Baffin Bay
zone = 'NAtl';
transect_xy = [

   376     1
   336    12
   348    29
   353    46
   353    54
   361    73
   379    91
   389   120
   404   155
   404   163
   396   180
   389   189
   390   209
   399   246
   401   276
   399   298
   353   325
   335   339
   292   342
   241   332
   202   334
   194   339
   163   343
   148   339
   134   324
   144   283
   142   264
   127   248
   124   236
   136   219
   134   206
   126   177
   137   152
   132   138
   103   139
    85   114
    63    93
    50    74
    38    47
    38    27
    37     9];
end

if 1
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
end

load /misc/1/home/laurenta/NEMO/workspace/matlab/mat/mask_nemo3oceans
for i0 = 1:numel(transect_xy(:,1))
    latp(i0,1) = lat(transect_xy(i0,1),transect_xy(i0,2));
    lonp(i0,1) = lon(transect_xy(i0,1),transect_xy(i0,2));
end

date_start = datenum(2017,1,1);
date_stop = datenum(2019,7,31);

t = date_start:date_stop;
count = nan(numel(t),1);
temp = zeros(numel(t),1);
do = zeros(numel(t),1);
no3 = zeros(numel(t),1);
ph = zeros(numel(t),1);
chl = zeros(numel(t),1);
bbp = zeros(numel(t),1);
for itime = 1:numel(t)
    disp(datestr(t(itime)))
    [ncfiles,lat,lon,tmp,floatid] = argofiles(latp,lonp,t(itime),'atlantic');
    if ~isempty(floatid)
        count(itime) = numel(unique(floatid));
        fnames = unique(ncfiles);
        for ifloat = 1:numel(fnames)
            rng('shuffle');
            tmpfile = ['argotmp',num2str(randi([1000 9999],1,1),'%d'),'.nc']; 
            % Download data file
            try
                urlwrite(fnames{ifloat},tmpfile);
            catch
                urlwrite(strrep(fnames{ifloat},'nodc_R','nodc_D'),tmpfile);
            end
            nc = netcdf.open(tmpfile,'NC_NOWRITE');
            try
                netcdf.inqVarID(nc,'TEMP')
                temp(itime) = temp(itime)+1;
            end
            try
                netcdf.inqVarID(nc,'DOXY')
                do(itime) = do(itime)+1;
            end
            try
                netcdf.inqVarID(nc,'NITRATE')
                no3(itime) = no3(itime)+1;
            end
            try
                netcdf.inqVarID(nc,'PH_IN_SITU_TOTAL')
                ph(itime) = ph(itime)+1;
            end
            try
                netcdf.inqVarID(nc,'CHLA')
                chl(itime) = chl(itime)+1;
            end
            try
                netcdf.inqVarID(nc,'BBP700')
                bbp(itime) = bbp(itime)+1;
            end
            netcdf.close(nc)
            delete(tmpfile) 
        end
    end
end

save(['ArgoCount_',zone,'_',datestr(date_start,'yyyymmdd'),'_',datestr(date_stop,'yyyymmdd')],'count','t','do','no3','ph','chl','bbp')
