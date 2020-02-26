
ncgrd = netcdf.open('/misc/1/input/laurenta/NEMO/CREG025E_grid.nc','NC_NOWRITE');
mask_t = netcdf.getVar(ncgrd,netcdf.inqVarID(ncgrd,'mask_t'),'double');
nav_lon = netcdf.getVar(ncgrd,netcdf.inqVarID(ncgrd,'nav_lon'),'double');
nav_lat = netcdf.getVar(ncgrd,netcdf.inqVarID(ncgrd,'nav_lat'),'double');
bathy = netcdf.getVar(ncgrd,netcdf.inqVarID(ncgrd,'bathy'),'double');
netcdf.close(ncgrd);

mask = mask_t(:,:,1);mask(mask==0) = nan;

latlim = [-5 5]; 
lonlim = [-170 -120];

% Need to loop over time
t = [datenum('may 1, 2009') datenum('may 31, 2009')]; 
f = argofiles(latlim,lonlim,t,'pacific');

[lat,lon,t,P,T,S,pn] = argodata(f);

% Need to loop over the profiles


plot(lon,lat,'ko','markerfacecolor','y') 
xlabel 'longitude'
ylabel 'latitude'