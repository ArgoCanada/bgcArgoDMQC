import bgcArgo as bgc

# setup for your system - these directories need to already exist!
argo_path = 'your/argo/data/path' # where to save Argo data
ncep_path = 'your/ncep/data/path' # where to save NCEP reanalysis data
woa_path  = 'your/woa18/data/path' # where to save WOA data

# download the data - this can take some time depending on connection
# Argo
wmos = [4902480, 6902905]
dacs = [bgc.get_dac(w) for w in wmos]
dacpath = '/ifremer/argo/dac'
fltpath = ['{}/{}/{}'.format(dacpath, d, w) for d, w in zip(dacs, wmos)]
# NCEP
bgc.io.get_ncep('pres', local_path=ncep_path)
bgc.io.get_ncep('land', local_path=ncep_path)
# WOA
bgc.io.get_woa18('O2sat', local_path=woa_path)

# tell the package where to look for data
bgc.set_dirs(argo_path=argo_path, ncep_path=ncep_path, woa_path=woa_path)
# load data from profiles for two floats
flts = bgc.profiles(wmos)
# calculate the dissolved oxygen gains
gains = flts.calc_gains()
print(gains)

# load a synthetic profile
syn = bgc.sprof(4902480)
# plot a time vs. depth section for the top 500m
g1 = syn.plot('cscatter', varname='DOXY', ylim=(0,500))
# plot the first 10 profiles for temperature, practical salinity,
# and adjusted oxygen
g2 = syn.plot('profiles', varlist=['TEMP','PSAL', 'DOXY'], Ncycle=10, Nprof=10, ylim=(0,500))