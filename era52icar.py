import xarray as xa
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import lib.gpcalc as gpcalc
import sys as sys

if len(sys.argv) < 3:
    print(' syntax: python era52icar.py [era5_file] [outname]')
    sys.exit(1)

e5_nc   = sys.argv[1]
outname = sys.argv[2]

try:
    era5_ds = xa.open_dataset('./{:s}'.format(e5_nc))
except:
    print(' error opening {:s}!'.format(e5_nc))
    sys.exit(1)

try:
    ml_df  = pd.read_csv('./data/model_level.csv',index_col='n')
except:
    print(' error opening model level data!')
    sys.exit(1)

#era5_ds = era5_ds.sel(level=slice(130,137)) # just for testing

gpcalc.set_data(era5_ds,ml_df,137)    # set the data required for calculations

# -
# create a dataset that contains the ak and bk coefficients of ERA5
# later needed to calculate pressure at each model level from surface pressure
# -
ab_ds = xa.Dataset(
    coords={
        'level'        : era5_ds.level
    },
    data_vars={
        'ak'     : (['level'],ml_df.loc[era5_ds.level.values,'a [Pa]'].values),
        'bk'     : (['level'],ml_df.loc[era5_ds.level.values,'b'].values),        
    }
)
era5_ds = era5_ds.merge(ab_ds)

Nt   = len(era5_ds.time)
Nlon = len(era5_ds.longitude)
Nlat = len(era5_ds.latitude)
Nlvl = len(era5_ds.level)

# first - calculate the pressure and geopotential (height) at every model level
# set pressure at level 0 (= top of the atmosphere) to nan. this shouldn't be a full level anymore
# and therefor we can't assign a pressure to it.

p           = np.zeros(Nt*Nlvl*Nlon*Nlat).reshape(Nt,Nlvl,Nlat,Nlon)  # pressure
ph          = np.zeros(Nt*Nlvl*Nlon*Nlat).reshape(Nt,Nlvl,Nlat,Nlon)  # geopotential height

print(' calculating pressure and geopotential at model levels...')
for n in range(0,Nlvl):
    lvl  = era5_ds.level.values[n]
    ak0  = era5_ds.ak.values[n]        # coef. for half level below (higher pressure)
    bk0  = era5_ds.bk.values[n]
    
    if n > 0:
        ak1  = era5_ds.ak.values[n-1]  # coef. for half level above (lower pressure)
        bk1  = era5_ds.bk.values[n-1]
    else:
        # for n = 0 we have to look up the ak and bk values from
        # the ml_df dataframe since we did not store those in the xarray
        # dataset. we do this only if lvl > 0, otherwise it doesn't matter
        # since the pressure at full-lvl = 0 cannot be specified anyways
        if lvl > 0:
            ak1  = ml_df.loc[lvl-1,'a [Pa]']
            bk1  = ml_df.loc[lvl-1,'b']
        
    
    ps  = era5_ds.sp[:].values
    print('  {:4n} {:4n} | {:10.4f} {:10.4f} | {:10.4f} {:10.4f} |'.format(n,lvl,ak0,bk0,ak1,bk1))
    
    ph[:,n] = gpcalc.get_phi(lvl)/9.81   # get_phi calculates geopotential at full model level
    
    if lvl == 0:
        p[:,n] = ps*np.nan
    else:
        p[:,n] = 0.5*((ak0+bk0*ps)+(ak1+bk1*ps))
    
era5_ds['p'] = (['time','level','latitude','longitude'],p)
era5_ds['ph'] = (['time','level','latitude','longitude'],ph)

print(' assigning data to forcing xarray dataset')
# prepare the data arrays
west_east = range(0,Nlon)
south_north = range(0,Nlon)
bottom_top  = (gpcalc.lvlmax-era5_ds.level)[::-1]                    # reverse order, 137 is lowest level in ERA5, is here now 0

xlong       = np.zeros(Nt*Nlon*Nlat).reshape(Nt,Nlat,Nlon)
xlat        = np.zeros(Nt*Nlon*Nlat).reshape(Nt,Nlat,Nlon)
Time        = era5_ds.time.values

HGT         = np.zeros(Nt*Nlon*Nlat).reshape(Nt,Nlat,Nlon)             # elevation at surface
PH          = np.zeros(Nt*Nlvl*Nlon*Nlat).reshape(Nt,Nlvl,Nlat,Nlon)   # geopotential height
U           = np.zeros(Nt*Nlvl*Nlon*Nlat).reshape(Nt,Nlvl,Nlat,Nlon)
V           = np.zeros(Nt*Nlvl*Nlon*Nlat).reshape(Nt,Nlvl,Nlat,Nlon)
QVAPOR      = np.zeros(Nt*Nlvl*Nlon*Nlat).reshape(Nt,Nlvl,Nlat,Nlon)
QCLOUD      = np.zeros(Nt*Nlvl*Nlon*Nlat).reshape(Nt,Nlvl,Nlat,Nlon)
QICE        = np.zeros(Nt*Nlvl*Nlon*Nlat).reshape(Nt,Nlvl,Nlat,Nlon)
QRAIN       = np.zeros(Nt*Nlvl*Nlon*Nlat).reshape(Nt,Nlvl,Nlat,Nlon)
QSNOW       = np.zeros(Nt*Nlvl*Nlon*Nlat).reshape(Nt,Nlvl,Nlat,Nlon)
PHB         = np.zeros(Nt*Nlvl*Nlon*Nlat).reshape(Nt,Nlvl,Nlat,Nlon)   # not used
PB          = np.zeros(Nt*Nlvl*Nlon*Nlat).reshape(Nt,Nlvl,Nlat,Nlon)   # not used
TSK         = np.zeros(Nt*Nlvl*Nlon*Nlat).reshape(Nt,Nlvl,Nlat,Nlon)   # not used
T           = np.zeros(Nt*Nlvl*Nlon*Nlat).reshape(Nt,Nlvl,Nlat,Nlon)   # potential temperature

# set values
U          = era5_ds.u[:,::-1,::-1,:]
V          = era5_ds.v[:,::-1,::-1,:]
P          = era5_ds.p[:,::-1,::-1,:]
QVAPOR     = era5_ds.q[:,::-1,::-1,:]/(1.0-era5_ds.q[:,::-1,::-1,:])
QCLOUD     = era5_ds.clwc[:,::-1,::-1,:]/(1.0-era5_ds.clwc[:,::-1,::-1,:])
QICE       = era5_ds.ciwc[:,::-1,::-1,:]/(1.0-era5_ds.ciwc[:,::-1,::-1,:])
QRAIN      = era5_ds.crwc[:,::-1,::-1,:]/(1.0-era5_ds.crwc[:,::-1,::-1,:])
QSNOW      = era5_ds.cswc[:,::-1,::-1,:]/(1.0-era5_ds.cswc[:,::-1,::-1,:])
PH         = era5_ds.ph[:,::-1,::-1,:]
HGT        = era5_ds.z[:,::-1,:]/9.81

T          = era5_ds.t[:,::-1,:]*((10.**5)/P)**(0.2854)

xlong[:,:] = era5_ds.longitude

for ny in range(Nlat):
    xlat[:,ny,:] = era5_ds.latitude[::-1].values[ny]
    
frc_ds = xa.Dataset(
    coords={
        'Time'        : Time,
      #  'west_east'   : west_east,
      #  'south_north' : south_north
    },
    data_vars={
        'XLONG'    : (['Time','south_north','west_east'],xlong),
        'XLAT'     : (['Time','south_north','west_east'],xlat),
        'HGT'      : (['Time','south_north','west_east'],HGT),
        'U'        : (['Time','bottom_top','south_north','west_east'],U),
        'V'        : (['Time','bottom_top','south_north','west_east'],V),
        'P'        : (['Time','bottom_top','south_north','west_east'],P),
        'PH'       : (['Time','bottom_top','south_north','west_east'],PH),
        'QVAPOR'   : (['Time','bottom_top','south_north','west_east'],QVAPOR),
        'QCLOUD'   : (['Time','bottom_top','south_north','west_east'],QCLOUD),
        'QICE'     : (['Time','bottom_top','south_north','west_east'],QICE),
        'QRAIN'    : (['Time','bottom_top','south_north','west_east'],QRAIN),
        'QSNOW'    : (['Time','bottom_top','south_north','west_east'],QSNOW),
        'PB'       : (['Time','bottom_top','south_north','west_east'],PB),
        'PHB'      : (['Time','bottom_top','south_north','west_east'],PHB),
        'TSK'      : (['Time','bottom_top','south_north','west_east'],TSK),
        'T'       : (['Time','bottom_top','south_north','west_east'],T)
    }
)

# copy the attributes of variables that have a correspondence in the era5 dataset

varmap=[
    ['u','U'],
    ['v','V'],
    ['q','QVAPOR'],
    ['clwc','QCLOUD'],
    ['ciwc','QICE'],
    ['crwc','QRAIN'],
    ['cswc','QSNOW']
]

frc_ds['P'].attrs['units']         = 'Pa'
frc_ds['P'].attrs['long_name']     = 'pressure'
frc_ds['P'].attrs['standard_name'] = 'pressure'

frc_ds['HGT'].attrs['units']         = 'm'
frc_ds['HGT'].attrs['long_name']     = 'geopotential height of orography surface'
frc_ds['HGT'].attrs['standard_name'] = 'surface_geopotential_height'

frc_ds['PH'].attrs['units']         = 'm'
frc_ds['PH'].attrs['long_name']     = 'geopotential height of grid cell'
frc_ds['PH'].attrs['standard_name'] = 'geopotential_height'

frc_ds['T'].attrs['units']         = 'K'
frc_ds['T'].attrs['long_name']     = 'potential temperature'
frc_ds['T'].attrs['standard_name'] = 'potential_temperature'

frc_ds['TSK'].attrs['units']         = ''
frc_ds['TSK'].attrs['long_name']     = 'unused variable'
frc_ds['TSK'].attrs['standard_name'] = ''

frc_ds['PB'].attrs['units']         = ''
frc_ds['PB'].attrs['long_name']     = 'unused variable'
frc_ds['PB'].attrs['standard_name'] = ''

frc_ds['PHB'].attrs['units']         = ''
frc_ds['PHB'].attrs['long_name']     = 'unused variable'
frc_ds['PHB'].attrs['standard_name'] = ''

frc_ds['Time'].encoding['units']         = 'hours since 1900-01-01 00:00:0.0'
frc_ds['Time'].encoding['dtype']         = 'f4'
frc_ds['Time'].encoding['calendar']      = 'gregorian'

for n in range(len(varmap)):
    row = varmap[n]
    eravar  = row[0]
    frcvar  = row[1]
    for key in era5_ds[eravar].attrs:
        val = era5_ds[eravar].attrs[key]
        frc_ds[frcvar].attrs[key] = val

comp     = dict(zlib=True, complevel=5)
encoding = {var: comp for var in frc_ds.data_vars}

print(' writing to disk and compressing...')
frc_ds.to_netcdf('./{:s}'.format(outname),encoding=encoding)
print(' finished')


