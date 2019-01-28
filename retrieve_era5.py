#!/usr/bin/env python
"""
Save as get-tp.py, then run "python get-tp.py".
  
Input file : None
Output file: tp_20170101-20170102.nc
"""
import cdsapi
import sys
from datetime import datetime, timedelta
from calendar import monthrange


c = cdsapi.Client()

if len(sys.argv) < 3:
    print(' not enough parameters supplied!')
    sys.exit(1)

strArea = sys.argv[1] # N/W/S/E
date    = sys.argv[2] # %Y-%m-%d/to/%Y-%m-%d
outfile = sys.argv[3] 

date_array	= date.split("/")
date_start	= datetime.strptime(date_array[0],"%Y-%m-%d")
date_end	= datetime.strptime(date_array[2],"%Y-%m-%d")

date_start_new 	= date_start
date_end_new	= date_end + timedelta(days=1)

date_start_string 	= str(date_start_new.year)+"-"+str(date_start_new.month).zfill(2)+"-"+str(date_start_new.day).zfill(2)
date_end_string 	= str(date_end_new.year)+"-"+str(date_end_new.month).zfill(2)+"-"+str(date_end_new.day).zfill(2)
date_string = date_start_string+"/to/"+date_end_string

print(' era5 reanalysis retrieval running...')
print('  region   : '+strArea)
print('  dates    : '+date_string)
print('  outfile  : '+outfile)

print(date_start_new.month)

months = list(range(date_start_new.month,date_end_new.month+1))
years  = list(range(date_start_new.year,date_end_new.year+1))
lvls   = list(range(0,138))

outnamesfc = '{:s}_sfc.nc'.format(outfile)
outnameatm = '{:s}_atm.nc'.format(outfile)

# ----------------------------------------------------------------
# parameter ids
# ----------------------------------------------------------------
# https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
#
# 129 ... geopotential
# 130 ... temperature
# 131 ... U component of wind
# 132 ... V component of wind
# 133 ... specific humidity
# 135 ... vertical velocity
# 155 ... divergence
# 52  ... ?
# 156 ... ?
# 246 ... specific cloud liquid water content
# 247 ... specific cloud ice water content


# query atmospheric data
'''
r = c.retrieve('reanalysis-era5-complete', {
    'class'   : 'ea',
    'expver'  : '1',
    'stream'  : 'oper',
    'type'    : 'an',
    'param'   : '129.128/130.128/131.128/132.128/133.128/135.128/155.128/246.128/247.128',
    'levtype' : 'ml',
    'levelist': '30/to/137',  # basically query all levels below 30 km
    'date'    : date_string,
    'area'    : strArea,
    'grid'    : '0.25/0.25',
    'time'    : '00/to/12',
    'format'  : 'netcdf'
})

r.download(outnameatm)
'''
'''
c.retrieve("reanalysis-era5-complete", {
        "product_type": "reanalysis",
        'area'        : strArea,
        "model_level" : lvls,
        'grid'        : '0.25/0.25',
        "variable":       ["129.128/130.128/131.128/132.128/133.128/135.128/155.128/52.162/156.128/246.128/247.128"],
        'year'        : years,
        'month'       : months,
        'day'         : ["01","02","03","04","05","06","07","08","09","10","11",
                       "12","13","14","15","16","17","18","19","20","21","22",
                       "23","24","25","26","27","28","29","30","31"],
        'time'        : [
            '00:00','01:00','02:00',
            '03:00','04:00','05:00',
            '06:00','07:00','08:00',
            '09:00','10:00','11:00',
            '12:00','13:00','14:00',
            '15:00','16:00','17:00',
            '18:00','19:00','20:00',
            '21:00','22:00','23:00'
        ],
        'format':         'netcdf'
    })
'''
#r.download(outnameatm)


# query surface data
r = c.retrieve(
    'reanalysis-era5-single-levels', {
            'variable'    : ['129.128/134.128'],          # query surface pressure
            'area'        : strArea,
            'grid'        : '0.25/0.25',
            'product_type': 'reanalysis',
            'year'        : years,
            'month'       : months,
            'day'         : ["01","02","03","04","05","06","07","08","09","10","11",
                           "12","13","14","15","16","17","18","19","20","21","22",
                           "23","24","25","26","27","28","29","30","31"],
            'time'        : [
                '00:00','01:00','02:00',
                '03:00','04:00','05:00',
                '06:00','07:00','08:00',
                '09:00','10:00','11:00',
                '12:00','13:00','14:00',
                '15:00','16:00','17:00',
                '18:00','19:00','20:00',
                '21:00','22:00','23:00'
            ],
            'format'      : 'netcdf'
    })
r.download(outnamesfc)

