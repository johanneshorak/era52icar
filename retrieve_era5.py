#!/usr/bin/env python
import cdsapi
import sys
import pandas as pd
from datetime import datetime, timedelta
from calendar import monthrange
import os.path
import time

c = cdsapi.Client()

def hms_string(sec_elapsed):
    h = int(sec_elapsed / (60 * 60))
    m = int((sec_elapsed % (60 * 60)) / 60)
    s = sec_elapsed % 60.
    return "{}:{:>02}:{:>05.2f}".format(h, m, s)

if len(sys.argv) < 3:
    print(' not enough parameters supplied!')
    sys.exit(1)

strArea = sys.argv[1] # N/W/S/E
date    = sys.argv[2] # %Y-%m-%d/to/%Y-%m-%d
outfile = sys.argv[3] 

date_array    = date.split('/')
date_start    = datetime.strptime(date_array[0],'%Y-%m-%d')
date_end    = datetime.strptime(date_array[2],'%Y-%m-%d')

date_start_new     = date_start
date_end_new    = date_end# + timedelta(days=1)

date_start_string     = str(date_start_new.year)+'-'+str(date_start_new.month).zfill(2)+'-'+str(date_start_new.day).zfill(2)
date_end_string     = str(date_end_new.year)+'-'+str(date_end_new.month).zfill(2)+'-'+str(date_end_new.day).zfill(2)
date_string         = date_start_string+'/to/'+date_end_string

print(' era5 reanalysis retrieval running...')
print('  region   : '+strArea)
print('  dates    : '+date_string)
print('  outfile  : '+outfile)

# for the surface variable sp we construct the request differently.
#
# this assumes that we query something where we need to download 3 months.
# e.g. 
#      26th of february to 28th of february 2009
#      1st  of march    to 31st of march    2009
#      1st  of april    to  2nd of april    2009
#
# this needs to be broken up in three requests.
# we don't need the entire month since we're currently forcing ICAR
# like this (e.g. a couple of days spinup, then the month we're interested
# in, and then one or two days more for cases where we're converting
# ICAR output to other timezones.

day0   = date_start_new.day
day2   = date_end_new.day
month0 = date_start_new.month
month2 = date_end_new.month
year0  = date_start_new.year
year2  = date_end_new.year


if month0 == 12:
    month1 = 1
    year1  = year0+1
    sameyear = False
else:
    month1 = month0+1
    year1  = year0
    sameyear = True

print('')
# construct temporal boundaries of the request
if month0 != month2:
    requests = [
        [day0,monthrange(year0,month0)[1],month0, year0],
        [1,   monthrange(year1,month1)[1],month1, year1],
        [1,   day2,                       month2, year2]
    ]   
    print('  standard request including spinup/spindown time and one complete month')
elif month0 == month2:
    requests = [
        [day0,day2,month0, year0]
    ]
    print('  request for slice of a month')



if not sameyear and not((month0 == 12) and (month1 == 1)):
    print('  request construction for multiple years not yet implemented!')
    sys.exit(1)

for nr in range(len(requests)):
    print('    ',requests[nr])

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

print('---------------------------------------------------------------')
print('  running request for ERA5 atmospheric data')
print('---------------------------------------------------------------')


atm_parameters = ['130.128','131.128','132.128','133.128','135.128','155.128','246.128','247.128']

n = 0

t0_total = time.time()

while n < len(atm_parameters):
    t0    = time.time()                       # clock the required wall-time of the request
    param = atm_parameters[n]
    
    print('   querying {:s}'.format(param))
    
    outnameatm = '{:s}{:s}{:s}-{:s}{:s}{:s}_{:s}_{:s}_atm.nc'.format(str(year0),str(month0).zfill(2),str(day0).zfill(2),str(year2),str(month2).zfill(2),str(day2).zfill(2),param,outfile)

    r = c.retrieve('reanalysis-era5-complete', {
        'class'   : 'ea',
        'expver'  : '1',
        'stream'  : 'oper',
        'type'    : 'an',
        'param'   : param,
        'levtype' : 'ml',
        'levelist': '1/to/137',  # basically query all levels below 30 km (starting from 30)
        'date'    : date_string,
        'area'    : strArea,
        'grid'    : '0.25/0.25',
        'time'    : '00/01/02/03/04/05/06/07/08/09/10/11/12/13/14/15/16/17/18/19/20/21/22/23',
        'format'  : 'netcdf'
    })

    r.download(outnameatm)
    
    # test whether the file downloaded. if not redo the current parameter
    if not os.path.isfile(outnameatm):
        print('   error querying {:s}, retrying'.format(param))
    else:
        t1 = time.time()
        print('   wall-time: {:s}'.format(hms_string(t1-t0)))
        n+=1


# query surface data
print('---------------------------------------------------------------')
print('  running request for ERA5 surface data')
print('---------------------------------------------------------------')
nr = 0
while nr < len(requests):
    t0    = time.time()                       # clock the required wall-time of the request
    req   = requests[nr]
    
    print('   working on request ',req)
    day0  = req[0]
    day1  = req[1]
    month = req[2]
    year  = req[3]
    
    days  = list(map(str,range(day0,day1+1)))
    for ns,s in enumerate(days):
        days[ns] = s.zfill(2) # need days in 02, 03,.. etc formatting
        
    outnamesfc = '{:s}{:s}{:s}-{:s}{:s}{:s}_{:s}_sfc.nc'.format(str(year),str(month).zfill(2),str(day0).zfill(2),str(year),str(month).zfill(2),str(day1).zfill(2),outfile)
    
    r = c.retrieve('reanalysis-era5-single-levels', {
        'grid'        : '0.25/0.25',
        'product_type': 'reanalysis',
        'variable'    : ['orography','surface_pressure'],          # query geopotential height and surface pressure
        'area'    : strArea,
        'grid'    : '0.25/0.25',
        'year':     str(year),
        'month':    str(month).zfill(2),
        'day':      days,
        'time':     ['00:00','01:00','02:00','03:00','04:00','05:00',
                     '06:00','07:00','08:00','09:00','10:00','11:00',
                     '12:00','13:00','14:00','15:00','16:00','17:00',
                     '18:00','19:00','20:00','21:00','22:00','23:00'],
        'format'  : 'netcdf'
    })
    r.download(outnamesfc)

    # test whether the file downloaded. if not redo the current parameter
    if not os.path.isfile(outnamesfc):
        print('   error querying request {:n}, retrying'.format(nr))
    else:
        t1 = time.time()
        print('   wall-time: {:s}'.format(hms_string(t1-t0)))
        nr+=1

t1 = time.time()
print('total wall-time: {:s}'.format(hms_string(t1-t_total)))
