import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

import os
import sys
print("Current Conda environment:", os.environ.get("CONDA_DEFAULT_ENV"))

import cartopy.crs as ccrs
from tqdm import tqdm
from apexpy import Apex
import datetime as dt
import importlib as il

sys.path.append('/home/pxv220016/prasoon/data/sat_interp_repo/repo2/prasoon_utility_programs')
sys.path.append('/home/pxv220016/prasoon/data/sat_interp_repo/repo2/Qingyu_Cesar_EIA_IHA')
import functions


from p_tqdm import p_map
from multiprocessing import Pool
import itertools




work = '/home/pxv220016/prasoon/data/sat_interp_repo/repo2/'
scratch = '/home/pxv220016/scratch/'

functions = il.reload(functions)
start_index = int(sys.argv[1])
end_index = int(sys.argv[2])
zone = str(sys.argv[3])        # '75W', '50W' - GLONS sectors
month = str(sys.argv[4])       # 'march','sept', 'dec', 'june' 


if month == 'sept' or month == 'march':
    phase = 'equinox'
else:
    phase = 'solstice'

print(month, zone, start_index, end_index)

'''
if month == 'dec':
    years_tot = list(range(2015, 2024))
else:
    years_tot = list(range(2012, 2025))
'''

#for year in tqdm(years_tot):
for year in tqdm(list(range(start_index, end_index+1))):
  
    print(year)

    glat_lim = 40

    if zone == '75W':
        glon_min = -85
        glon_max = -60
        mlon_min = -10
        mlon_max = 10
        zone_mlon = '0W'
    elif zone == '50W':
        glon_min = -65
        glon_max = -40
        mlon_min = 15
        mlon_max = 25
        zone_mlon = '20E'
        
    
    # Reading Madrigal Cedar data for +- 21 days around March equinox of 2010-2024
    # Files in below scratch folder are obtained by using multiple file download 
    # command in ASCII format and then doing `gunzip file.gz`  
    path = f'{scratch}Qingyu_Cesar_EIA/{month}_data/{str(year)}_{month}_{phase}/'
    files = os.listdir(path)
    files = [path + i for i in files]

    for f in files:
        if f[-4:] != '.txt':
            files.remove(f)
            
    #print('No. of files to process:', len(files))
    
    tec_g = []
    columns = ['GDLAT', 'GLON', 'TEC', 'DT']
    grnd_tec = pd.DataFrame(columns=columns)
    
    def process_file(f):
        # Read function passed continuously during multiple processing to quicken the process
        df = pd.read_csv(f, sep=r'\s+')
        d = [dt.datetime(y, m, d, h, mi, s) for y, m, d, h, mi, s in zip(df.YEAR, df.MONTH, df.DAY, df.HOUR, df.MIN, df.SEC)]
        df['DT'] = d
        # Dropping unnecassary columns from the Dataframe
        df = df.drop(['RECNO', 'KINDAT', 'KINST', 'UT1_UNIX', 'UT2_UNIX', 'YEAR', 'MONTH', 'DAY', 'HOUR', 'MIN', 'SEC', 'DTEC'], axis=1)
        if 'GDALT' in df.columns:
            df = df.drop(['GDALT'], axis=1)

        df = df[(df.GDLAT > -glat_lim) & (df.GDLAT < glat_lim) & (df.GLON > glon_min) & (df.GLON < glon_max)].reset_index(drop=True)
        return df    

    
    # Speeding the process by using parallel processing
    tec_g = p_map(process_file, files)  # Parallel processing with progress bar
    print('1')
    grnd_tec = pd.concat(tec_g, axis=0).reset_index(drop=True)
    grnd_tec0 = grnd_tec.sort_values(by=['DT', 'GDLAT'], ascending=[True, True])
    #print('grnd_tec0\n', grnd_tec0)
    

    # Reading the Kp index values for all the days and filtering undesired points where Kp > 3
    file = f'{work}Qingyu_Cesar_EIA_IHA/kp3_index_values/kp_{str(year)}_{month}.txt'
    kp = pd.read_csv(file,sep=r'\s+')
    date_kp = [functions.day_to_date(i, year) for i in kp.DOY]
    m, d = zip(*date_kp)
    kp['DT'] = [dt.datetime(year, j, i, k, 0, 0) for i,j,k in zip(d,m,kp.Hour)]
    kp['kp'] = [i/10 for i in kp.Kp]
    kp = kp.drop(['Year', 'Hour','Kp'], axis = 1)
    # the below command is messing up the kp values as it is trying to fit the kp value 
    # based on weighted average of all the known points surrounding it
    grnd_tec0 = pd.merge_asof(grnd_tec0.sort_values('DT'), kp.sort_values('DT'), on='DT', direction='nearest')
    #grnd_tec0 = merged_data.drop(['date'], axis = 1)
    print('2')
    
    
    # Calculation of magnetic coordinates by using Apex library and Parallel prcoessing 
    t_start = dt.datetime.now() # just a timer
    i = 1
    with Pool(40) as pool:
        print(i)
        i += 1
        p = pool.starmap(functions.magnetic_coords_parallel, zip(grnd_tec0.DT, grnd_tec0.GDLAT, grnd_tec0.GLON, grnd_tec0.TEC))
    pool.close()
    pool.join()
    # Separating the data from output list 
    sat_date, sat_glat, sat_glon, sat_tec, sat_mlat, sat_mlon, sat_mlt = zip(*p)


    # Reordering the outputs and applying further conditions on magnetic coordinates
    grnd_temp = pd.DataFrame({'DT': sat_date, 'DOY': grnd_tec0.DOY, 'GDLAT': sat_glat, 'GLON': sat_glon, 'TEC': sat_tec, 'MLAT': sat_mlat, 'MLON': sat_mlon, 'MLT': sat_mlt, 'kp':grnd_tec0.kp, 'F10.7':grnd_tec0['F10.7']})
    
    
    grnd_tec1 = grnd_temp.sort_values(by=['DT', 'GDLAT'], ascending=[True, True]).reset_index()
    grnd_tec1 = grnd_tec1[(grnd_tec1.MLON <= mlon_max) & (grnd_tec1.MLON >= mlon_min)].reset_index(drop=True)
    grnd_tec1 = grnd_tec1[(grnd_tec1.MLAT <= 40) & (grnd_tec1.MLAT >= -40)].reset_index(drop=True)
    grnd_tec2 = grnd_tec1.drop(['GDLAT', 'GLON','MLON'], axis = 1)
    print('grnd_tec2\n',grnd_tec2)

    
    # Writing the output into csv files for easy post processing
    grnd_tec2.to_csv(f'{scratch}Qingyu_Cesar_EIA/{zone}/outputs/{month}/{year}_{month}_{phase}.csv', index=False)

    t_total = dt.datetime.now() - t_start
    print(t_total)
    
