import numpy as np
import xarray as xr
import dask
import pandas as pd
import time
import argparse
import os


t1 = time.time()
parser = argparse.ArgumentParser()
parser.add_argument('--rgiid',type=str,required = True,help = 'input RGIId')
args = parser.parse_args()
Id = args.rgiid
print(Id)



root = '/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing'
glaciers = pd.read_csv(f'{root}/Glacier_list.csv',header = 2)

glaciers = glaciers[glaciers['RGI index'] == Id]
height = xr.load_dataset(f'{root}/ERA5L/geopotential.nc')


g = 9.80665
M = 0.02896968
R = 8.314462618



height['z'] = height.z/9.80665
lapse_rate = -6.5 #deg/km


def barometric_p(P,dh,T):
    a = -lapse_rate*dh/(1000*T)    #note lapse rate in atmospheric is +ve, hence change sign in all formulas
    b = (g*M)/(R*lapse_rate/1000)

    return P*(1+a)**b


years = [2019,2020]
vars = ['Ta','DTa','Lwin','Precip','Psfc','Swin','u10','v10']



    

points_df = pd.read_csv(f'{root}/All_glaciers/{Id}/coords_out_{Id}.csv', header = 1)
#points_df = pd.read_csv(f'All_glaciers/RGI60-15.07886/coords_out_RGI60-15.07886.csv', header = 1)

lats = points_df['lat']
lons = points_df['lon']
index = 1
data_paths =[]


num_points = len(points_df)
expected_paths = []
for i in range (num_points):
     expected_paths.append(f'{root}/All_glaciers/{Id}/ERA5L_at_gridpoint/{i}.parquet')
print('number of points EXTRAct')
print(num_points)
if not all(os.path.isfile(p) for p in expected_paths):


    for year in years:
            for var in vars:
                data_paths.append(f'{root}/ERA5L/ERA5Land_HMA_{var}_{year}.nc')
    ds = xr.open_mfdataset( data_paths,chunks={'time': 17544, 'latitude':1, 'longitude':1})

    
    for lon,lat in zip(lons,lats):

        
        point_pos = [lat,lon]  #lat,lon
        h_ERA5L = height.sel(latitude=point_pos[0], longitude=point_pos[1], method='nearest').z.values
        h_dem = points_df.elev_m.iloc[index-1]
        delta_T = ((h_dem-h_ERA5L)/1000)*lapse_rate
    
        
        
        
        out_ds = ds.sel(latitude=point_pos[0], longitude=point_pos[1], method='nearest')
        df_out = out_ds.to_dataframe().reset_index()

        df_out['sp'] = barometric_p(df_out['sp'],h_dem-h_ERA5L,df_out.t2m)
        df_out['Ws'] = (df_out['u10']**2 + df_out['v10']**2)**(0.5)
        
        df_out['t2m'] = df_out['t2m']+delta_T
        df_out['d2m'] = df_out['d2m']+delta_T
        df_out = df_out.drop(['u10','v10','longitude','latitude'],axis=1)


        df_out.to_parquet(f'{root}/All_glaciers/{Id}/ERA5L_at_gridpoint/{index}.parquet')
        index+=1 
        #print(index)
    print(time.time()-t1)
    print('ERA5 extraction complete')
else:
    print('all_files_exist')
