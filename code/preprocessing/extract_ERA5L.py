import numpy as np
import xarray as xr
import dask
import pandas as pd
import time
import argparse
import os
import datetime

t1 = time.time()
parser = argparse.ArgumentParser()
parser.add_argument('--rgiid',type=str,required = True,help = 'input RGIId')
parser.add_argument('--datemin',type=str,required = True,help = 'input start date')
parser.add_argument('--datemax',type=str,required = True,help = 'input end date')
args = parser.parse_args()
mindate = args.datemin
maxdate = args.datemax
Id = args.rgiid
print(Id)

#parse dates into datetime objects
mindate_split = np.array(mindate.split(':')).astype(int)
mindate_dt = datetime.datetime(year=mindate_split[0],month=mindate_split[1],day = mindate_split[2],
                               hour=mindate_split[3],minute=mindate_split[4],second=mindate_split[5])
mindate_dt = mindate_dt - datetime.timedelta(hours=1)  #required for accummulated variables
maxdate_split = np.array(maxdate.split(':')).astype(int)

maxdate_dt = datetime.datetime(year=maxdate_split[0],month=maxdate_split[1],day = maxdate_split[2],
                               hour=maxdate_split[3],minute=maxdate_split[4],second=maxdate_split[5])

minyr = int((mindate_dt-datetime.timedelta(hours=12)).year)  #-1 day to make sure to cover for if timezone gmt +12
maxyr = int((maxdate_dt+datetime.timedelta(hours=12)).year)

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


years = np.arange(minyr,maxyr+1,1)
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
    ds = xr.open_mfdataset( data_paths,chunks={'time': -1, 'latitude':1, 'longitude':1},parallel = False)  #the chunking is extremely important, 
                                                                                                            #reduce overhead by chunking to single points
                                                                                                            # so in next step, only one chunk needs to be loaded to extract the data

    
    for lon,lat in zip(lons,lats):

        deltaGMT = np.ceil(lon / 15.0)
        mindate_dt_GMT = mindate_dt - datetime.timedelta(hours = deltaGMT)
        maxdate_dt_GMT = maxdate_dt - datetime.timedelta(hours = deltaGMT)
        point_pos = [lat,lon]  #lat,lon
        h_ERA5L = height.sel(latitude=point_pos[0], longitude=point_pos[1], method='nearest').z.values
        h_dem = points_df.elev_m.iloc[index-1]
        delta_T = ((h_dem-h_ERA5L)/1000)*lapse_rate
    
        
        
        
        out_ds = ds.sel(latitude=point_pos[0], longitude=point_pos[1], method='nearest').sel(time=slice(mindate_dt_GMT,maxdate_dt_GMT)).compute()
        df_out = out_ds.to_dataframe().reset_index()

        df_out['sp'] = barometric_p(df_out['sp'],h_dem-h_ERA5L,df_out.t2m)
        df_out['Ws'] = (df_out['u10']**2 + df_out['v10']**2)**(0.5)
        
        df_out['t2m'] = df_out['t2m']+delta_T
        df_out['d2m'] = df_out['d2m']+delta_T
        
        #############deal with accummulated variables (https://confluence.ecmwf.int/pages/viewpage.action?pageId=197702790)
        df_out['strd_mod'] = df_out['strd'].diff()/3600
        mask = df_out.time.dt.hour == 1
        df_out.loc[mask, 'strd_mod'] = df_out.loc[mask, 'strd'] / 3600
        df_out.strd = df_out.strd_mod

        df_out['ssrd_mod'] = df_out['ssrd'].diff()/3600
        mask = df_out.time.dt.hour == 1
        df_out.loc[mask, 'ssrd_mod'] = df_out.loc[mask, 'ssrd'] / 3600
        df_out.ssrd = df_out.ssrd_mod

        df_out['tp_mod'] = df_out['tp'].diff()*1000
        mask = df_out.time.dt.hour == 1
        df_out.loc[mask, 'tp_mod'] = df_out.loc[mask, 'tp'] *1000
        df_out.tp = df_out.tp_mod

        df_out = df_out.drop(['u10','v10','longitude','latitude','ssrd_mod','strd_mod','tp_mod'],axis=1)
        df_out.time = df_out.time + datetime.timedelta(hours=deltaGMT)
        df_out = df_out.iloc[1:]    # drop first row, has Nans for accummulated variables
        df_out.to_parquet(f'{root}/All_glaciers/{Id}/ERA5L_at_gridpoint/{index}.parquet')
        index+=1 
        #print(index)
    print(time.time()-t1)
    print('ERA5 extraction complete')
else:
    print('all_files_exist')
