'''
Extract ERA5L timeseries at glacier points
adjust P,T to elevation of glacier point

INPUTS
must define rgiid,datemin, datemax when launching script
need ERA5L data grouped in yearly folders for vars ['Ta','DTa','Lwin','Precip','Psfc','Swin','u10','v10'] 
                note: Precip,Swin,Lwin are all accummulated variables
need geopotential height for region (ERA5L)
OUTPUTS


TO DO
use a local monthly averaged hourly lapse rate!
Maybe similar with a ERA5L native tp lapse rate? 
'''


import numpy as np
import xarray as xr
import pandas as pd
import time
import argparse
import os
import datetime

t1 = time.time()
#### parse inputs
parser = argparse.ArgumentParser()
parser.add_argument('--rgiid',type=str,required = True,help = 'input RGIId')
parser.add_argument('--datemin',type=str,required = True,help = 'input start date')  #dates in UTC!!!!
parser.add_argument('--datemax',type=str,required = True,help = 'input end date')
parser.add_argument('--outlocation',type=str,required = False,help = 'outlocation')
args = parser.parse_args()
mindate = args.datemin
maxdate = args.datemax
Id = args.rgiid
outlocation=args.outlocation


input_path='/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing'
if outlocation:
    root=f'{outlocation}'
else:
    root = '/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing'

#parse dates into datetime objects
mindate_split = np.array(mindate.split(':')).astype(int)
mindate_dt = datetime.datetime(year=mindate_split[0],month=mindate_split[1],day = mindate_split[2],
                               hour=mindate_split[3],minute=mindate_split[4],second=mindate_split[5])
mindate_dt = mindate_dt - datetime.timedelta(hours=1)  #required for accummulated variables (first row dropped at end)
maxdate_split = np.array(maxdate.split(':')).astype(int)

maxdate_dt = datetime.datetime(year=maxdate_split[0],month=maxdate_split[1],day = maxdate_split[2],
                               hour=maxdate_split[3],minute=maxdate_split[4],second=maxdate_split[5])

minyr = int((mindate_dt-datetime.timedelta(hours=12)).year)  #-1 day to make sure to cover for if timezone gmt +12
maxyr = int((maxdate_dt+datetime.timedelta(hours=12)).year)  


glaciers = pd.read_csv(f'{input_path}/Glacier_list.csv',header = 2)

glaciers = glaciers[glaciers['RGI index'] == Id]
height = xr.load_dataset(f'{input_path}/ERA5L/geopotential.nc')
g = 9.80665      # grav accell
M = 0.02896968   #molar mass of dry air (in kg)
R = 8.314462618  #gas constant J/(mol*K)
height['z'] = height.z/9.80665
lapse_rate = -6.5 #deg/km

def barometric_p(P,dh,T):
    a = lapse_rate*dh/(1000*T)    #convert to m
    b = -(g*M)/(R*lapse_rate/1000)
    return P*(1+a)**b


years = np.arange(minyr,maxyr+1,1)
vars = ['Ta','DTa','Lwin','Precip','Psfc','Swin','u10','v10']


points_df = pd.read_csv(f'{root}/Forcing_data/{Id}/coords_out_{Id}.csv', header = 1)
lats = points_df['lat']
lons = points_df['lon']
index = 1  #keeps track of number of loops  (maybe should use enumerate() for cleaner version)



num_points = len(points_df)
expected_paths = []
for i in range (num_points):
     expected_paths.append(f'{root}/Forcing_data/{Id}/{i}.parquet') #create list of expected filepaths
print('number of points EXTRAct')
print(num_points)

data_paths =[]
#if not all(os.path.isfile(p) for p in expected_paths):   #check if this step had already been completed
lapse_rates_map = xr.open_dataset(f'{input_path}/ERA5L/ERA5_Monthly_Hourly_Lapse_Rates_HMA_2000-2010.nc')
if True:
    for year in years:   # 
            for var in vars:
                data_paths.append(f'{input_path}/ERA5L/{year}/ERA5Land_HMA_{var}_{year}.nc')
    ds = xr.open_mfdataset( data_paths,chunks={'time': -1, 'latitude':1, 'longitude':1},parallel = False)  #the chunking is extremely important, 
                                                                                                            #reduce overhead by chunking to single points
                                                                                                            # so in next step, only one chunk needs to be loaded to extract the data

    
    for lon,lat in zip(lons,lats):   #iterate for each point
        point_pos = [lat,lon]  #lat,lon
        lapse_rates = lapse_rates_map.sel(latitude=point_pos[0], longitude=point_pos[1], method='nearest')
        df_lapse = lapse_rates['lapse rates mean'].to_dataframe()*1000

        deltaGMT = np.ceil(lon / 15.0)
        mindate_dt_GMT = mindate_dt - datetime.timedelta(hours = deltaGMT)  #get GMT times, to use for filtering ERA5L data
        maxdate_dt_GMT = maxdate_dt - datetime.timedelta(hours = deltaGMT)
        
        h_ERA5L = height.sel(latitude=point_pos[0], longitude=point_pos[1], method='nearest').z.values
        h_dem = points_df.elev_m.iloc[index-1]
        
        
        out_ds = ds.sel(latitude=point_pos[0], longitude=point_pos[1], method='nearest').sel(time=slice(mindate_dt_GMT,maxdate_dt_GMT)).compute()
        df_out = out_ds.to_dataframe().reset_index()
        df_out['hour'] = df_out.time.dt.hour
        df_out['month'] = df_out.time.dt.month
        df_out = pd.merge(df_out,df_lapse,on=['hour','month'])   #merge on the lapse rates based on month and hour!
        df_out['delta_T'] =  ((h_dem-h_ERA5L))*df_out['lapse rates mean']

           
        out_ds = ds.sel(latitude=point_pos[0], longitude=point_pos[1], method='nearest').sel(time=slice(mindate_dt_GMT,maxdate_dt_GMT)).compute()
        df_out = out_ds.to_dataframe().reset_index()

        df_out['sp'] = barometric_p(df_out['sp'],h_dem-h_ERA5L,df_out.t2m)
        df_out['Ws'] = (df_out['u10']**2 + df_out['v10']**2)**(0.5)  #wind speed
        
        df_out['t2m'] = df_out['t2m']+df_out['delta_T']   #apply T lapse rate  
        df_out['d2m'] = df_out['d2m']+df_out['delta_T']

        
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

        df_out = df_out.drop(['u10','v10','longitude_x','latitude_x','longitude_y','latitude_y','ssrd_mod',
                              'strd_mod','tp_mod','hour','month','lapse rates mean','delta_T'],axis=1)
        df_out.time = df_out.time + datetime.timedelta(hours=deltaGMT)
        df_out = df_out.iloc[1:]    # drop first row, has Nans for accummulated variables ()
        df_out.to_parquet(f'{root}/Forcing_data/{Id}/{index}.parquet',index=False)
        index+=1 
        #print(index)
    print(time.time()-t1)
    print('ERA5 extraction complete')
else:
    print('all_files_exist')