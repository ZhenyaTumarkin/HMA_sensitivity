import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import argparse
import os

print('python start')
data_root = '/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data'

glacier_ids = pd.read_csv(f'{data_root}/preprocessing/Glacier_list.csv',header = 2)
glacier_ids = glacier_ids['RGI index']

for Id in glacier_ids:
    region_num = Id.split('-')[1][0:2]   # extract RGI region number from Id
    Id_num = Id.split('-')[1]
    region_names = {
        '13':'CentralAsia',
        '14':'SouthAsiaWest',
        '15':'SouthAsiaEast'
    }

    run_name = 'run_1.000'

    ########now load model output
    #first check if all glacier points ran
    points_all = pd.read_csv(f'{data_root}/preprocessing/All_glaciers/{Id}/coords_out_{Id}.csv',header=1)
    num_points = len(points_all)

    paths_out = np.array([f'{data_root}/Outputs/{Id}/{run_name}/{point_i+1}_results.parquet' for point_i in range(num_points)])
    print(paths_out[0])
    #calculate mass change per m^2 (using weights)
    #then multiply by area total

    exist_files = np.array([os.path.isfile(p) for p in paths_out])
    if all(exist_files) == False:
        print(f"points_not computed {np.arange(1,num_points+1,1)[~exist_files].astype(int)}")
        
        

    else:
        print('all files exist')

    #drop not computed, then recompute weights
    points_computed_list = np.arange(1,num_points+1,1)[exist_files].astype(int)
    points_computed_df = points_all[exist_files]
    points_computed_df.weight = points_computed_df.weight/np.nansum(points_computed_df.weight)


    print(points_computed_list)
    if not any(points_computed_list):
        print('any')
        print('No outputs')
    else:
        for i,point_id in enumerate(points_computed_list):
            
            out_file = pd.read_parquet(f'{data_root}/Outputs/{Id}/{run_name}/{point_id}_results.parquet')
            #remove spin up
            date_start = out_file.Date.iloc[0]
            date_spin_complete = date_start + datetime.timedelta(days = 30)
            date_end =date_spin_complete + datetime.timedelta(days = 365*4+1)
            out_file = out_file[(out_file.Date>date_spin_complete)&(out_file.Date<=date_end)]


            if i==0:
                ice_we = (out_file.ICE-out_file.ICE.iloc[0])*points_computed_df.weight.iloc[i]
                snow_we = (out_file.SWE-out_file.SWE.iloc[0])*points_computed_df.weight.iloc[i]
                
            else:
                ice_we += (out_file.ICE-out_file.ICE.iloc[0])*points_computed_df.weight.iloc[i]
                snow_we += (out_file.SWE-out_file.SWE.iloc[0])*points_computed_df.weight.iloc[i]


        total_mb = ice_we+snow_we   #(kg per m^2)

        hugonnet_path = f'{data_root}/postprocessing/Hugonnet_2021/dh_{region_num}_rgi60_pergla_rates.parquet'
        glacier_data = pd.read_parquet(hugonnet_path,filters = [('rgiid','==',Id),('period','==','2000-01-01_2010-01-01')],engine="pyarrow")
        model_mb = (total_mb.iloc[-1]*1/4) 
        geodetic_mb = (glacier_data.dmdt/glacier_data.area).iloc[0] * 1e12
        print(run_name)
        print('model: ' + str(model_mb) + ' mm we yr^-1')
        print('geodetic: ' + str(geodetic_mb) + ' mm we yr^-1')
        


        #### now read in current previous run files
        if os.path.isfile(f'{data_root}/Outputs/{Id}/tp_calib.csv'):
            in_df = pd.read_csv(f'{data_root}/Outputs/{Id}/tp_calib.csv')
            if len(in_df)>=1:   # maybe add check for == for safety? (or try except)
                #### secant (adapted) search method, Now not purel;y iterative, but chooses two closest points
                try:
                    p_lower_bound = np.nanmax(in_df.p_fact[(in_df.model_clim_mb < geodetic_mb)])
                    p_upper_bound = np.nanmin(in_df.p_fact[in_df.model_clim_mb > geodetic_mb])
                    print(p_lower_bound)
                    print(p_upper_bound)
                    mb_lower = in_df.model_clim_mb[in_df.p_fact==p_lower_bound].values
                    mb_upper = in_df.model_clim_mb[in_df.p_fact==p_upper_bound].values
                    print(mb_lower,mb_upper)
                    flag_bounds = True
                except ValueError:    # if both values are still below (or above) use normal secant method
                    p_lower_bound =  out_df.p_fact.iloc[-1]
                    p_upper_bound =  out_df.p_fact.iloc[-2]
                    mb_lower      =  out_df.model_clim_mb.iloc[-1]
                    mb_upper      =  out_df.model_clim_mb.iloc[-2]
                    flag_bounds = False
                new_row = pd.DataFrame({'p_fact':np.array([in_df.next_p.iloc[-1]]),'model_clim_mb':model_mb,'next_p':np.array([0]).astype(float)})
                out_df = pd.concat([in_df, new_row], ignore_index=True)

                p_new = p_lower_bound +(geodetic_mb- mb_lower)* (p_lower_bound - p_upper_bound
                                                                            )/(mb_lower -mb_upper)
                
                #check if p_new in between p_lower and p_upper
                if flag_bounds:
                    if (p_new < p_lower_bound) | (p_new >p_upper_bound):
                        ######bisection method!
                        p_new = (p_lower_bound + p_upper_bound)/2

                print(p_new)
                out_df.loc[out_df.index[-1],'next_p'] = p_new
        else:    
            if model_mb < geodetic_mb:
                p_new = 2
                
            elif model_mb > geodetic_mb:
                p_new = 0.5
            out_df = pd.DataFrame({'p_fact':np.array([1]),'model_clim_mb':model_mb,'next_p':np.array([p_new]).astype(float)})
        print(out_df)
        # out_df.to_csv(f'{data_root}/Outputs/{Id}/tp_calib.csv',index=False)

