import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument('--run_name',type=str,required = True,help = 'input RGIId')

args = parser.parse_args()
run_name = args.run_name



print('python start')
data_root = '/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data'

glacier_ids = pd.read_csv(f'{data_root}/preprocessing/Glacier_list.csv',header = 2)
glacier_ids = glacier_ids['RGI index'].iloc[:1]

print(glacier_ids)
for Id in glacier_ids:
    region_num = Id.split('-')[1][0:2]   # extract RGI region number from Id
    Id_num = Id.split('-')[1]
    region_names = {
        '13':'CentralAsia',
        '14':'SouthAsiaWest',
        '15':'SouthAsiaEast'
    }

    p_mod_path = f'/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/Outputs/{Id}/tp_calib.csv'
    

    ########now load model output
    #first check if all glacier points ran
    points_all = pd.read_csv(f'{data_root}/N_points_experiment/{run_name}/coords_out_{Id}.csv',header=1)
    num_points = len(points_all)

    paths_out = np.array([f'{data_root}/N_points_experiment/{run_name}/{point_i+1}_results.parquet' for point_i in range(num_points)])
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


    print(points_computed_df)
    if not any(points_computed_list):
        print('any')
        print('No outputs')



    
    else:
        print(points_computed_df.weight.values)
        for i,point_id in enumerate(points_computed_list):
            #print(i,point_id)
            out_file = pd.read_parquet(f'{data_root}/N_points_experiment/{run_name}/{point_id}_results.parquet')
            #remove spin up
            date_start = out_file.Date.iloc[0]
            date_spin_complete = date_start + datetime.timedelta(days = 365)
            date_end =date_spin_complete + datetime.timedelta(days = 365)
            out_file = out_file[(out_file.Date>date_spin_complete)&(out_file.Date<=date_end)]


            if i==0:
                ice_we = (out_file.ICE-out_file.ICE.iloc[0])*points_computed_df.weight.iloc[i]
                snow_we = (out_file.SWE-out_file.SWE.iloc[0])*points_computed_df.weight.iloc[i]
                
            else:
                ice_we += (out_file.ICE-out_file.ICE.iloc[0])*points_computed_df.weight.iloc[i]
                snow_we += (out_file.SWE-out_file.SWE.iloc[0])*points_computed_df.weight.iloc[i]


        total_mb = ice_we+snow_we   #(kg per m^2)

        ####now load current N_points results
        print((total_mb))
        if os.path.isfile(f'{data_root}/N_points_experiment_{Id}.csv'):
            df_exp = pd.read_csv(f'{data_root}/N_points_experiment_{Id}.csv',index_col=False)
            num_cols = len(df_exp.columns)
          
            df_exp[f'run_{num_cols+1}'] = total_mb.values
        else:
            df_exp = pd.DataFrame()
            df_exp['run_1'] = total_mb.values
        print('saving_fikle')
        print(total_mb.iloc[-1])
        df_exp.to_csv(f'{data_root}/N_points_experiment_{Id}.csv',index=False)