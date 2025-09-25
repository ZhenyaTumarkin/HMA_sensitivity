import pandas as pd
import numpy as np 
import os


root_data = '/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing'
glacier_list = pd.read_csv(f'{root_data}/Glacier_list.csv',header = 2)


run_list = pd.DataFrame({'glacier_id':[],'point_id':[]})
for glacier in glacier_list['RGI index'].iloc[:2]:
    point_list =  pd.read_csv(f'{root_data}/All_glaciers/{glacier}/coords_out_{glacier}.csv',header = 1)
    num_points = len(point_list)
    points = np.arange(1,num_points +1,1)
    print(points)
    point_df = pd.DataFrame({'point_id':points})
    
    point_df['glacier_id'] = glacier
    
    run_list = pd.concat([point_df,run_list])
run_list.point_id = run_list.point_id.astype(int)
run_list_reord = pd.DataFrame({'glacier_id':run_list.glacier_id,'point_id':run_list.point_id.astype(int)})
run_list_reord.to_csv(f'{root_data}/run_list.csv')