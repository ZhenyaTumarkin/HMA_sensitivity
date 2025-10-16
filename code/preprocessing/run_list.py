import pandas as pd
import numpy as np 
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--glaciers',type=int,required = True,help = 'input number of glaciers')
parser.add_argument('--outlocation',type=str,required = False,help = 'outlocation')
args = parser.parse_args()
glaciers = args.glaciers
outlocation=args.outlocation


input_data = '/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing'
if outlocation:
    root = f'{outlocation}'
else:
    root = input_data

glacier_list = pd.read_csv(f'{input_data}/Glacier_list.csv',header = 2)


run_list = pd.DataFrame({'glacier_id':[],'point_id':[]})
if glaciers == 1:
    glacier_list =  [glacier_list['RGI index'].iloc[0]]
else:
    glacier_list =  glacier_list['RGI index'].iloc[:glaciers-1]


for glacier in glacier_list:
    point_list =  pd.read_csv(f'{root}/Forcing_data/{glacier}/coords_out_{glacier}.csv',header = 1)
    num_points = len(point_list)
    points = np.arange(1,num_points +1,1)
    
    point_df = pd.DataFrame({'point_id':points})
    
    point_df['glacier_id'] = glacier
    
    run_list = pd.concat([point_df,run_list])
run_list.point_id = run_list.point_id.astype(int)
run_list_reord = pd.DataFrame({'glacier_id':run_list.glacier_id,'point_id':run_list.point_id.astype(int)})
run_list_reord.to_csv(f'{root}/run_list.csv')