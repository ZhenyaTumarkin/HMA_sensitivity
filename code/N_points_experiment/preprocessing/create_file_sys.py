import pandas as pd
import numpy as np
import os
import sys

root = '/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/N_points_experiment'

N_point_txt = open(f'{root}/list_n.txt').read().split(' ')
N_point_list = np.array(N_point_txt).astype(int)

glaciers = pd.read_csv(f'/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing/Glacier_list.csv',header = 2)
subfolders = ['ERA5L_at_gridpoint','ERA5L_radiation_partitioned','ERA5L_RH']   # make most of these on the fly/ write on top of itself after testing
glacier_ids = glaciers['RGI index']

if os.path.isdir(f'{root}/preprocessing/All_glaciers')==False:
    os.mkdir(f'{root}/preprocessing/All_glaciers')

for glacier in glacier_ids:
    for n_points in N_point_list:
        if os.path.isdir(f'{root}/preprocessing/All_glaciers/{glacier}/{n_points}')==False:
            os.mkdir(f'{root}/preprocessing/All_glaciers/{glacier}/{n_points}')
        
        for folder in subfolders:
            if os.path.isdir(f'{root}/preprocessing/All_glaciers/{glacier}/{n_points}/{folder}')==False:
                os.mkdir(f'{root}/preprocessing/All_glaciers/{glacier}/{n_points}/{folder}')

print('folders made')