import pandas as pd
import numpy as np
import os
import sys

root = '/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing'
glaciers = pd.read_csv(f'{root}/Glacier_list.csv',header = 2)
subfolders = ['ERA5L_at_gridpoint','ERA5L_radiation_partitioned','ERA5L_RH']   # make most of these on the fly/ write on top of itself after testing
glacier_ids = glaciers['RGI index']

if os.path.isdir(f'{root}/All_glaciers')==False:
    os.mkdir(f'{root}/All_glaciers')

for glacier in glacier_ids:
    if os.path.isdir(f'{root}/All_glaciers/{glacier}')==False:
        os.mkdir(f'{root}/All_glaciers/{glacier}')
    
    for folder in subfolders:
        if os.path.isdir(f'{root}/All_glaciers/{glacier}/{folder}')==False:
            os.mkdir(f'{root}/All_glaciers/{glacier}/{folder}')

print('folders made')

