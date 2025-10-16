import pandas as pd
import numpy as np
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--outlocation',type=str,required = False,help = 'outlocation')
args = parser.parse_args()
outlocation=args.outlocation

if outlocation:
    root = f'{outlocation}'
else:
    root = '/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing'


print(root)
glaciers = pd.read_csv(f'/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing/Glacier_list.csv',header = 2)
#subfolders = ['ERA5L_at_gridpoint','ERA5L_radiation_partitioned','ERA5L_RH']   # make most of these on the fly/ write on top of itself after testing
glacier_ids = glaciers['RGI index']


for glacier in glacier_ids:
    os.makedirs(f'{root}/Forcing_data/{glacier}',exist_ok=True)

    
    # for folder in subfolders:
    #     if os.path.isdir(f'{root}/Forcing_data/{glacier}/{folder}')==False:
    #         os.mkdir(f'{root}/Forcing_data/{glacier}/{folder}')

print('folders made')

