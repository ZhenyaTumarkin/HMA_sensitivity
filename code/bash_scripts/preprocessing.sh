#!/bin/bash
#  source /nfs/scistore18/pelligrp/etumarki/.conda/envs

module load miniforge3
module load matlab

datemin="2000:09:01:00:00:00"   #"$2"   #"2021:02:28:05:00:00"   YYYY:MM:DD:HH:MM:SS   TIME UTC
datemax="2004:12:25:23:00:00"    #"$3"

rgiid="$1"

conda activate HMA




python /nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/preprocessing/create_file_sys.py
python /nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/preprocessing/extract_points_kmeans.py --rgiid $rgiid
python /nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/preprocessing/extract_ERA5L.py --rgiid $rgiid --datemin $datemin --datemax $datemax

matlab -nodesktop -nojvm -nosplash -r "clear all;Id = '$rgiid'; \
run('/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/preprocessing/temp_tdew_press/relative_humidity_calc.m');\
clear all;glacier_id = '$rgiid';\
run('/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/preprocessing/Radiation_split/Radiation_Partition_ERA5_point.m')\
;exit;"

exit

