#!/bin/bash
#  source /nfs/scistore18/pelligrp/etumarki/.conda/envs

module load miniforge3
module load matlab

datemin="2000:09:01:00:00:00"   #"$2"   #"2021:02:28:05:00:00"   YYYY:MM:DD:HH:MM:SS   TIME UTC
datemax="2002:10:01:23:00:00"    #"$3"

rgiid="$1"
num_points="$2"
outlocation="$3"


    

echo "$num_points"
conda activate HMA




code_root="/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/preprocessing"


if [ -n "$outlocation" ]; then
    echo "outlocation provided : $outlocation"
    python ${code_root}/create_file_sys.py --outlocation $outlocation
    python ${code_root}/extract_points_kmeans.py --rgiid $rgiid --N_points $num_points --outlocation $outlocation
    python ${code_root}/extract_ERA5L.py --rgiid $rgiid --datemin $datemin --datemax $datemax --outlocation $outlocation
    matlab -nodesktop -nojvm -nosplash -r "clear all;Id = '$rgiid'; outlocation='$outlocation'; \
    run('/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/preprocessing/temp_tdew_press/relative_humidity_calc.m');\
    clear all;glacier_id = '$rgiid';outlocation='$outlocation';\
    run('/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/preprocessing/Radiation_split/Radiation_Partition_ERA5_point.m')\
    ;exit;"
else
    echo "no outlocation"
    python ${code_root}/create_file_sys.py
    python ${code_root}/extract_points_kmeans.py --rgiid $rgiid --N_points $num_points 
    python ${code_root}/extract_ERA5L.py --rgiid $rgiid --datemin $datemin --datemax $datemax
    matlab -nodesktop -nojvm -nosplash -r "clear all;Id = '$rgiid'; \
    run('/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/preprocessing/temp_tdew_press/relative_humidity_calc.m');\
    clear all;glacier_id = '$rgiid';\
    run('/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/preprocessing/Radiation_split/Radiation_Partition_ERA5_point.m')\
    ;exit;"
    
fi





exit

