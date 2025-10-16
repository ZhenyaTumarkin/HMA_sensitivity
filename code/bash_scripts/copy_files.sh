#!/bin/bash

###list of source destination to copy (for kyzylsu test)
remote_sensing="/fs3/group/pelligrp/Remote_Sensing_Data/Regional_Global"
data_path="/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data"



# cp  /fs3/group/pelligrp/PERSONAL_FOLDERS/Zhenya/HMA_test_glaciers/data/ERA5L/geopotential.nc ${data_path}/preprocessing/ERA5L/
# cp /fs3/group/pelligrp/REANALYSIS/ERA5_Lvls/Monthly_Hourly_LapseRates/ERA5_Monthly_Hourly_Lapse_Rates_HMA_2000-2010.nc ${data_path}/preprocessing/ERA5L/
# cp -r /fs3/group/pelligrp/REANALYSIS/ERA5LAND/2000 ${data_path}/preprocessing/ERA5L/
# cp -r /fs3/group/pelligrp/REANALYSIS/ERA5LAND/2001 ${data_path}/preprocessing/ERA5L/
# cp -r /fs3/group/pelligrp/REANALYSIS/ERA5LAND/2002 ${data_path}/preprocessing/ERA5L/

#KYZYLSU
# cp -r ${remote_sensing}/ASTER_GDEM3/ASTGTMV003_N39E071 ${data_path}/preprocessing/ASTER_DEM/
# cp -r ${remote_sensing}/Rounce2021_131415/hd_tifs/13/13.19847_hdts_m.tif ${data_path}/preprocessing/Rounce2021_deb_thick/13/

#ALL HMA
# cp -r ${remote_sensing}/Rounce2021_131415/hd_tifs/13 ${data_path}/preprocessing/Rounce2021_deb_thick/
# cp -r ${remote_sensing}/Rounce2021_131415/hd_tifs/14 ${data_path}/preprocessing/Rounce2021_deb_thick/
cp -r ${remote_sensing}/Rounce2021_131415/hd_tifs/15 ${data_path}/preprocessing/Rounce2021_deb_thick/


# echo "Copying ASTER tiles: 20 < N < 45, 65 < E < 110"
# copied=0
# for lat in $(seq 21 44); do
#     echo "$lat"
#     for lon in $(seq 66 109); do
#         lat_str=$(printf "N%02d" $lat)
#         lon_str=$(printf "E%03d" $lon)
        
#         # Use find to locate all matching files/directories
#         find ${remote_sensing}/ASTER_GDEM3/ -maxdepth 1 -name "ASTGTMV003_${lat_str}${lon_str}*" | while read item; do
#             if [ -e "$item" ]; then
#                 cp -rv "$item" ${data_path}/preprocessing/ASTER_DEM/
#                 #echo "Copied: $(basename "$item")"
#                 ((copied++))
#             fi
#         done
#     done
#     echo "copied total $copied"
# done

# echo "Total copied: $copied files/directories"