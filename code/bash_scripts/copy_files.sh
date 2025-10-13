#!/bin/bash

###list of source destination to copy (for kyzylsu test)
remote_sensing="/fs3/group/pelligrp/Remote_Sensing_Data/Regional_Global"
data_path="/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data"


# cp -r ${remote_sensing}/ASTER_GDEM3/ASTGTMV003_N39E071 ${data_path}/preprocessing/ASTER_DEM/
# cp -r ${remote_sensing}/Rounce2021_131415/hd_tifs/13/13.19847_hdts_m.tif ${data_path}/preprocessing/Rounce2021_deb_thick/13/
# cp -r /fs3/group/pelligrp/REANALYSIS/ERA5LAND/2000 ${data_path}/preprocessing/ERA5L/
# cp  /fs3/group/pelligrp/PERSONAL_FOLDERS/Zhenya/HMA_test_glaciers/data/ERA5L/geopotential.nc ${data_path}/preprocessing/ERA5L/


cp -r /fs3/group/pelligrp/REANALYSIS/ERA5LAND/2001 ${data_path}/preprocessing/ERA5L/
cp -r /fs3/group/pelligrp/REANALYSIS/ERA5LAND/2002 ${data_path}/preprocessing/ERA5L/