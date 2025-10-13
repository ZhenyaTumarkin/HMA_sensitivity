#!/bin/bash

####this is for debugging and finding issues with single point runs (that failed previously)

n_point=95
point_id=25



root="/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code"
outloc="/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/N_points_experiment/${n_point}"
jid1=$(sbatch --parsable --array=1 --export=num_points=$n_point,outlocation=$outloc ${root}/bash_scripts/launch_array_preprocess.sh)
rm /nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing/run_list.csv
echo  ",glacier_id,point_id"> "/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing/run_list.csv"
echo  ",RGI60-13.19847,$point_id">> "/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing/run_list.csv"

sbatch --parsable --dependency=afterok:$jid1 --array=1 --export=outlocation=$outloc ${root}/bash_scripts/array_run_model.sh