#!/bin/bash


n_points=($(seq 5 5 100))

echo "${n_points[@]}"
prev_jid=""

root="/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code"

rm /nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/N_points_experiment_RGI60-13.19847.csv #clear previous
#
for n_point in "${n_points[@]}"; do
    outloc="/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/N_points_experiment/${n_point}"
    echo "$outloc"
    if [ -z "$prev_jid" ]; then
        jid1=$(sbatch --parsable --array=1 --export=num_points=$n_point,outlocation=$outloc ${root}/bash_scripts/launch_array_preprocess.sh)
    else
        jid1=$(sbatch --parsable --array=1 --dependency=afterany:$prev_jid --export=num_points=$n_point,outlocation=$outloc ${root}/bash_scripts/launch_array_preprocess.sh)
    fi
    jid2=$(sbatch --parsable --dependency=afterany:$jid1 ${root}/bash_scripts/run_list.sh)
    jid3=$(sbatch --parsable --dependency=afterok:$jid2 --export=outlocation=$outloc --array=1-$n_point:1 ${root}/bash_scripts/array_run_model.sh)
  
      #postprocess; store glacierwide mass balance to a file
   
    jid4=$(sbatch --parsable --dependency=afterany:$jid3 --export=n_point=$n_point ${root}/N_points_experiment/postprocess_mb.sh   )    

    prev_jid=$jid4

    echo "$n_point"

done


