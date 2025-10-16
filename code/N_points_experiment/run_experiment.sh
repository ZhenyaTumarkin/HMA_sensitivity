#!/bin/bash


n_points=($(seq 15 5 100))
# n_points=(5)
echo "${n_points[@]}"
jobs_pooled=50

root="/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code"

# rm /nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/N_points_experiment_RGI60-13.19847.csv #clear previous
#
tot_jobs=0
for n_point in "${n_points[@]}"; do
    num_iter=$((($n_point*15+$jobs_pooled+1)/($jobs_pooled)))   # ceiling function (round up)
    echo "$num_iter"
    outloc="/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/N_points_experiment/${n_point}"
    echo "$outloc"
    # jid1=$(sbatch --parsable --array=1-15%5 --export=num_points=$n_point,outlocation=$outloc ${root}/bash_scripts/launch_array_preprocess.sh)
    jid2=$(sbatch --parsable  --export=outlocation=$outloc ${root}/bash_scripts/run_list.sh)
    jid3=$(sbatch --parsable --dependency=afterok:$jid2 --export=outlocation=$outloc,jobs=$jobs_pooled --array=1-$num_iter:1 ${root}/bash_scripts/array_run_model.sh)
  
    #   #postprocess; store glacierwide mass balance to a file
   
    # # jid4=$(sbatch --parsable --dependency=afterany:$jid3 --export=n_point=$n_point ${root}/N_points_experiment/postprocess_mb.sh   )    

    # # prev_jid=$jid

    # echo "$n_point"
    tot_jobs=$(($tot_jobs+$num_iter))

done

echo "total $tot_jobs"
