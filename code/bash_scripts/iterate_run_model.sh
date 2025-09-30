#!/bin/bash

bash_path="/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/bash_scripts"
prev_jid=""
num_iterations=5
for iter in $(seq 1 $num_iterations); do
    if [ -z "$prev_jid" ]; then
        jid1=$(sbatch --parsable --array=1-30:1 ${bash_path}/array_run_model.sh)
    else
        jid1=$(sbatch --parsable --dependency=afterok:$prev_jid --array=1-30:1 ${bash_path}/array_run_model.sh)
    fi

    jid2=$(sbatch --parsable --dependency=afterok:$jid1 ${bash_path}/postprocessing.sh $iter)
    prev_jid=$jid2
    echo "$iter"
done
