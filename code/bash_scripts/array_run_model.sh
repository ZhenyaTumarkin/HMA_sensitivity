#!/bin/bash

##----------------------------------------------------------------
# running a multiple independent jobs
#----------------------------------------------------------------

##  Defining options for slurm how to run
#----------------------------------------------------------------
#SBATCH --job-name=Glacier_arr_preprocess_test
#SBATCH --output=/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/outlog/run_model_%A-%a.log 
# %A and %a are placeholders for the jobid and taskid, resp.
# Number of CPU cores to use within one node
#SBATCH -c 1

#Send emails when a job starts, it is finished or it exits
#SBATCH --mail-user=evgeny.tumarkin@ist.ac.at
#SBATCH --mail-type=SUBMIT,END,FAIL

# Define the number of hours the job should run. 
# Maximum runtime is limited to 10 days, ie. 240 hours
#SBATCH --time=00-00:15
#SBATCH --constraint=matlab

# Define the amount of RAM used by your job in GigaBytes
# In shared memory applications this is shared among multiple CPUs
#SBATCH --mem=4G

# Do not requeue the job in the case it fails.
#SBATCH --no-requeue

# Do not export the local environment to the compute nodes
#SBATCH --export=NONE
set -euo pipefail
unset SLURM_EXPORT_ENV
source /etc/profile

#read csv first

module load matlab/R2024b



csv_file="/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing/run_list.csv"

first_col=()
second_col=()
while IFS=, read -r col1 col2 col3 _; do
    first_col+=("$col2")
    second_col+=("$col3")
done < <(tail -n +2 "$csv_file")

echo "${first_col[@]}"
echo "${SLURM_ARRAY_TASK_ID}"

rgiid=${first_col[${SLURM_ARRAY_TASK_ID}-1]}
point_id=${second_col[${SLURM_ARRAY_TASK_ID}-1]}



if [ -n "${outlocation:-}" ]; then   #check if outlocation set previously
    echo "outlocation exists!"
    
    srun --cpu_bind=verbose matlab -nodesktop -nojvm -nosplash -r "clear all;\
    glacier_id='$rgiid'; point_id='$point_id'; outlocation='$outlocation';\
    run('/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/Run_model/Launcher_point.m');\
    exit;"
else
echo "no outlocation"
    srun --cpu_bind=verbose matlab -nodesktop -nojvm -nosplash -r "clear all;\
    glacier_id='$rgiid'; point_id='$point_id';\
    run('/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/Run_model/Launcher_point.m');\
    exit;"
fi
