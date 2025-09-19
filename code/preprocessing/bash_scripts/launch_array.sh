#!/bin/bash
##----------------------------------------------------------------
# running a multiple independent jobs
#----------------------------------------------------------------

##  Defining options for slurm how to run
#----------------------------------------------------------------
#SBATCH --job-name=Glacier_arr_preprocess_test
#SBATCH --output=/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/outlog/Glaconetes_%A-%a.log 
# %A and %a are placeholders for the jobid and taskid, resp.
# Number of CPU cores to use within one node
#SBATCH -c 1

# Define the number of hours the job should run. 
# Maximum runtime is limited to 10 days, ie. 240 hours
#SBATCH --time=00-00:10
#SBATCH --constraint=matlab

# Define the amount of RAM used by your job in GigaBytes
# In shared memory applications this is shared among multiple CPUs
#SBATCH --mem=4G

# Do not requeue the job in the case it fails.
#SBATCH --no-requeue

# Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

# load the respective software module(s) you intend to use
#----------------------------------------------------------------
# none needed for this example

# define sequence of jobs to run as you would do in a BASH script
# use variable $SLURM_ARRAY_TASK_ID to address individual behaviour
# in different iteration of the script execution
#----------------------------------------------------------------

csv_file="/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing/Glacier_list.csv"

# Skip header if needed: tail -n +2
first_col=()
while IFS= read -r value; do
    first_col+=("$value")
done < <(tail -n +4 "$csv_file" | cut -d',' -f2)

rgiid=${first_col[${SLURM_ARRAY_TASK_ID}]}

srun --cpu_bind=verbose  bash /nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/preprocessing/bash_scripts/preprocessing.sh $rgiid
