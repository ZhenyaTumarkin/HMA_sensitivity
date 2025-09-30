#!/bin/bash
#  source /nfs/scistore18/pelligrp/etumarki/.conda/envs

#SBATCH --job-name=Postprocess
#SBATCH --output=/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/outlog/post-%j.log   
#            %j is a placeholder for the jobid
#

#Define the number of hours the job should run. 
#Maximum runtime is limited to 10 days, ie. 240 hours
#SBATCH --time=00-00:10
#
#Define the amount of RAM used by your job in GigaBytes
#SBATCH --mem=4G
#
#Send emails when a job starts, it is finished or it exits
#SBATCH --mail-user=evgeny.tumarkin@ist.ac.at
#SBATCH --mail-type=SUBMIT,END,FAIL
#
#Pick whether you prefer requeue or not. If you use the --requeue
#option, the requeued job script will start from the beginning, 
#potentially overwriting your previous progress, so be careful.
#For some people the --requeue option might be desired if their
#application will continue from the last state.
#Do not requeue the job in the case it fails.
#SBATCH --no-requeue
#
#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
#
#for single-CPU jobs make sure that they use a single thread
export OMP_NUM_THREADS=1
#
#load the respective software module you intend to use
module load miniforge3
conda activate HMA
#
#
#run the respective binary through SLURM's srun
srun --cpu_bind=verbose  python /nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/postprocessing/compare_MB.py