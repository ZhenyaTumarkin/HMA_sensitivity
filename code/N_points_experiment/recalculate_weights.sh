#!/bin/bash

module load miniforge3
conda activate HMA

code_root="/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/preprocessing"

csv_file="/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing/Glacier_list.csv"

# Skip header if needed: tail -n +2
first_col=()
while IFS= read -r value; do
    first_col+=("$value")
done < <(tail -n +4 "$csv_file" | cut -d',' -f2)    #

echo "${first_col[@]}"

for glacier in {0..0}
do
rgiid=${first_col[$glacier]}

for num_points in $(seq 5 5 100)
do
outlocation="/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/N_points_experiment/${num_points}"


python ${code_root}/extract_points_kmeans.py --rgiid $rgiid --N_points $num_points --outlocation $outlocation
done

done


