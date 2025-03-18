#!/bin/bash -l
#$ -l h_rt=11:59:99
#$ -j y
#$ -o jobs/out.txt
#$ -N edge

mkdir -p jobs # creates folder for qsub out files
mkdir -p out_files # creates folder for SAC out files

echo "Start $JOB_NAME-$JOB_ID: $(date)"

module load julia

julia sac_edge.jl

echo "End $JOB_NAME-$JOB_ID: $(date) "