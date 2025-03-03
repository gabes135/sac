#!/bin/bash -l
#$ -l h_rt=11:59:99
#$ -j y
#$ -o /dev/null
#$ -N p_scan

p_0=$1
delta_p=$2

ac=$3
ar=$4

p=$(echo "scale=2; $p_0 + $delta_p * ($SGE_TASK_ID-1)" | bc)


echo "Start $JOB_NAME-$JOB_ID: $(date) $ac $ar $p"

module load julia
julia sac_edge.jl $ac $ar $p >> jobs/p_scan/task_$SGE_TASK_ID.txt 2>&1
wait
echo "End $JOB_NAME-$JOB_ID: $(date) $ac $ar $p"
