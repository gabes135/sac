#!/bin/bash -l
#$ -l h_rt=11:59:99
#$ -j y
#$ -o /dev/null
#$ -N ac_scan



ac_0=$1
delta_ac=$2

p=$3
ar=$4

ac=$(echo "scale=2; $ac_0 + $delta_ac * ($SGE_TASK_ID-1)" | bc)


echo "Start $JOB_NAME-$JOB_ID: $(date) $ac $ar $p"

module load julia
julia sac_edge.jl $ac $ar $p >> jobs/ac_scan/task_$SGE_TASK_ID.txt 2>&1
wait
echo "End $JOB_NAME-$JOB_ID: $(date) $ac $ar $p"
