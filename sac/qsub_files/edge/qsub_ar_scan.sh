#!/bin/bash -l
#$ -l h_rt=11:59:99
#$ -j y
#$ -o /dev/null
#$ -N ar_scan

ar_0=$1
delta_ar=$2

p=$3
ac=$4

ar=$(echo "scale=2; $ar_0 + $delta_ar * ($SGE_TASK_ID-1)" | bc)


echo "Start $JOB_NAME-$JOB_ID: $(date) $ac $ar $p"

module load julia
julia sac_edge.jl $ac $ar $p >> jobs/ar_scan/task_$SGE_TASK_ID.txt 2>&1
wait
echo "End $JOB_NAME-$JOB_ID: $(date) $ac $ar $p"
