#!/bin/bash -l
#$ -l h_rt=11:59:99
#$ -o logs/$JOB_NAME.txt
#$ -j y
#$ -t 1-1

a0=$1
delta_a=$2

a=$(echo "scale=2; $a0 + $delta_a * ($SGE_TASK_ID-1)" | bc)


echo "Start $JOB_NAME-$JOB_ID: $(date) $@"

module load julia

shift
# julia run.jl edge_p $p
julia run.jl $mode $a $@
wait
echo "End $JOB_NAME-$JOB_ID: $(date) $@"
