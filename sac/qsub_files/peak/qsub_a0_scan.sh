#!/bin/bash -l
#$ -l h_rt=11:59:99
#$ -j y
#$ -o /dev/null
#$ -N A0_scan

A0_0=$1
delta_A0=$2

Np=$3

A0=$(echo "scale=2; $A0_0 + $delta_A0 * ($SGE_TASK_ID-1)" | bc)


echo "Start $JOB_NAME-$JOB_ID: $(date) $A_0 $Np"

module load julia
julia sac_peak.jl $A0 $Np >> jobs/A0_scan/task_$SGE_TASK_ID.txt 2>&1
wait
echo "End $JOB_NAME-$JOB_ID: $(date) $A_0 $Np"
