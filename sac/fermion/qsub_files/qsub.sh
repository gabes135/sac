#!/bin/bash -l
#$ -l h_rt=11:59:99
#$ -o jobs/out.txt
#$ -e jobs/err.txt



echo "Start $JOB_NAME-$JOB_ID: $(date)"

module load julia

julia run.jl $@
wait
echo "End $JOB_NAME-$JOB_ID: $(date) "