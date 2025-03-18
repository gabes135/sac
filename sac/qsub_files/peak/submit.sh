#!/bin/bash -l


mkdir -p out_files 

### Uncomment the single # commented lines corresponding to the
### type of run you want to perform
## For the scans, edit the parameters above 'qsub ...'
## Make sure there are no spaces when setting variable values!

### Plain Run ###
#qsub qsub.sh


## A0 Scan ###

mkdir -p jobs/A0_scan
A0_0=0.1 # initial value
delta_A0=0.1 # spacing in scan
N_p=9 #number of values in scan

Np=1

qsub -t 1-$N_p qsub_a0_scan.sh $A0_0 $delta_A0 $Np
