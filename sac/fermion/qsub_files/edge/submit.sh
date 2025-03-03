#!/bin/bash -l


mkdir -p out_files 

### Uncomment the single # commented lines corresponding to the
### type of run you want to perform
## For the scans, edit the parameters above 'qsub ...'
## Make sure there are no spaces when setting variable values!

### Plain Run ###
#qsub qsub.sh


### p Scan ###

# mkdir -p jobs/p_scan
# p_0=0.3 # initial value
# delta_p=0.1 # spacing in scan
# N_p=5 #number of values in scan

# ac=0.0
# ar=0.5

# qsub -t 1-$N_p qsub_p_scan.sh $p_0 $delta_p $ac $ar


## a_c Scan ###

# mkdir -p jobs/ac_scan
# ac_0=0.0 # initial value
# delta_ac=0.1 # spacing in scan
# N_p=7 #number of values in scan

# p=0.5
# ar=0.7

# qsub -t 1-$N_p qsub_ac_scan.sh $ac_0 $delta_ac $p $ar


### a_r Scan ###

mkdir -p jobs/ar_scan
ar_0=0.8 # initial value
delta_ar=0.1 # spacing in scan
N_p=2 #number of values in scan

p=0.5
ac=0.0

qsub -t 1-$N_p qsub_ar_scan.sh $ar_0 $delta_ar $p $ac
