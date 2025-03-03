# sac - Julia code to run the Stochastic Analytic Continuation Method
Currently only supports unconstrained sampling and the monotonic edge constrained parameterizations for fermionic spectral functions.

## Instructions for submit batch jobs using `qsub` for the monotonic edge constrained parameterization


### Edit the `submit.sh` file

Uncomment the lines correpsonding to the type of job you want to run:
* Plain run using the parameters you set in `in_edge.in`
* Scan over $p$
* Scan over $A_c$
* Scan over $A_r$

For the scans, set the values scanned over and the values of the two parameters fixed during the scan using the lines above `qsub -t ...`. For example, if running a scan over the weight $A_c$ from 0.0 to 0.7, set `ac_0=0`, `delta_ac=0.1`, and `N_p=8`. Also set `p` and `ar` to the desired values (`p=0.5` and `ar=0.5`, for example). *Note, do not include any spaces before and after equal sign when assigning variable values.*




### Subsmit the batch job
To submit the batch job, simply run `sh submit.sh`. The scan will be run as a job array and the qsub output files will be written to a folder `jobs/[parameter]_scan/task_[#].txt`, where `parameter` is the what is scanned over and `#` is the job array number corresponding to each value in the scan.

