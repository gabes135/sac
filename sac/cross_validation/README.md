## Instructions for running the cross validation procedure for model selection

The functions contained in this folder can be used to run a cross validation procedure as a tool for model selection. This can be read about in depth in our [PRE article](https://arxiv.org/pdf/2406.06763).

## Instructions for running the cross validation procedure

* Collect your QMC-generated bins in a `cor.dat` and `tgrid.dat` file as detailed [here](../../README.md) and copy to the `process_G` folder in this directory.
* Generate the sampling and validation data files from `cor.dat` using the `generate_sets.sh` script. This script takes two inputs:
	1. `folder`: the name of the folder containing the input/output files
	2. `reps`: the number of sampling and validation data sets to generate
For each rep, the bins in `cor.dat` are split into two mutually exclusive sets and a `t.in` file is geneated for each (written to the directory `in_files/{folder}` and labeled by the rep number and a/b).
* Run `cross_val.jl` for each rep and for each SAC parameterization you want to test. This program takes two inputs:
	1. `rep`: the specific rep number
	2. `param`: the SAC parameterization (`free`, `peak`, or `edge`)
Before running the program, modify the `in_{param}.in` file with the appropriate input/output folder and your desired settings (`A0` for `peak`, for example). This program runs a SAC anneal for each of the two in_files (a/b) and ouputs the SAC $G(\tau)$ at each step in the anneal to file (`GSAC.csv`).
* Process the output SAC $G(\tau)$ results to calculate the cross validation $\chi^2$ value for each parameterization, averaged over all reps. This is done using the `calc_cv.jl` program, which takes three inputs:
	1. `folder`: the name of the folder containing the input/output files
	2. `reps`: the number of sampling and validation data sets
	3. `param`: the SAC parameterization (`free`, `peak`, or `edge`). If `peak` or `edge`, include the parameters that are used in the ouput directory, i.e. `peak/Np_01/A0_0.700`
For each rep, the validation $\chi^2$ is calculated treating set a as the sampling set and b as the validation, then the roles are reversed. The final values for this rep are the average of these two 'rotations.' The results are written to the file `out_files/{folder}/chi2/{param}.csv`, and contains the validation $\chi^2$ at each temperature step in the annealing run (columns) for each rep (rows). The results can be analyzed in the included [Jupyter notebook](./plot_results.ipynb).






