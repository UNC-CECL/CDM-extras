## Instructions
#### Generation of a single parameter file and associated `DATA` directory
This script will generate a `params.par` file and an empty `DATA` folder in the same directory as `CDM_v2_0_parameter_file_initiation.py` is located.
1. Place `CDM_v2_0_parameter_file_initiation.py` in the directory in which you want to generate the parameter file.
2. Modify parameters in `CDM_v2_0_parameter_file_initiation.py` as necessary
3. Execute `CDM_v2_0_parameter_file_initiation.py`. 
4. To start the CDM, run `DUNE params.par`

#### Generation of a set of parameter files and `DATA` directories for running multiple simulations
This script will generate a series of directories containing params.par files andsempty DATA folders. For example, this `CDM_v2_0_batch_parameter_file_initiation.py` is currently set to create five CDM simulations, each with a unique `params.par` file. When run, this script will create five directories (`model_iter1`, `model_iter2`, ... , `model_iter5`). Each directory will contain the CDM program (`DUNE`), a unique `params.par` file, and an empty `DATA` folder. In addition, the parent directory will contain an `iter_metadata` file containing information about each simulation.
1. Place `CDM_v2_0_batch_parameter_file_initiation.py` in the parent directory in which you want to generate a series of CDM simulations.
2. Place `CDM_v2_0` in this same parent directory.
3. Modify parameters in `CDM_v2_0_batch_parameter_file_initiation.py` as necessary
4. Execute `CDM_v2_0_batch_parameter_file_initiation.py`. 
5. To start multiple CDM instances simulateously:
   1. Place `batch_CDM_run.sh` in parent directory
   2. Open file and specify number of cores to use (currently set to 3 cores). Save.
   3. Run `batch_CDM_run.sh`

