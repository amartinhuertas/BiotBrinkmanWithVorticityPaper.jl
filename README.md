

[![DOI](https://zenodo.org/badge/658182544.svg)](https://zenodo.org/badge/latestdoi/658182544)



# BiotBrinkmanWithVorticityPaper

Repository that holds the Julia software and computational results reported in Section 5.3 of the paper _"Robust finite element methods and solvers for the Biot-Brinkman equations in vorticity form"_ by _R. Caraballo, C. W. In, A. F. Martín, and R. Ruiz-Baier_. The results in Section 5.1 of the paper were not actually obtained with the software in this repository. However, the repository also contains Julia scripts to reproduce the experiments in this section in the particular case of mixed boundary conditions (not actually reported in the paper for the sake of brevity).

In particular:

* The accuracy test in Section 5.1 for mixed boundary conditions can be generated with the script available [here](https://github.com/amartinhuertas/BiotBrinkmanWithVorticityPaper.jl/blob/main/test/ConvergenceTests.jl).
* The plots in Section 5.3 can be generated using the [Pluto.jl](https://github.com/fonsp/Pluto.jl) interactive/reactive notebook available [here](https://github.com/amartinhuertas/BiotBrinkmanWithVorticityPaper.jl/blob/main/notebooks/BiotBrinkmanVorticityPreconditionerTests.jl). This notebook loads a set of BSON data files available [here](https://github.com/amartinhuertas/BiotBrinkmanWithVorticityPaper.jl/tree/main/data/BiotBrinkmanVorticityPreconditionerTests/915f21471c248e0371a237cf7cd0d833904de63b). This data files were in turn automatically generated using the [DrWatson.jl](https://github.com/JuliaDynamics/DrWatson.jl) script available [here](https://github.com/amartinhuertas/BiotBrinkmanWithVorticityPaper.jl/blob/main/experiments/BiotBrinkmanVorticityPreconditionerTests/run_experiment.jl). These scripts are designed such that the user might seamlessly run new combination of physical and discretization parameter values. See below for instructions. 

## Instructions to generate and visualize Riesz mapping preconditioner evaluation results (Section 5.3)

**IMPORTANT NOTE**: _As a pre-requisite to follow the instructions below, `DrWatson.jl` must be installed in the main julia 
environment, e.g., `v1.9` if you have Julia 1.9. Please install it in the 
main julia environment before following the instructions in the sequel._

1. Select those combinations of parameter-values to be run by editing the dictionary created by the  `generate_param_dicts()` function in the script available [here]().
2. Run the script:  
   ```bash
   cd experiments/BiotBrinkmanVorticityPreconditionerTests/
   julia run_experiment.jl
   ```
   Note that the `--project=XXX` flag is not required in the call to the `julia` command. 
   The script is smart enough  in order to locate the `Project.toml` of the project in an ancestor directory. Upon completion, the results are generated in the `data` directory 
of the project, in particular in the `data/BiotBrinkmanVorticityPreconditionerTests/commit-ID` folder, where `commit-ID` denotes the latest commit in the repository. If you have changes in your local repository, then the BSON data files will be generated at `data/BiotBrinkmanVorticityPreconditionerTests/commit-ID-dirty`. The script is intelligent enough so that if it is interrupted and re-run, it will only run those combinations which are pending.
3. Visualize the results with the Pluto notebook as follows. Go to the root folder of the project and run the following command:
   ```
   julia --project=.
   ``` 
   Then, in the Julia REPL, run:
   ```julia
   import Pluto 
   Pluto.run(notebook="notebooks/BiotBrinkmanVorticityPreconditionerTests.jl") 
   ```
   this will trigger a web browser tab with the contents of the notebook. Finally, in the notebook, select the data directory with the results you would like to visualize using the _"Select commit ID data directory"_ drop-down list.
