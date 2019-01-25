# LuxUS
A tool for differential methylation analysis. The tool is based generalized linear mixed model with spatial correlation structure, which is fitted using probabilistic programming language Stan. Savage-Dickey Bayes factor estimates are used for statistical testing of a covariate of interest. LuxUS supports both continuous and binary variables. The model takes into account the experimental parameters, such as bisulfite conversion efficiency.

The needed Python scripts and Stan model files for running the tool and example input and output files are stored in this GitHub repository.

## Outline
* Requirements
* Running preanalysis
* Running LuxUS analysis
* Combining results files
* References

## Requirements
- Python 2.7
- Numpy
- Scipy
- PyStan
- CmdStan (for running ADVI) https://github.com/stan-dev/cmdstan

The versions used were: Numpy 1.14.5, PyStan 2.17.1.0, CmdStan 2.12.0. The tool has been tested in Linux environment.

## Running preanalysis

The script *prepare_data_for_luxus_publishing.py* can be used to prepare BS-seq data for LuxUS analysis. Before LuxUS analysis the data should be aligned and the methylation proportions for each cytosine in each sample should be calculated. The input file for a BS-seq experiment with N samples in total should be a headerless tab-separated proportion table with the following format

```
<Chromosome name> <Start> <End> <Total count, sample 1> ... <Totanl count, sample N> <Methylation level, sample 1> ... <Methylation level, sample N>
```
The proportion table should contain cytosines from one chromosome only and be sorted based on the genomic location of the cytosines. If available, the experimental parameters bisulfite conversion efficiency, incorrect bisulfite conversion efficiency and sequencing error rates should also be provided to the script.

The corresponding design matrix for the experiment should also be given as a headerless text file. The rows of the file correspond to the samples, which should be given in the same order with the input proportion table. The first column represents the intercept in the model, and should be set to all ones. The next columns correspond to the other covariates in the model. By default, the covariate to be tested is the second column (column 1, as the numbering starts from 0). This covariate is used as a test covariate in the F-test in the preanalysis filtering step. The parameter *INSERT PARAMETER NAME* defines the p-value cutoff-value for the F-test. This covariate is assumed to be a binary variable, although both binary and continuous variables can be tested with the LuxUS model. The variables in the design matrix should be either binary or continuous.  

In the preanalysis each cytosine is first evaluated separately. The cytosine has to be located at maximum in the defined window width range from the start of the window. There also has to be samples with minimum coverage of *REQUIRED_COVERAGE* at least in as many samples as defined by *N_REQUIRED_SAMPLES* in both case and control groups, which are defined by the covariate *T_COVARIATE*. The sorted cytosines are evaluated one by one, until either the maximum width of a genomic window or maximum number of cytosines *N_CYTOSINES* is reached. Then the genomic window goes through the preanalysis step. The mean coverage over the window should be at least *REQUIRED_COVERAGE* in at least as many samples as defined by *N_REQUIRED_SAMPLES* in both case and control groups, which are again defined by the covariate *T_COVARIATE*. Then the F-test is performed using the log-transformed average methylation states of the samples as data. 

For the windows that pass the preanalysis step, the script produces input files for the *run_luxus.py* script. Each input file contains the data of one genomic window. The number of cytosines and mean coverages for each window are stored into separate text files. The information of whether a cytosine passed the preanalysis step and the index of the possible genomic window into which it belongs is stored into file *INSERT FILE NAME HERE*. This file is later used when combining the calculated Bayes factors and the original input file into the final result file.


```
usage: prepare_data_for_luxus_publishing.py [-h] -i INPUT_NAME -o
                                            OUTPUT_FOLDER -d DESIGN_MATRIX
                                            [-w WIDTH] [-c N_CYTOSINES] -r
                                            N_REPLICATES [-b SIGMAB2]
                                            [-p ALPHA_L] [-q BETA_L]
                                            [-m ALPHA_E] [-n BETA_E]
                                            [-k ALPHA_R] [-l BETA_R]
                                            [-g ALPHA_C] [-f BETA_C]
                                            [-s N_REQUIRED_SAMPLES] [-a BSEFF]
                                            [-e BSBEFF] [-j SEQERR]
                                            [-t T_COVARIATE]
                                            [-u REQ_METH_DIFF]
                                            [-v REQUIRED_COVERAGE] -y
                                            MEANCOVFILE -z CYTNFILE

Takes in the data in text file format and then prepares it for running LuxUS.
The data objects will be saved to the specified location.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_NAME, --BS_data_file_name INPUT_NAME
                        The input data file.
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        The output location.
  -d DESIGN_MATRIX, --design_matrix DESIGN_MATRIX
                        The design matrix file (and file path if needed),
                        where there is one row for each sample. Include also
                        the column of ones corresponding to the intercept.
  -w WIDTH, --window_width WIDTH
                        Maximum width of the genomic window in basepairs. If
                        not specified, default value 1000 is used.
  -c N_CYTOSINES, --N_cytosines N_CYTOSINES
                        Maximum number of cytosines in a window. If not
                        specified 20 is used as default.
  -r N_REPLICATES, --N_replicates N_REPLICATES
                        Number of replicates.
  -b SIGMAB2, --sigmab2 SIGMAB2
                        Prior variance for the coefficients b. If not
                        specified, default value 15 is used.
  -p ALPHA_L, --alpha_l ALPHA_L
                        Hyperparameter alpha for prior distribution of l. If
                        not specified, default value 38 is used.
  -q BETA_L, --beta_l BETA_L
                        Hyperparameter beta for prior distribution of l. If
                        not specified, default value 1 is used.
  -m ALPHA_E, --alpha_e ALPHA_E
                        Hyperparameter alpha for prior distribution of e. If
                        not specified, default value 2 is used.
  -n BETA_E, --beta_e BETA_E
                        Hyperparameter beta for prior distribution of e. If
                        not specified, default value 2 is used.
  -k ALPHA_R, --alpha_r ALPHA_R
                        Hyperparameter alpha for prior distribution of
                        sigmaR2. If not specified, default value 2 is used.
  -l BETA_R, --beta_r BETA_R
                        Hyperparameter beta for prior distribution of sigmaR2.
                        If not specified, default value 2 is used.
  -g ALPHA_C, --alpha_c ALPHA_C
                        Hyperparameter alpha for prior distribution of
                        sigmaC2. If not specified, default value 2 is used.
  -f BETA_C, --beta_c BETA_C
                        Hyperparameter beta for prior distribution of sigmaC2.
                        If not specified, default value 2 is used.
  -s N_REQUIRED_SAMPLES, --N_required_samples N_REQUIRED_SAMPLES
                        Number of samples (from both case and control groups)
                        each cytosine must be present for it to be included in
                        the analysis. If not specified, default value 1 is
                        used.
  -a BSEFF, --bsEff BSEFF
                        Bisulphite conversion efficiency in format
                        [bsEff_1,...,bsEff_N_replicates]. If not specified,
                        default value 1 is used for all samples.
  -e BSBEFF, --bsbEff BSBEFF
                        Incorrect bisulphite conversion efficiency in in
                        format [bsbEff_1,...,bsbEff_N_replicates]. If not
                        specified, default value 0 is used for all samples.
  -j SEQERR, --seqErr SEQERR
                        Sequencing error in format
                        [seqErr_1,...,seqErr_N_replicates]. If not specified,
                        default value 0 is used for all samples.
  -t T_COVARIATE, --test_covariate T_COVARIATE
                        Give integer value indicating the column which is the
                        covariate to be tested. Assumed to be a binary
                        covariate. Indexing starts from 0. If not specified,
                        default value 1 is used.
  -u REQ_METH_DIFF, --required_meth_diff REQ_METH_DIFF
                        Required p-value for the coefficient for the variable
                        given as --test_covariate for the cytosine window to
                        be included in the analysis. Set to 1 if no p-value
                        restriction is desired.
  -v REQUIRED_COVERAGE, --required_coverage REQUIRED_COVERAGE
                        Required average coverage over window for a replicate.
                        If not specified, default value 5 is used.
  -y MEANCOVFILE, --meanCovFile MEANCOVFILE
                        File name for storing mean coverages for each window
  -z CYTNFILE, --cytNFile CYTNFILE
                        File name for storing the number of cytosines in each
                        window.
```
## Running LuxUS analysis

The *run_luxus.py* script can be used to run the LuxUS analysis for the desired input data, which has first been prepared with the *prepare_data_for_LuxUS.py* script. The script loads the input data, fits the model with using the chosen model fitting method, calculates the Bayes factor and saves it into the result file. The computation time can also be saved if desired. As the script is run on one genomic window at a time, it is possible to parallelize the computation of Bayes factors if desired. 

With parameter *-a* the algorithm to be used for fitting the model parameter can be chosen, with Hamiltonian Monte Carlo sampling being the default option. When performing model fitting with HMC, the standard Stan fit summary can be found from the Python log file. If desired, the sample chains and histograms of the samples for the variance parameters for the cytosine and replicate random effects and the noise term, fixed effect coefficients and lengthscale parameter can be plotted. When using ADVI for model fitting, the input files are first saved as temporal [R dump files](https://pystan.readthedocs.io/en/latest/conversion.html). Then CmdStan is called using [subprocess package](https://docs.python.org/2/library/subprocess.html) in Python and temporal sample files are produced as a result. The samples are extracted from these files and used for Bayes factor calculation. After this the temporal files are removed by the script. Before running this script with ADVI chosen as the model fitting method, the Stan model should be built using CmdStan. This can be done by running [*make* command](https://github.com/stan-dev/cmdstan/wiki/Getting-Started-with-CmdStan) for the Stan model file (.stan file) in the CmdStan folder. The result will be stored in the same folder where the original .stan file was located. This should be the same folder from which the *run_luxus.py* script is run. 

The script uses a precompiled Stan model if possible by pickling the Stan model. Please see the description on the [PyStan manual page](https://pystan.readthedocs.io/en/latest/avoiding_recompilation.html).

```
usage: run_luxus.py [-h] -d INPUT_DATA -o OUTPUTFOLDER -i INPUTFOLDER -j
                    OUTPUTFILE -c TEST_COVARIATE [-b SIGMAB2] -a ALGORITHM
                    [-p DIAGNOSTIC_PLOTS] [-g N_GRADSAMPLES]
                    [-e N_ELBOSAMPLES] [-v N_OUTPUTSAMPLES_VI]
                    [-m N_OUTPUTSAMPLES_HMC] [-w WINDOW_INDEX] [-t TIMEFILE]

Runs LuxUS model for the given input data and returns a BF for the whole
window.

optional arguments:
  -h, --help            show this help message and exit
  -d INPUT_DATA, --input_data INPUT_DATA
                        Name of the data file.
  -o OUTPUTFOLDER, --outputFolder OUTPUTFOLDER
                        Folder where to store the results. Format
                        /level1/level2
  -i INPUTFOLDER, --inputFolder INPUTFOLDER
                        Folder where the input data is stored. Format
                        /level1/level2
  -j OUTPUTFILE, --outputFile OUTPUTFILE
                        File into which the BFs are written. Will be located
                        in folder specified in -o.
  -c TEST_COVARIATE, --test_covariate TEST_COVARIATE
                        Covariate to be tested. Give index (in design matrix)
                        starting from 0.
  -b SIGMAB2, --sigmaB2 SIGMAB2
                        Variance for B. Default value is used if not
                        specified.
  -a ALGORITHM, --algorithm ALGORITHM
                        Give value 0 (use HMC) or 1 (use VI).
  -p DIAGNOSTIC_PLOTS, --diagnostic_plots DIAGNOSTIC_PLOTS
                        Give value 0 (do not plot sample diagnostics for HMC)
                        or 1 (plot sample diagnostics for HMC). Default value
                        is 0.
  -g N_GRADSAMPLES, --N_gradsamples N_GRADSAMPLES
                        Number of gradient samples used in VI. Default value
                        is used if not specified.
  -e N_ELBOSAMPLES, --N_elbosamples N_ELBOSAMPLES
                        Number of gradient samples used in VI. Default value
                        is used if not specified.
  -v N_OUTPUTSAMPLES_VI, --N_outputsamples_VI N_OUTPUTSAMPLES_VI
                        Number of posterior samples used in VI. Default value
                        is used if not specified.
  -m N_OUTPUTSAMPLES_HMC, --N_outputsamples_HMC N_OUTPUTSAMPLES_HMC
                        Number of posterior samples per chain used in HMC.
                        Default value is used if not specified.
  -w WINDOW_INDEX, --window_index WINDOW_INDEX
                        The index of the window being analysed. If value is
                        not given the BF is saved without window index into
                        the defined output file.
  -t TIMEFILE, --timeFile TIMEFILE
                        File name (and path) for storing computation time. If
                        no file name is given the computation times will not
                        be stored into a file.
```


## Combining results files

If desired, the script *final_results.py* can be used to combine the Bayes factor result file and the original proportion table file into a final result file. This script produces a new column, in which either the calculated Bayes factor or a "\*\" is presented for each cytosine. The needed input files are the proportion table file, which was given as an input for the *prepare_data_for_luxus.py* script, the file into which the Bayes factor values and corresponding window indices were stored and the window index file, which is an output file from the *prepare_data_for_luxus.py* script. A folder in which the resulting file should be stored should be specified.

```
usage: final_results_publishing.py [-h] -i INPUT_FILE [-f OUTPUT_FOLDER] -o
                                   OUTPUT_FILE -b BF_FILE -w WINDOW_INDEX

Generates a final results file where the calculated Bayes factors are included
as columns into the original data file.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --BS_data_file_name INPUT_FILE
                        The input data file name.
  -f OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        The output folder where the results files will be
                        stored.
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Name for the final output file.
  -b BF_FILE, --BF_file BF_FILE
                        The BF file name.
  -w WINDOW_INDEX, --window_index_file WINDOW_INDEX
                        The window index file.
```
## References
