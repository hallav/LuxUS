# LuxUS
A tool for differential methylation analysis.

The needed Python scripts, Stan model files and example input and output files are stored in this GitHub repository.


## Outline
* Requirements
* Running preanalysis
* Running LuxUS analysis
* Combining results files
* References

## Requirements
- Python 2.7
- Numpy
- PyStan
- CmdStan (for running ADVI)

## Running preanalysis

prepare_data_for_luxus_publishing.py
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

run_luxus.py

```
usage: run_luxus.py [-h] -d INPUT_DATA -o OUTPUTFOLDER -i INPUTFOLDER -j
                    OUTPUTFILE -c TEST_COVARIATE [-b SIGMAB2] -a ALGORITHM
                    [-p DIAGNOSTIC_PLOTS] [-g N_GRADSAMPLES]
                    [-e N_ELBOSAMPLES] [-v N_OUTPUTSAMPLES_VI]
                    [-m N_OUTPUTSAMPLES_HMC] [-w WINDOW_INDEX] -x SIGMAR2_FILE
                    -y SIGMAC2_FILE -z SIGMAE2_FILE [-t TIMEFILE]

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
  -x SIGMAR2_FILE, --sigmaR2_file SIGMAR2_FILE
                        File for storing sigmaR2 posterior means. REMOVE FROM
                        FINAL VERSION.
  -y SIGMAC2_FILE, --sigmaC2_file SIGMAC2_FILE
                        File for storing sigmaC2 posterior means. REMOVE FROM
                        FINAL VERSION.
  -z SIGMAE2_FILE, --sigmaE2_file SIGMAE2_FILE
                        File for storing sigmaE2 posterior means. REMOVE FROM
                        FINAL VERSION.
  -t TIMEFILE, --timeFile TIMEFILE
                        File name (and path) for storing computation time. If
                        no file name is given the computation times will not
                        be stored into a file.
```

## Combining results files

final_results.py
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
