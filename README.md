# LuxUS
*Lux*GLM *U*sing *S*patial correlation (LuxUS) is a tool for differential methylation analysis. The tool is based on generalized linear mixed model with spatial correlation structure. The model parameters are fitted using probabilistic programming language Stan [2]. Savage-Dickey Bayes factor estimates are used for statistical testing of a covariate of interest. LuxUS supports both continuous and binary variables. The model takes into account the experimental parameters, such as bisulfite conversion efficiency.

The needed Python scripts and Stan model files for running the tool and example input and output files are stored in this GitHub repository.

## Outline
* Requirements
* Running preanalysis
* Running LuxUS analysis
* Combining results files
* Simulating data from the LuxUS model
* References

## Requirements
- Python 2.7
- NumPy
- SciPy
- Matplotlib
- PyStan [3]
- CmdStan [4] (for running ADVI)

The versions used were: NumPy 1.14.5, SciPy 1.1.0, Matplotlib 2.2.2, PyStan 2.17.1.0, CmdStan 2.12.0. The tool has been tested in Linux environment.

## Running preanalysis

The script *prepare_data_for_luxus_publishing.py* can be used to prepare BS-seq data for LuxUS analysis. Before LuxUS analysis the data should be aligned and the methylation proportions for each cytosine in each sample should be calculated. The input file for a BS-seq experiment with N samples in total should be a headerless tab-separated proportion table with the following format

```
<Chromosome name> <Genomic location start> <Genomic location end> <Total count, sample 1> <Methylated read count, sample 1> ... <Total count, sample N> <Methylated read count, sample N>
```
The proportion table should contain cytosines from one chromosome only and be sorted based on the genomic location of the cytosines. If available, the experimental parameters bisulfite conversion efficiency, incorrect bisulfite conversion efficiency and sequencing error rates should also be provided to the script.

The corresponding design matrix for the experiment should also be given as a headerless text file. The rows of the file correspond to the samples, which should be given in the same order with the input proportion table. The first column represents the intercept in the model, and should be set to all ones. The next columns correspond to the other covariates in the model. By default, the covariate to be tested is the second column (column 1, as the numbering starts from 0). This covariate is used as a test covariate in the F-test in the preanalysis filtering step. The argument *-u REQ_PVAL* defines the p-value cutoff-value for the F-test. This covariate is assumed to be a binary variable, although both binary and continuous variables can be tested with the LuxUS model. The variables in the design matrix should be either binary or continuous.  

In the preanalysis phase each cytosine is first evaluated separately. The cytosine has to be located at maximum in the defined window width range from the start of the window. There also has to be samples with minimum coverage of *-v REQUIRED_COVERAGE* at least in as many samples as defined by *-s N_REQUIRED_SAMPLES* in both case and control groups, which are defined by the covariate *-t T_COVARIATE*. The sorted cytosines are evaluated one by one, until either the maximum width of a genomic window or maximum number of cytosines *-c N_CYTOSINES* is reached. Then the genomic window goes through another filtering step. The mean coverage over the window should be at least *-v REQUIRED_COVERAGE* in at least as many samples as defined by argument *-s N_REQUIRED_SAMPLES* in both case and control groups, which are again defined by covariate *-t T_COVARIATE*. Then the F-test is performed using the log-transformed average methylation states of the samples as data. 

For the windows that pass the preanalysis step, the script produces input files for the *run_luxus.py* script. Each input file contains the data of one genomic window. The number of cytosines and mean coverages for each window can be stored into separate text files. The information of whether a cytosine passed the preanalysis step and the index of the possible genomic window into which it belongs is stored into file with name *INPUT_NAME_in_analysis_indicator.txt*, where *INPUT NAME* is the name of the proportion table file. This file is later used when combining the calculated Bayes factors and the original input file into the final results file.

```
usage: prepare_data_for_luxus.py [-h] -i INPUT_NAME -o OUTPUT_FOLDER -d
                                 DESIGN_MATRIX [-w WIDTH] [-c N_CYTOSINES] -r
                                 N_REPLICATES [-b SIGMAB2] [-p ALPHA_L]
                                 [-q BETA_L] [-m ALPHA_E] [-n BETA_E]
                                 [-k ALPHA_R] [-l BETA_R] [-g ALPHA_C]
                                 [-f BETA_C] [-s N_REQUIRED_SAMPLES]
                                 [-a BSEFF] [-e BSBEFF] [-j SEQERR]
                                 [-t T_COVARIATE] -u REQ_PVAL
                                 [-v REQUIRED_COVERAGE] [-y MEANCOVFILE]
                                 [-z CYTNFILE]

Takes in the data in text file format and then prepares it for running LuxUS.
The data objects will be saved to the specified location.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_NAME, --BS_data_file_name INPUT_NAME
                        The input data file. The input data file should
                        contain only cytosines in one chromosome.
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
                        not specified, default value 5 is used.
  -n BETA_E, --beta_e BETA_E
                        Hyperparameter beta for prior distribution of e. If
                        not specified, default value 5 is used.
  -k ALPHA_R, --alpha_r ALPHA_R
                        Hyperparameter alpha for prior distribution of
                        sigmaR2. If not specified, default value 98 is used.
  -l BETA_R, --beta_r BETA_R
                        Hyperparameter beta for prior distribution of sigmaR2.
                        If not specified, default value 143 is used.
  -g ALPHA_C, --alpha_c ALPHA_C
                        Hyperparameter alpha for prior distribution of
                        sigmaC2. If not specified, default value 6 is used.
  -f BETA_C, --beta_c BETA_C
                        Hyperparameter beta for prior distribution of sigmaC2.
                        If not specified, default value 3 is used.
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
  -u REQ_PVAL, --required_pval REQ_PVAL
                        Required maximum p-value for the coefficient for the
                        variable given as --test_covariate for the cytosine
                        window to be included in the analysis. Set to 1 if no
                        p-value restriction is desired.
  -v REQUIRED_COVERAGE, --required_coverage REQUIRED_COVERAGE
                        Required average coverage over window for a replicate.
                        If not specified, default value 5 is used.
  -y MEANCOVFILE, --meanCovFile MEANCOVFILE
                        File name for storing mean coverages for each window
  -z CYTNFILE, --cytNFile CYTNFILE
                        File name for storing the number of cytosines in each
                        window.
```

Test input data files and all the result files have been provided in the *data* folder in this GitHub repository. The command for running the *prepare_data_for_luxus.py* script (with default parameters and the p-value limit for the preanalysis F-test set to 0.1) for the test input file *proportion_table_test_data_diff1.txt* and corresponding design matrix *design_matrix_test_data_diff1.txt* (both stored in folder defined in $INPUT_FOLDER) could be for example the following
```
python prepare_data_for_luxus.py -i "$INPUT_FOLDER"/proportion_table_test_data_diff1.txt -d "$INPUT_FOLDER"/design_matrix_test_data_diff1.txt -o $OUTPUT_FOLDER -r 12 -t 1 -u 0.1 -y "$OUTPUT_FOLDER"/window_mean_coverage_test_data_diff1.txt -z "$OUTPUT_FOLDER"/window_number_of_cytosines_test_data_diff1.txt
```
This command produces the following output files
- Stan input file *input_for_luxus_1.txt*. As the test data covers only 10 cytosines in a 1000bp region, only one genomic window was found, resulting in one Stan input file. The file is stored in the folder defined in *$OUTPUT_FOLDER*. 
- Window index file *proportion_table_test_data_diff1_in_analysis_indicator.txt*, where an indicator value is stored for each cytosine. Value 0 means that the cytosine was not included in a genomic window (or the window did not pass the preanalysis phase) and other values indicate the genomic window into which the cytosine belongs to. 
- *window_mean_coverage_test_data_diff1.txt* which contains the mean coverages for each sample (over the genomic window) for each genomic window that passed the preanalysis filtering phase. 
- *window_number_of_cytosines_test_data_diff1.txt* which contains the number of cytosines in the genomic window for each genomic window that passed the preanalysis filtering phase. 

The test input file *proportion_table_test_data_diff0.txt* contains a test data set in which no differential methylation is present and it will not pass the preanalysis filtering.

## Running LuxUS analysis

The *run_luxus.py* script can be used to run the LuxUS analysis for the desired input data, which has first been prepared with the *prepare_data_for_LuxUS.py* script. The script loads the input data, fits the model with using the chosen model fitting method, calculates the Bayes factor and saves it into the result file. The computation time can also be saved if desired. As the script is run on one genomic window at a time, it is possible to parallelize the computation of Bayes factors if desired. 

With argument *-a ALGORITHM* the algorithm to be used for fitting the model parameter can be chosen between Hamiltonian Monte Carlo (HMC) sampling and Automatic Differentation Variational Inference (ADVI), with Hamiltonian Monte Carlo sampling being the default option. When performing model fitting with HMC, the standard Stan fit summary can be found from the Python log file. If desired, the sample chains and histograms of the samples for the variance parameters for the cytosine and replicate random effects and the noise term, fixed effect coefficients and lengthscale parameter can be plotted. When using ADVI for model fitting, the input files are first saved as temporal [R dump files](https://pystan.readthedocs.io/en/latest/conversion.html). Then CmdStan is called using [subprocess package](https://docs.python.org/2/library/subprocess.html) in Python and temporal sample files are produced as a result. The samples are extracted from these files and used for Bayes factor calculation. After this the temporal files are removed by the script. Before running this script with ADVI chosen as the model fitting method, the Stan model should be built using CmdStan. This can be done by running [*make* command](https://github.com/stan-dev/cmdstan/wiki/Getting-Started-with-CmdStan) for the Stan model file (.stan file) in the CmdStan folder. The result will be stored in the same folder where the original .stan file was located. This should be the same folder from which the *run_luxus.py* script is run. 

The script uses a precompiled Stan model if possible by pickling the compiled Stan model after first usage. Please see the detailed description on the [PyStan manual page](https://pystan.readthedocs.io/en/latest/avoiding_recompilation.html).

```
usage: run_luxus.py [-h] -d INPUT_DATA -o OUTPUTFOLDER -i INPUTFOLDER -j
                    OUTPUTFILE -x TEST_COVARIATE [-y TEST_COVARIATE2]
                    [-b SIGMAB2] -a ALGORITHM [-p DIAGNOSTIC_PLOTS]
                    [-g N_GRADSAMPLES] [-e N_ELBOSAMPLES]
                    [-v N_OUTPUTSAMPLES_VI] [-m N_OUTPUTSAMPLES_HMC]
                    [-c N_HMC_CHAINS] [-w WINDOW_INDEX] [-t TIMEFILE]

Runs LuxUS model for the given input data and returns a BF for the whole
window.

optional arguments:
  -h, --help            show this help message and exit
  -d INPUT_DATA, --input_data INPUT_DATA
                        Name of the data file.
  -o OUTPUTFOLDER, --outputFolder OUTPUTFOLDER
                        Folder where to store the results.
  -i INPUTFOLDER, --inputFolder INPUTFOLDER
                        Folder where the input data is stored.
  -j OUTPUTFILE, --outputFile OUTPUTFILE
                        File into which the BFs are written. Will be located
                        in folder specified in -o.
  -x TEST_COVARIATE, --test_covariate TEST_COVARIATE
                        Covariate to be tested. Give index (in design matrix)
                        starting from 0.
  -y TEST_COVARIATE2, --test_covariate2 TEST_COVARIATE2
                        Type 2 test: the covariate to be compared to the
                        covariate defined by argument -x. If not provided,
                        type 1 test will be performed. Give index (in design
                        matrix) starting from 0.
  -b SIGMAB2, --sigmaB2 SIGMAB2
                        Variance for B. Default value 15 is used if not
                        specified.
  -a ALGORITHM, --algorithm ALGORITHM
                        Give value 0 (use HMC, default) or 1 (use VI).
  -p DIAGNOSTIC_PLOTS, --diagnostic_plots DIAGNOSTIC_PLOTS
                        Give value 0 (do not plot sample diagnostics for HMC)
                        or 1 (plot sample diagnostics for HMC). Default value
                        is 0.
  -g N_GRADSAMPLES, --N_gradsamples N_GRADSAMPLES
                        Number of gradient samples used in VI. Default value
                        10 is used if not specified.
  -e N_ELBOSAMPLES, --N_elbosamples N_ELBOSAMPLES
                        Number of ELBO samples used in VI. Default value 200
                        is used if not specified.
  -v N_OUTPUTSAMPLES_VI, --N_outputsamples_VI N_OUTPUTSAMPLES_VI
                        Number of posterior samples used in VI. Default value
                        2000 is used if not specified.
  -m N_OUTPUTSAMPLES_HMC, --N_outputsamples_HMC N_OUTPUTSAMPLES_HMC
                        Number of posterior samples per chain used in HMC (the
                        burn-in will be removed from this sample number).
                        Default value 1000 is used if not specified.
  -c N_HMC_CHAINS, --N_HMC_chains N_HMC_CHAINS
                        Number of chains in HMC sampling. Default value 4 is
                        used if not specified.
  -w WINDOW_INDEX, --window_index WINDOW_INDEX
                        The index of the window being analysed. If value is
                        not given the BF is saved without window index into
                        the defined output file.
  -t TIMEFILE, --timeFile TIMEFILE
                        File name (and path) for storing computation time. If
                        no file name is given the computation times will not
                        be stored into a file.

```

Continuing the analysis of the test input data set, the produced Stan input file *input_for_luxus_1.txt* can now be given as input to the *run_luxus.py* script. The following command will run LuxUS analysis (with type 1 test) on the input data with default parameters and by saving the diagnostics plot
```
python run_luxus.py -d input_for_luxus_1.txt -o $OUTPUT_FOLDER -i $OUTPUT_FOLDER -j test_data_diff1_result.txt -x 1 -p 1 -w 1
```
This will produce an output file *test_data_diff1_result.txt* (stored in folder defined by argument *-o OUTPUTFOLDER*), where the computed Bayes factor for the input data is stored along with the window index defined by argument *-w WINDOW_INDEX*. The diagnostics plot *input_for_luxus_1_diagnostic_plots_HMC.png* is also saved in the same folder. Input file name (the Stan input file) should be given as argument *-d INPUT_DATA* and the input file folder as argument *-i INPUTFOLDER*. 

 
## Combining results files

If desired, the script *final_results.py* can be used to combine the Bayes factor result file and the original proportion table file into a final result file. This script produces a new column, in which either the calculated Bayes factor or a "\*\" (for cytosines which were filtered in the preanalysis phase) is presented for each cytosine. The needed input files are the proportion table file, which was given as an input for the *prepare_data_for_luxus.py* script, the file into which the Bayes factor values and corresponding window indices were stored and the window index file, which is an output file from the *prepare_data_for_luxus.py* script. A folder in which the resulting file should be stored should be specified. The results file can then be used and filtered according to the user's needs. 

```
usage: final_results.py [-h] -i INPUT_FILE [-f OUTPUT_FOLDER] -o OUTPUT_FILE
                        -b BF_FILE -w WINDOW_INDEX

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

Again, continuing the analysis of the test input data, we combine the input file and the calculated Bayes factors into a final results file with command
```
python  final_results.py -i "$INPUT_FOLDER"/proportion_table_test_data_diff1.txt  -f $OUTPUT_FOLDER -o test_data_diff1_final_results.txt -b "$OUTPUT_FOLDER"/test_data_diff1_result.txt -w "$OUTPUT_FOLDER"/proportion_table_test_data_diff1_in_analysis_indicator.txt
```
The script produces the file *test_data_diff1_final_results.txt*. 


## Simulating data from the LuxUS model

The script *simulate_data_LuxUS.py* creates *-n N_DATASETS* datasets with differential methylation and *-n N_DATASETS* with no differential methylation between the case and control samples, each with *-c N_CYTOSINES* cytosines in them. The locations for the cytosines are generated randomly from 1 to *-w WIDTH*. The simulation function uses experimental design with an intercept term and case/control indicator variable, and the coefficients *b* are defined using parameter *-m MEAN_B*. The example command would create data with total read count *-q READS* for every cytosine. The number of samples simulated would be 6 + 6 (case replicates + control replicates), if 6 is given as parameter *-r REPLICATES*. The script saves the generated data sets to the folder defined by *-f FOLDER*. The user can define the formats in which the data is saved by using parameters *-x SAVE_LUXUS*, *-y SAVE_PROPORTION_TABLE* and *-z SAVE_LUXUS_SEP*. When specifying *-x 1* and/or *-z 1*, the data is saved in a format that can be given as input to the LuxUS Stan model. The parameter *-j SIGMAB2_STAN* defines the variance of the coefficients *b* used in the LuxUS analysis. If the data is saved only in proportion table format (*-x 0 -y 1 -z 0*), the parameter *-j* is not needed and arbitrary value can be given as input to the parameter. Defining *-z 1* saves data for every cytosine separately in a format that can be used with LuxUS Stan model.

```
usage: simulate_data_LuxUS.py [-h] -e SIGMAE2 -q READS -r REPLICATES -c
                              N_CYTOSINES -l L -w WIDTH -n N_DATASETS -f
                              FOLDER -g ID -v SIGMAB2 -d SIGMAR2 -a SIGMAC2 -m
                              MEAN_B -j SIGMAB2_STAN -x SAVE_LUXUS -y
                              SAVE_PROPORTION_TABLE -z SAVE_LUXUS_SEP

Simulates bisulphite sequencing data from the LuxUS model.

optional arguments:
  -h, --help            show this help message and exit
  -e SIGMAE2, --sigmae2 SIGMAE2
                        Variance for residual term e.
  -q READS, --reads READS
                        Number of BS calls i.e. (total) reads.
  -r REPLICATES, --replicates REPLICATES
                        Number of replicates for cases and controls, the total
                        number of samples will be 2x replicates.
  -c N_CYTOSINES, --cytosines N_CYTOSINES
                        Number of cytosines.
  -l L, --lengthscale L
                        Lengthscale parameter for covariance matrix for the
                        simulations.
  -w WIDTH, --width WIDTH
                        Width of the genomic region where the N_cytosines lie
                        (in basepairs).
  -n N_DATASETS, --datasets N_DATASETS
                        Number of datasets to be created for both cases (no
                        differential methylation, differential methylation).
                        Total number of resulting datasets is 2xN_datasets.
  -f FOLDER, --folder FOLDER
                        Folder where to store the generated data files. Format
                        /level1/level2
  -g ID, --ID ID        Identifier used in the resulting file names.
  -v SIGMAB2, --sigmaB2 SIGMAB2
                        Variance for B, which is inserted into covariance
                        matrix diagonal (SigmaB). Used for simulations.
  -d SIGMAR2, --sigmar2 SIGMAR2
                        sigma_r2 used for simulations.
  -a SIGMAC2, --sigmac2 SIGMAC2
                        sigma_c2 used for simulations.
  -m MEAN_B, --mean_B MEAN_B
                        Mean for B0 (coefficient for the intercept term) and
                        value for B1 (coefficient for case/control covariate)
                        for the cases with differential methylation. Should be
                        given as string in format [B1, B2]
  -j SIGMAB2_STAN, --sigmaB2_stan SIGMAB2_STAN
                        Variance for B, the value which will be used in LuxUS
                        analysis. (Can be different from what is used for
                        simulations).
  -x SAVE_LUXUS, --save_LuxUS SAVE_LUXUS
                        0 or 1 indicating whether the seimulated data should
                        be saved in a format supported by LuxUS.
  -y SAVE_PROPORTION_TABLE, --save_proportion_table SAVE_PROPORTION_TABLE
                        0 or 1 indicating whether the seimulated data should
                        be saved in proportion table format, which can be used
                        with eg. Radmeth. Also a design matrix is saved into a
                        text file.
  -z SAVE_LUXUS_SEP, --save_LuxUS_sep SAVE_LUXUS_SEP
                        0 or 1 indicating whether the seimulated data should
                        be saved in a format supported by LuxUS analysis where
                        each cytosine is analysed separately.
```

Usage example for the simulation function. This command would create 100 + 100 datasets with 6 + 6 samples. The total read count is 12 for each 10 cytosine in each data set. The simulated data is saved in folder /folder1/folder2 in LuxUS input format. The mean of the normal distribution from which the intercept term coefficients are drawn is -1.4 and the variance is set to 0.25. When simulating a data set with differential methylation, the coefficient for the case/control variable is set to 0.9 and when simulating a data set with no differential methylation, the coefficient for the case/control variable is set to 0. The variance terms and the region width are set as described in [5].

```
FOLDER=/folder1/folder2
SIGMAE2=1
SIGMAR2=0.69
SIGMAC2=2
SIGMAB2=0.25
SIGMAB2_STAN=15
N_CYT=10
READS=12
REPS=6
N_DATASETS=100
L=38
WIDTH=1000

python simulate_data_LuxUS.py -e $SIGMAE2 -q $READS -r $REPS -c $N_CYT -l $L -w $WIDTH -n $N_DATASETS -f $FOLDER -g "testing" -v $SIGMAB2 -d $SIGMAR2 -a $SIGMAC2 -m '[-1.4, 0.9]' -j $SIGMAB2_STAN -x 1 -y 0 -z 0
```

## References

[1] Äijö, T., Yue, X., Rao, A., Lähdesmäki, H (2016) LuxGLM: a probabilistic covariate model for quantification of DNA methylation modifications with complex experimental designs. *Bioinformatics*, 32(17), i511-i519.

[2] Carpenter, B., Gelman, A., Hoffman,M. D., Lee, D., Goodrich, B., Betancourt, M., Brubaker, M., Guo, J., Li, P., Riddell, A. (2017) Stan: A probabilistic programming language. *Journal of Statistical Software*, 76(1).

[3] Stan Development Team (2018) CmdStan: the command-line interface to Stan, Version 2.18.0.   http://mc-stan.org

[4] Stan Development Team (2017) PyStan: the Python interface to Stan, Version 2.16.0.0.   http://mc-stan.org

[5] Halla-aho, V. and Lähdesmäki, H. (2019) LuxUS: Detecting differential DNA methylation using generalized linear mixed model with spatial correlation structure. doi: https://doi.org/10.1101/536722
