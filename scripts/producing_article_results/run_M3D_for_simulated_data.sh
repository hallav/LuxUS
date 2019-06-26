#!/bin/sh -l

cd /scratch/work/hallav1/LuxUS/ROC/scripts

module load r/3.5.0-python-2.7.14


#THE LIST OF THE INPUT PARAMETERS BELOW IS NOT UP TO DATE
#The input arguments should be 
#1: The input data file name identifier (with path!)
#2: The number of simulated data files to be used as input.
#4: Output filename identifier (with path!)
#5: Maximum number of cytosines per window


INPUTFLDR=/path/to/input/folder
OUTPUTFLDR=/path/to/output/folder

N_DATASETS=100
NCYT=10

MU_B_IDENTIFIER=mu_1_4_minus1

for REPS in 3 6 12
do

    for READS in 6 12 24
    do

        OUTPUTIDENTIFIER="$OUTPUTFOLDER"/"$MU_B_IDENTIFIER"/M3D/M3D_C10_Q"$READS"_R"$REPS"_W1000
        INPUTIDENTIFIER="$INPUTFOLDER"/"$MU_B_IDENTIFIER"/proportion_table_C10_Q"$READS"_R"$REPS"_W1000_set

        srun Rscript --no-save run_M3D_for_simulated_data.R $INPUTIDENTIFIER $N_DATASETS $OUTPUTIDENTIFIER $NCYT $REPS
        
    done
done
