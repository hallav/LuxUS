#!/bin/sh -l
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --job-name=run_M3D_for_simulated_data
#SBATCH --mail-type=ALL
#SBATCH --mail-user=viivi.halla-aho@aalto.fi



cd /path/to/scripts

module load r/3.5.0-python-2.7.14

#The input arguments should be
#1: The input data file name identifier (with path!)
#2: The total number of simulated data files to be used as input.
#3: Output filename identifier (with path!)
#4: Number of cytosines
#5: Number of replicates (total number is double of this)
#6: Number of data files for diff=0
#7: Number of data files for diff=1

N_DATASETS=1000
NCYT=10
N_D0=950
N_D1=50

for REPS in 6 12
do

    for READS in 6 12 24
    do

        OUTPUTFOLDER=/path/to/results/M3D_output_C10_Q"$READS"_R"$REPS"_W1000
        INPUTFILE=/path/to/generated/data/proportion_table_C10_Q"$READS"_R"$REPS"_W1000_set

        srun Rscript --no-save  run_M3D_for_simulated_data_confounding_covs.R $INPUTFILE $N_DATASETS $OUTPUTFOLDER $NCYT $REPS $N_D0 $N_D1
        
    done
done
