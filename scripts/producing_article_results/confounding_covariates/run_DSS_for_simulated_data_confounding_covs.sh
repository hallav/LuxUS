#!/bin/sh -l
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=500MB
#SBATCH --job-name=run_DSS_comparison
#SBATCH --mail-type=ALL
#SBATCH --mail-user=viivi.halla-aho@aalto.fi
#SBATCH --array=0-949

cd /path/to/scripts

module load Perl/5.28.0-GCCcore-7.3.0
module load sqlite/3.28.0
module load gcc/9.2.0-cuda-nvptx

module load r/3.6.1-python3

#The input arguments should be
#1: Folder to be set as environment
#2: Input file name (proportion table)
#3: Input file name (design matrix)
#4: Output file name (+ path if different than environment)
#5: Number of columns in proportion table = number of replicates


INPUTFOLDER=/path/to/generated/data
OUTPUTFOLDER=/path/to/results
DESIGNMATRIX="design_matrix"
PTABLE="proportion_table"
RESULTF="DSS_results"
RUNTIMEF="runtimes"
WIDTH=1000

for REPS in 6 12
do
    TOTREPS=$((REPS*2))

    for READS in 6 12 24
    do
        for DIFF in 0 #for DIFF=1 set --array=0-99
        do

            srun Rscript --no-save run_DSS_for_simulated_data_confounding_covs.R  $INPUTFOLDER "$PTABLE"_C10_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt "$DESIGNMATRIX"_C10_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt "$OUTPUTFOLDER"/"$RESULTF"_C10_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt $TOTREPS

        done
    done
done
