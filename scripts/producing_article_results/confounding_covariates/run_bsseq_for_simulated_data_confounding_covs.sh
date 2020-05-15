#!/bin/sh -l
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=run_bsseq_comparison
#SBATCH --mail-type=ALL
#SBATCH --mail-user=viivi.halla-aho@aalto.fi
#SBATCH --array=0-949

cd /scratch/work/hallav1/LuxUS/ROC/scripts

#module load GCCcore/6.3.0
#module load GSL/2.4-goolf-triton-2017a

module load r/3.5.0-python-2.7.14

#The input arguments should be
#1: Folder to be set as environment
#2: Input file name (proportion table)
#3: Output file name (+ path if different than environment)
#4: Number of columns in proportion table = number of replicates


INPUTFOLDER=/path/to/generated/data/
OUTPUTFOLDER=/path/to/results/
DESIGNMATRIX="design_matrix"
PTABLE="proportion_table"
RESULTF="bsseq_results"
RUNTIMEF="runtimes"
WIDTH=1000

for REPS in 6 12
do

    TOTREPS=$((REPS*2))

    for READS in 6 12 24
    do
        for DIFF in 0
        do

            srun Rscript --no-save run_bsseq_for_simulated_data.R $INPUTFOLDER "$PTABLE"_C10_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt "$OUTPUTFOLDER"/"$RESULTF"_C10_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF" $TOTREPS

        done
    done
done
