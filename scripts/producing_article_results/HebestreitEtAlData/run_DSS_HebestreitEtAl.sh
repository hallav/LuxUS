#!/bin/sh -l
#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=run_DSS
#SBATCH --mail-type=ALL
#SBATCH --array=1-3000

cd /path/to/scripts

module load r/3.5.0-python-2.7.14

#The input arguments should be
#1: Folder to be set as environment
#2: Input file name (proportion table)
#3: Input file name (design matrix) not really needed!
#4: Output file name (+ path if different than environment)
#5: Number of columns in proportion table = number of replicates


DMRSTATUS="nonDMR"
INPUTFOLDER=/path/to/preprocessed/input/files
OUTPUTFOLDER=/path/to/results/folder/DSS/"$DMRSTATUS"
DESIGNMATRIX="/path/to/data/files/design_simple_radmeth.txt" #Not really needed in this run
RESULTF="DSS_results"
RUNTIMEF="runtimes"
TOTREPS=12

srun Rscript --no-save run_DSS_KleinHebestreit.R $INPUTFOLDER TEMP_WINDOW_proportion_table_KleinHebestreit_"$DMRSTATUS"_wHeader_"$SLURM_ARRAY_TASK_ID".txt "$DESIGNMATRIX" "$OUTPUTFOLDER"/"$RESULTF"_"$DMRSTATUS"_"$SLURM_ARRAY_TASK_ID" $TOTREPS

DMRSTATUS="DMR" 
OUTPUTFOLDER=/path/to/results/folder/DSS/"$DMRSTATUS"
srun Rscript --no-save run_DSS_KleinHebestreit.R $INPUTFOLDER TEMP_WINDOW_proportion_table_KleinHebestreit_"$DMRSTATUS"_wHeader_"$SLURM_ARRAY_TASK_ID".txt "$DESIGNMATRIX" "$OUTPUTFOLDER"/"$RESULTF"_"$DMRSTATUS"_"$SLURM_ARRAY_TASK_ID" $TOTREPS
