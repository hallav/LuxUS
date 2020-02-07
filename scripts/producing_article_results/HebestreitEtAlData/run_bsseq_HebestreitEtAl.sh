#!/bin/sh -l
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --job-name=run_bsseq_KleinHebestreit_allData
#SBATCH --mail-type=ALL


cd /path/to/scripts/scripts

module load r/3.5.0-python-2.7.14

#The input arguments should be
#1: Folder to be set as environment
#2: Input file name (proportion table)
#3: Input file name (design matrix) NOT ACTUALLY NEEDED
#4: Output file name (+ path if different than environment)
#5: Number of columns in proportion table = number of replicates


DMRSTATUS="DMR"
INPUTFOLDER=/path/to/preprocessed/input/file/folder
OUTPUTFOLDER=/path/to/results/folder/bsseq
DESIGNMATRIX="/path/to/input/data/folder/design_simple_radmeth.txt" #Not really needed in this run
PTABLE="KleinHebestreit_allWindows_sorted_wHeader_v3.txt"
RESULTF="bsseq_results_allWindows"
RUNTIMEF="runtimes"
TOTREPS=12

srun Rscript --no-save run_bsseq_KleinHebestreit.R $INPUTFOLDER "$PTABLE" "$DESIGNMATRIX" "$OUTPUTFOLDER"/"$RESULTF" $TOTREPS

