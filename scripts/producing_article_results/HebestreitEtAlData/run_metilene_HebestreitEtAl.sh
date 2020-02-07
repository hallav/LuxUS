#!/bin/sh -l
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --job-name=run_metilene
#SBATCH --mail-type=ALL

cd /path/to/metilene_v0.2-8

INPUTFOLDER=/path/to/preprocessed/metilene/input/folder
OUTPUTFOLDER=/scratch/work/hallav1/LuxUS/ROC/results/KleinHebestreit/fixed/metilene
MTABLE=metilene_input_allWindows_sorted_wHeader.txt
RESULTF1="metilene_results_mode1"
WIDTH=2000

./metilene_linux64 -M $WIDTH -m 1 -d 0 -a control -b case "$INPUTFOLDER"/"$MTABLE" > "$OUTPUTFOLDER"/"$RESULTF1"_allWindows.txt

