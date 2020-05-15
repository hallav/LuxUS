#!/bin/sh -l
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_metilene_comparison
#SBATCH --mail-type=ALL
#SBATCH --mail-user=viivi.halla-aho@aalto.fi
#SBATCH --array=0-949


cd /scratch/work/hallav1/metilene_v0.2-8

INPUTFOLDER=/path/to/generated/data
OUTPUTFOLDER=/path/to/results
MTABLE="metilene_input"
MREGIONS="metilene_input_regions"
RESULTF1="metilene_results_mode1"
RESULTF2="metilene_results_mode2"
RESULTF3="metilene_results_mode3"
RUNTIMEF="runtimes"
WIDTH=1000

for REPS in 6 12
do

    TOTREPS=$((REPS*2))

    for READS in 6 12 24
    do
        for DIFF in 0
        do

            ./metilene_linux64 -M 1000 -m 1 -d 0 -a control -b case "$INPUTFOLDER"/"$MTABLE"_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt > "$OUTPUTFOLDER"/"$RESULTF1"_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt

            ./metilene_linux64 -M 1000 -m 1 -d 0 -f 2 -a control -b case -B "$INPUTFOLDER"/"$MREGIONS"_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".bed "$INPUTFOLDER"/"$MTABLE"_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt > "$OUTPUTFOLDER"/"$RESULTF2"_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt

            ./metilene_linux64 -M 1000 -m 1 -d 0 -f 3 -a control -b case "$INPUTFOLDER"/"$MTABLE"_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt > "$OUTPUTFOLDER"/"$RESULTF3"_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt

        done
    done
done
