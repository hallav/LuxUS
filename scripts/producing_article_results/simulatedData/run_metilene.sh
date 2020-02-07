#!/bin/sh -l
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_bsseq_comparison
#SBATCH --mail-type=ALL
#SBATCH --array=0-99

cd /scratch/work/hallav1/metilene_v0.2-8

INPUTFOLDER=/path/to/simulated/data/mu_minus1_4_1
OUTPUTFOLDER=/path/to/results/folder/mu_minus1_4_1
MTABLE="metilene_input"
MREGIONS="metilene_input_regions"
RESULTF1="metilene_results_mode1"
RESULTF2="metilene_results_mode2"
RESULTF3="metilene_results_mode3"
RUNTIMEF="runtimes"
WIDTH=1000

for REPS in 3 6 12
do

    TOTREPS=$((REPS*2))

    for READS in 6 12 24
    do
        for DIFF in 0 1
        do

            ./metilene_linux64 -M 1000 -m 1 -d 0 -a control -b case "$INPUTFOLDER"/"$MTABLE"_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt > "$OUTPUTFOLDER"/"$RESULTF1"_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt

            ./metilene_linux64 -M 1000 -m 1 -d 0 -f 2 -a control -b case -B "$INPUTFOLDER"/"$MREGIONS"_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".bed "$INPUTFOLDER"/"$MTABLE"_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt > "$OUTPUTFOLDER"/"$RESULTF2"_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt

            ./metilene_linux64 -M 1000 -m 1 -d 0 -f 3 -a control -b case "$INPUTFOLDER"/"$MTABLE"_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt > "$OUTPUTFOLDER"/"$RESULTF3"_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt

        done
    done
done
