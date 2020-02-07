#!/bin/bash -l
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=200M
#SBATCH --array=0-99

#An array job computed with the computation cluster. Bash variable $SLURM_ARRAY_TASK_ID contains the array job ID.

cd /scratch/work/hallav1/LuxUS/ROC/scripts

module load GSL/2.4-goolf-triton-2017a

MU_B_IDENTIFIER=mu_minus1_4_2_3

RADMETHPATH=/path/to/methpipe/tool

OUTPUTFOLDER=/path/to/output/folder/"$MU_B_IDENTIFIER"
INPUTFOLDER=/path/to/input/folder/"$MU_B_IDENTIFIER"
DESIGNFOLDER=/path/to/design/matrix/folder/"$MU_B_IDENTIFIER"
TIMERESULTFOLDER=/path/to/output/folder/"$MU_B_IDENTIFIER"
TIMEFILEIDENTIFIER="radmeth_runtimes"
WIDTH=1000
echo "$SLURM_ARRAY_TASK_ID"

for REPS in 3 6 12
do
        for READS in 6 12 24 
        do
                for DIFF in 0 1
                do
                        START=$(date +%s.%N)
                        "$RADMETHPATH"/methpipe-3.4.3/bin/radmeth regression -factor IsCase "$DESIGNFOLDER"/design_matrix_C10_Q"$READS"_R"$REPS"_W1000_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt "$INPUTFOLDER"/proportion_table_C10_Q"$READS"_R"$REPS"_W1000_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt > "$OUTPUTFOLDER"/cpgs_Q"$READS"_R"$REPS"_diff"$DIFF"_"$SLURM_ARRAY_TASK_ID".bed
                        "$RADMETHPATH"/methpipe-3.4.3/bin/radmeth adjust -bins 1:200:1 "$OUTPUTFOLDER"/cpgs_Q"$READS"_R"$REPS"_diff"$DIFF"_"$SLURM_ARRAY_TASK_ID".bed > "$OUTPUTFOLDER"/cpgs_Q"$READS"_R"$REPS"_diff"$DIFF"_"$SLURM_ARRAY_TASK_ID".adjusted.bed
                        RUNTIME=$(echo "$(date +%s.%N) - $START" | bc)
                        printf "%s\t0\tradmeth\n" "$RUNTIME" >> "$TIMERESULTFOLDER"/"$TIMEFILEIDENTIFIER"_C10_Q"$READS"_R"$REPS"_D"$DIFF"_l38_w"$WIDTH"_"$SLURM_ARRAY_TASK_ID".txt

                done
        done
done
