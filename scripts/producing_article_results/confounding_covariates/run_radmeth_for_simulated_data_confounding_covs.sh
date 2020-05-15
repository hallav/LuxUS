#!/bin/bash -l
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=200M
#SBATCH --job-name=run_radmeth_regression_comparison
#SBATCH --mail-type=ALL
#SBATCH --mail-user=viivi.halla-aho@aalto.fi
#SBATCH --array=0-949

cd /path/to/scripts

module load GSL/2.4-goolf-triton-2017a
OUTFOLDER=/path/to/results
CHUNKFOLDER=/path/to/generated/data
DESIGNFOLDER=/path/to/generated/data
TIMERESULTFOLDER=/path/to/results
TIMEF="runtimes"
WIDTH=1000
echo "$SLURM_ARRAY_TASK_ID"

for REPS in 6 12
do
        for READS in 6 12 24
        do
                for DIFF in 0
                do

                        START=$(date +%s.%N)
                        /scratch/cs/csb/users/hallav1/software/methpipe-3.4.3/bin/radmeth regression -factor IsCase "$DESIGNFOLDER"/design_matrix_binary_C10_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt "$CHUNKFOLDER"/proportion_table_C10_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt > "$OUTFOLDER"/cpgs_Q"$READS"_R"$REPS"_diff"$DIFF"_"$SLURM_ARRAY_TASK_ID".bed
                        /scratch/cs/csb/users/hallav1/software/methpipe-3.4.3/bin/radmeth adjust -bins 1:200:1 "$OUTFOLDER"/cpgs_Q"$READS"_R"$REPS"_diff"$DIFF"_"$SLURM_ARRAY_TASK_ID".bed > "$OUTFOLDER"/adjusted/cpgs_Q"$READS"_R"$REPS"_diff"$DIFF"_"$SLURM_ARRAY_TASK_ID".adjusted.bed
                        #RUNTIME=$SECONDS
                        RUNTIME=$(echo "$(date +%s.%N) - $START" | bc)
                        printf "%s\t0\tradmeth\n" "$RUNTIME" >> "$TIMERESULTFOLDER"/"$TIMEF"_C10_Q"$READS"_R"$REPS"_D"$DIFF"_l38_w"$WIDTH"_"$SLURM_ARRAY_TASK_ID".txt

                done
        done
done

