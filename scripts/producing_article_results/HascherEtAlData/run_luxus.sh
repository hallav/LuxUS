#!/bin/bash -l
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --job-name=LuxUS_for_NSCLC_data
#SBATCH --mail-type=ALL
#SBATCH --array=1-500

cd /path/to/scripts/folder

source /path/to/virtual/environment/bin/activate

CHR=4

WINDOW=$(sed ''"$SLURM_ARRAY_TASK_ID"'q;d' /path/to/preanalysis/results/chr"$CHR"/random_window_inds_N_cytosines_10_or_more_one_based_chr"$CHR".txt)

echo "${WINDOW}"

INPUT_FOLDER=/path/to/preanalysis/results/chr"$CHR"
OUTPUT_FOLDER=/path/to/LuxUS/results/chr"$CHR"
OUTPUT_FILE=NSCLC_data_chr"$CHR"_BF.txt
TIMEFILE="$OUTPUT_FOLDER"/computation_times_chr"$CHR".txt

ALGORITHM=0
TEST_COV_IND=1
SIGMAB2=15 #By changing this parameter LuxUS can be run with the different sigmab2 values. 5 and 25 were used.
D_PLOTS=1

N_SAMPLES=1500

SIGMAR2_FILE="$OUTPUT_FOLDER"/sample_means_sigmaR2_chr"$CHR".txt
SIGMAC2_FILE="$OUTPUT_FOLDER"/sample_means_sigmaC2_chr"$CHR".txt
SIGMAE2_FILE="$OUTPUT_FOLDER"/sample_means_sigmaE2_chr"$CHR".txt


python run_luxus.py -d input_for_luxus_"$WINDOW".txt -o $OUTPUT_FOLDER -i $INPUT_FOLDER -j $OUTPUT_FILE -c $TEST_COV_IND -b 15 -a $ALGORITHM -p $D_PLOTS -m $N_SAMPLES -w $WINDOW -t $TIMEFILE -x $SIGMAR2_FILE -y $SIGMAC2_FILE -z $SIGMAE2_FILE

