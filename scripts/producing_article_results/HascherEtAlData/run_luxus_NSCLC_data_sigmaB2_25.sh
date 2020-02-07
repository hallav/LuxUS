#!/bin/bash -l

#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --job-name=LuxUS_for_NSCLC_data_sigmaB2_25
#SBATCH --mail-type=ALL
#SBATCH --mail-user=viivi.halla-aho@aalto.fi
#SBATCH --array=1-500

cd /scratch/work/hallav1/LuxUS/ROC/scripts/

source /scratch/work/hallav1/uusipython/bin/activate

CHR=1

WINDOW=$(sed ''"$SLURM_ARRAY_TASK_ID"'q;d' /scratch/work/hallav1/LuxUS/ROC/data/NSCLC/preanalysis/chr"$CHR"/random_window_inds_N_cytosines_10_or_more_one_based_chr"$CHR".txt)

#let "WINDOW = $WINDOW_0base + 1 "
echo "${WINDOW}"

INPUT_FOLDER=/scratch/work/hallav1/LuxUS/ROC/data/NSCLC/preanalysis_sigmaB2_25/chr"$CHR"
OUTPUT_FOLDER=/scratch/work/hallav1/LuxUS/ROC/results/NSCLC/sigmaB2_25/chr"$CHR"
OUTPUT_FILE=NSCLC_data_chr"$CHR"_BF_sigmaB2_25.txt
TIMEFILE="$OUTPUT_FOLDER"/computation_times_chr"$CHR"_sigmaB2_25.txt

ALGORITHM=0
TEST_COV_IND=1
SIGMAB2=25
D_PLOTS=1

N_SAMPLES=1500

SIGMAR2_FILE="$OUTPUT_FOLDER"/sample_means_sigmaR2_chr"$CHR".txt
SIGMAC2_FILE="$OUTPUT_FOLDER"/sample_means_sigmaC2_chr"$CHR".txt
SIGMAE2_FILE="$OUTPUT_FOLDER"/sample_means_sigmaE2_chr"$CHR".txt


python run_luxus.py -d input_for_luxus_"$WINDOW".txt -o $OUTPUT_FOLDER -i $INPUT_FOLDER -j $OUTPUT_FILE -c $TEST_COV_IND -b $SIGMAB2 -a $ALGORITHM -p $D_PLOTS -m $N_SAMPLES -w $WINDOW -t $TIMEFILE -x $SIGMAR2_FILE -y $SIGMAC2_FILE -z $SIGMAE2_FILE

