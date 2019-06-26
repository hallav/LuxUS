#!/bin/sh -l
#SBATCH --time=00:45:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --array=1-4728

#An array job computed with the computation cluster. Bash variable $SLURM_ARRAY_TASK_ID contains the array job ID.

source /path/to/virtual/environment/bin/activate

INPUTFOLDER=/path/to/input/file/folder

#Find the "$SLURM_ARRAY_TASK_ID"th randomly picked genomic window index
WINDOW=$(sed ''"$SLURM_ARRAY_TASK_ID"'q;d' "$INPUTFOLDER"/random_indices_chr22.txt)
echo $WINDOW

SCRIPTFOLDER=/path/to/scripts
OUTPUTFOLDER=/path/to/output/folder
OUTPUTFILE=LuxUS_BFs_Hansen_data_chr22.txt
TIMEFILE=/path/to/output/folder/LuxUS_computation_times_chr22.txt

ALGORITHM=0
TEST_COV_IND=1
SIGMAB2=15
D_PLOTS=0

N_SAMPLES=1500

SIGMAR2_FILE="$OUTPUTFOLDER"/sample_means_sigmaR2_chr22.txt
SIGMAC2_FILE="$OUTPUTFOLDER"/sample_means_sigmaC2_chr22.txt
SIGMAE2_FILE="$OUTPUTFOLDER"/sample_means_sigmaE2_chr22.txt


python "$SCRIPTFOLDER"/run_luxus.py -d input_for_luxus_v3_"$WINDOW".txt -o $OUTPUT_FOLDER -i $INPUT_FOLDER -j $OUTPUT_FILE -c $TEST_COV_IND -b 15 -a $ALGORITHM -p $D_PLOTS -m $N_SAMPLES -w $SLURM_ARRAY_TASK_ID -x $SIGMAR2_FILE -y $SIGMAC2_FILE -z $SIGMAE2_FILE  -t $TIMEFILE
