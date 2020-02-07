#!/bin/bash -l
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=6GB
#SBATCH --array=1-3000

#An array job computed with the computation cluster. Bash variable $SLURM_ARRAY_TASK_ID contains the array job ID.

source /path/to/virtual/environment/bin/activate

INPUTFOLDER=/path/to/input/data

WINDOW=$(sed ''"$SLURM_ARRAY_TASK_ID"'q;d' "$INPUTFOLDER"/HebestreitEtAl_DMR_saved_window_indices.txt)
echo "${WINDOW}"

SCRIPTFOLDER=/path/to/scripts

OUTPUTFOLDER=/path/to/output/folder
OUTPUTFILE=LuxUS_BFs_DMR.txt
TIMEFILE=/path/to/output/folder/DMR_LuxUS_computation_times.txt

ALGORITHM=0
TEST_COV_IND=1
SIGMAB2=15
D_PLOTS=1

N_SAMPLES=1500

SIGMAR2_FILE="$OUTPUTFOLDER"/DMR_estimated_sigmaR2.txt
SIGMAC2_FILE="$OUTPUTFOLDER"/DMR_estimated_sigmaC2.txt
SIGMAE2_FILE="$OUTPUTFOLDER"/DMR_estimated_sigmaE2.txt


python "$SCRIPTFOLDER"/run_luxus.py -d input_for_luxus_"${WINDOW}".txt -o $OUTPUT_FOLDER -i $INPUT_FOLDER -j $OUTPUT_FILE -c $TEST_COV_IND -b 15 -a $ALGORITHM -p $D_PLOTS -m $N_SAMPLES -w $SLURM_ARRAY_TASK_ID -x $SIGMAR2_FILE -y $SIGMAC2_FILE -z $SIGMAE2_FILE  -t $TIMEFILE

