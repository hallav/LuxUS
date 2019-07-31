#!/bin/sh -l
#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu=1GB
#SBATCH --job-name=Prepare_data_for_LuxUS


source /path/to/virtual/environment/bin/activate

INPUT_FILE=proportion_table_HebestreitEtAl_DMR_only_wLocation.txt
INPUT_FOLDER=/path/to/HebestreitEtAl/data
OUTPUT_FOLDER=/path/to/output/folder
DESIGN_MATRIX=/path/to/HebestreotEtAl/data/design_matrix.txt
MEANCOVFILE=/path/to/output/folder/mean_coverage_DMRs.txt
CYTNFILE=/path/to/output/folder/number_of_cytosines_in_window_DMRs.txt

#The experimental parameters have been set to the "perfect" values as the correct values are not known
BSEFF='[1,1,1,1,1,1,1,1,1,1,1,1]'
BSBEFF='[0,0,0,0,0,0,0,0,0,0,0,0]'
SEQERR='[0,0,0,0,0,0,0,0,0,0,0,0]'
TEST_COV_IND=1
WINDOW_WIDTH=2000
REQUIRED_DIFF=1
REQUIRED_COVERAGE=1
N_REQUIRED_SAMPLES=1

python /path/to/script/folder/prepare_data_for_luxus_HebestreitEtAl_saveRandomWindows.py -i "$INPUT_FOLDER"/"$INPUT_FILE" -d $DESIGN_MATRIX -o $OUTPUT_FOLDER -r 12  -t $TEST_COV_IND -w $WINDOW_WIDTH -u $REQUIRED_DIFF -y $MEANCOVFILE -z $CYTNFILE -v $REQUIRED_COVERAGE -s $N_REQUIRED_SAMPLES -x 3000
