#!/bin/sh -l
#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu=1GB
#SBATCH --job-name=Prepare_data_for_LuxUS

source /path/to/virtual/environment/bin/activate

INPUT_FILE=wgbs_methylationlevels_chr21.txt
INPUT_FOLDER=/path/to/Hansen/data
OUTPUT_FOLDER=/path/to/output/folder/chr21
DESIGN_MATRIX=/path/to/Hansen/data/design_matrix.txt
MEANCOVFILE=/path/to/output/folder/chr21/mean_coverage_chr21.txt
CYTNFILE=/path/to/output/folder/chr21/number_of_cytosines_in_window_chr21.txt

#The BSEFF values have been estimated as a mean of the two flowcell rates from Hansen et al. 2011, supplement table 14. Rounded to 4 decimals.
BSEFF='[0.9978,0.9977,0.9979,0.9975,0.9979,0.9979]'
BSBEFF='[0,0,0,0,0,0]'
SEQERR='[0,0,0,0,0,0]'
TEST_COV_IND=1
WINDOW_WIDTH=2000
REQUIRED_DIFF=0.1

python /path/to/script/folder/prepare_data_for_luxus_Hansen.py -i "$INPUT_FOLDER"/"$INPUT_FILE" -d $DESIGN_MATRIX -o $OUTPUT_FOLDER -r 6  -t $TEST_COV_IND -w $WINDOW_WIDTH -u $REQUIRED_DIFF -y $MEANCOVFILE -z $CYTNFILE
