#!/bin/bash -l
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --job-name=run_LuxUS_preanalysis_NSCLC
#SBATCH --mail-type=ALL

cd /path/to/NSCLC/meth_files

source /path/to/virtual/environment/bin/activate

while read i
do
  echo $i

  INPUTFOLDER=/path/to/NSCLC/meth_files
  PROP_TABLE="$INPUTFOLDER"/NSCLC_proportion_table_"$i"_luxusformat.txt
  OUTPUTFOLDER=/path/to/NSCLC/preanalysis/"$i"
  DESIGNMATRIX=NSCLC_design_matrix_noCol1.txt
  REPS=8
  N_REQUIRED_SAMPLES=1
  TESTCOV=1
  PVALREQ=0.1
  REQUIRED_COVERAGE=5
  MEANCOVFILE="$OUTPUTFOLDER"/mean_coverage_in_window_"$i".txt
  CYTNFILE="$OUTPUTFOLDER"/number_of_cytosines_in_window_"$i".txt
  WINDOW_WIDTH=2000


  ALPHA_E=2
  ALPHA_R=2
  ALPHA_C=2

  BETA_E=2
  BETA_C=2
  BETA_R=2


  python /path/to/scripts/prepare_data_for_luxus.py -i $PROP_TABLE -o $OUTPUTFOLDER -d $DESIGNMATRIX -r $REPS -s $N_REQUIRED_SAMPLES -t $TESTCOV -u $PVALREQ -y $MEANCOVFILE -z $CYTNFILE -v $REQUIRED_COVERAGE -w $WINDOW_WIDTH -m $ALPHA_E -n $BETA_E -k $ALPHA_R -l $BETA_R -g $ALPHA_C -f $BETA_C

done <chromosome_list_unique.txt
