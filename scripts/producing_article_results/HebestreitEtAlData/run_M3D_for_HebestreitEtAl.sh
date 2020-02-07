#!/bin/sh -l

cd $SCRIPTFOLDER

module load r/3.5.0-python-2.7.14

#The input arguments should be 
#1: The input data file name identifier (with path!) for DMR
#2: The input data file name identifier (with path!) for nonDMR
#3: The number of simulated data files to be used as input.
#4: Output filename identifier (with path!)
#5: Maximum number of cytosines per window


N_DATASETS=3000
NCYT=20

INPUTDATAFOLDER=/path/to/input/data/folder
OUTPUTFOLDER=/path/to/output/folder

OUTPUTIDENTIFIER="$OUTPUTFOLDER"/M3D_results_KleinHebestreit_addedRequirements2_realChr_eiPurkkaa
INPUTFILEDMR="$INPUTDATAFOLDER"/TEMP_WINDOW_proportion_table_HebestreitEtAl_DMR_wHeader_
INPUTFILENONDMR="$INPUTDATAFOLDER"/TEMP_WINDOW_proportion_table_HebestreitEtAl_nonDMR_wHeader_

srun Rscript --no-save run_M3D_for_HebestreitEtAl.R $INPUTFILEDMR $INPUTFILENONDMR $N_DATASETS $OUTPUTIDENTIFIER $NCYT
