#!/bin/bash -l
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=200M
#SBATCH --array=1-3000

#An array job computed with the computation cluster. The bash variable $SLURM_ARRAY_TASK_ID contains the array run ID.

module load GSL/2.4-goolf-triton-2017a

SCRIPTFOLDER=/path/to/scripts
INPUTFOLDER=/path/to/input/data
RESULTFOLDER=/path/to/DMR/results/folder


/scratch/cs/csb/users/hallav1/software/methpipe-3.4.3/bin/radmeth regression -factor IsCase "$INPUTFOLDER"/design_simple_radmeth.txt "$INPUTFOLDER"/proportion_table_HebestreitEtAl_DMR_wHeader_"$SLURM_ARRAY_TASK_ID".txt > "$RESULTFOLDER"/radmeth_cpgs_HebestreitEtAl_"$SLURM_ARRAY_TASK_ID".bed

