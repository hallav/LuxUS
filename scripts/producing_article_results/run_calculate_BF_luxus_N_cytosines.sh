#!/bin/sh -l 
#SBATCH --time=15:0:00 
#SBATCH --mem-per-cpu=15G 
#SBATCH --array=0-199

#An array job computed with the computation cluster. Bash variable $SLURM_ARRAY_TASK_ID contains the array job ID.

cd /path/to/scripts

if (($#==7)); then

    #INPUT PARAMETERS

    DIFF=$1
    REPS=$2
    READS=$3
    RUN_V3_SEP=$4
    RUN_V3=$5
    RUN_V1=$6
    N_CYT=$7

    #Call the following bash script
    srun  srun_calculate_BF_luxus_N_cytosines.sh $DIFF $REPS $READS $RUN_V3_SEP $RUN_V3 $RUN_V1 $N_CYT
else

    echo 'Wrong amount of arguments'

fi
