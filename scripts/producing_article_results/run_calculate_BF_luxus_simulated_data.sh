#!/bin/sh -l 
#SBATCH --time=10:0:00 
#SBATCH --mem-per-cpu=15G 
#SBATCH --job-name=LuxUS
#SBATCH --array=0-200

cd /path/to/scripts/folder


if (($#==6)); then

    #INPUT PARAMETERS

    DIFF=$1
    REPS=$2
    READS=$3
    MU=$4
    RUN_SEP=$5
    RUN_LUXUS=$6

    srun  srun_calculate_BF_simulated_data.sh $DIFF $REPS $READS $MU $RUN_SEP $RUN_LUXUS
else

    echo 'Wrong amount of arguments'

fi
