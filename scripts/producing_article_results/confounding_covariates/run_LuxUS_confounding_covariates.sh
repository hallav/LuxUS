#!/bin/sh -l 
#SBATCH --time=05:00:00 
#SBATCH --mem-per-cpu=5G
#SBATCH --job-name=LuxUS_confounding_covariates
#SBATCH --mail-type=ALL
#SBATCH --mail-user=viivi.halla-aho@aalto.fi
#SBATCH --array=0-949

cd /path/to/scripts

if (($#==5)); then

    #INPUT PARAMETERS

    DIFF=$1
    REPS=$2
    READS=$3
    RUN_LUXUS_SEP=$4
    RUN_LUXUS=$5

    srun  srun_LuxUS_confounding_covs.sh $DIFF $REPS $READS $MU $RUN_LUXUS_SEP $RUN_LUXUS
else

    echo 'Wrong amount of arguments'

fi
