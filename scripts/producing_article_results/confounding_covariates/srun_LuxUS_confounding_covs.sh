#!/bin/sh -l

module load GCC/5.4.0-2.25
cd /path/to/scripts
source /path/to/virtual/environment/bin/activate

if (($#==6)); then

    #INPUT PARAMETERS

    DIFF=$1
    REPS=$2
    READS=$3
    RUN_LUXUS_SEP=$4
    RUN_LUXUS=$5

    RESULTFOLDER=/path/to/results
    INPUTFOLDER=/path/to/generated/data

    if [[ "$4" == 1 ]]
    then

        python run_LuxUS_confounding_covs.py -d $SLURM_ARRAY_TASK_ID -w 1000 -l 38 -p 2 -c 10 -r $REPS -q $READS -f $RESULTFOLDER -i $INPUTFOLDER -b 15 -a $DIFF -x 0 -y 0 -z 0 -t 0 -s 1

    fi

    if [[ "$5" == 1 ]]
    then

        python run_LuxUS_confounding_covs.py -d $SLURM_ARRAY_TASK_ID -w 1000 -l 38 -p 2 -c 10 -r $REPS -q $READS -f $RESULTFOLDER -i $INPUTFOLDER -b 15 -a $DIFF -x 0 -y 0 -z 0 -t 1 -s 0

    fi


else

    echo 'Wrong amount of arguments'

fi

