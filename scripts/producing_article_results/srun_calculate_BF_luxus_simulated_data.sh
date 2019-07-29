#!/bin/sh -l

module load GCC/5.4.0-2.25
cd /path/to/scripts/folder
source /path/to/virtualenvironment/bin/activate

if (($#==6)); then

    #INPUT PARAMETERS

    DIFF=$1
    REPS=$2
    READS=$3
    MU=$4
    RUN_SEP=$5
    RUN_LUXUS=$6

    if [[ "$5" == 1 ]]
    then

        python luxus_calculate_BF.py -d $SLURM_ARRAY_TASK_ID -w 1000 -l 38 -p 2 -c 10 -r $REPS -q $READS -f /path/to/results/folder/"$MU" -i /path/to/generated/data/"$MU" -b 15 -a $DIFF -x 0 -y 0 -z 0 -t 0 -s 1

    fi

    if [[ "$6" == 1 ]]
    then

        python luxus_calculate_BF.py -d $SLURM_ARRAY_TASK_ID -w 1000 -l 38 -p 2 -c 10 -r $REPS -q $READS -f /path/to/results/folder/"$MU" -i /path/to/generated/data/"$MU" -b 15 -a $DIFF -x 0 -y 0 -z 0 -t 1 -s 0

    fi

else

    echo 'Wrong amount of arguments'

fi

