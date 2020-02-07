#!/bin/sh -l

module load GCC/5.4.0-2.25
cd /path/to/scripts
source /path/to/virtual/environment/bin/activate

OUTPUTFOLDER=/path/to/output/folder
INPUTFOLDER=/path/to/input/folder

if (($#==6)); then

    #INPUT PARAMETERS

    DIFF=$1
    REPS=$2
    READS=$3
    RUN_V3_SEP=$4
    RUN_V3=$5
    N_CYT=$6

    if [[ "$4" == 1 ]]
    then

        python luxus_calculate_BF.py -d $SLURM_ARRAY_TASK_ID -w 1000 -l 38 -p 2 -c $N_CYT -r $REPS -q $READS -f "$OUTPUFOLDER"/N_cytosines_"$N_CYT" -i "$INPUTFOLDER"/N_cytosines_"$N_CYT" -b 15 -a $DIFF -x 0 -y 0 -z 0 -t 0 -s 1

    fi

    if [[ "$5" == 1 ]]
    then

        python luxus_calculate_BF.py -d $SLURM_ARRAY_TASK_ID -w 1000 -l 38 -p 2 -c $N_CYT -r $REPS -q $READS -f "$OUTPUFOLDER"/N_cytosines_"$N_CYT" -i "$INPUTFOLDER"/N_cytosines_"$N_CYT" -b 15 -a $DIFF -x 0 -y 0 -z 0 -t 1 -s 0

    fi

else

    echo 'Wrong amount of arguments'

fi

