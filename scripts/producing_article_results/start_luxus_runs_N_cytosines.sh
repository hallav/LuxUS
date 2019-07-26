#This script calls run_calculate_BF_luxus_N_cytosines.sh, which then calls srun_calculate_BF_luxus_N_cytosines.sh

RUN_SEP=0 #Run the LuxUS analysis separately for each cytosine
RUN_LUXUS=1 #Run the LuxUS analysis for the whole window


for N_CYT in {1..20}
do 
    for REPS in 6
    do

        for READS in 12

        do

            for DIFF in 0 1
            do
                #Start computation cluster runs for the following script
                sbatch run_calculate_BF_luxus_N_cytosines.sh $DIFF $REPS $READS $RUN_V3_SEP $RUN_V3 $N_CYT

            done

        done
    done
done
    
