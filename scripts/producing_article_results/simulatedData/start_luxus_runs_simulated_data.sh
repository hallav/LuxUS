#This script calls run_calculate_BF_luxus_simulated_data.sh, which in turn calls srun_calculate_BF_luxus_simulated_data.sh

RUN_SEP=1 #Run the LuxUS analysis separately for each cytosine
RUN_LUXUS=0 #Run the LuxUS analysis for the whole window

declare -a muB=("mu_minus1_4_1" "mu_minus1_4_2_3" "mu_minus1_4_6" "mu_1_4_minus1" "mu_1_4_minus2_3" "mu_1_4_minus6") #Go through all the simulation settings for beta mean, each stored in their own folder

for MU in "${muB[@]}"
do 
    for REPS in 6 12 24
    do

        for READS in 6 12 24
        do

            for DIFF in 0 1
            do

                sbatch run_calculate_BF_luxus_simulated_data.sh $DIFF $REPS $READS $MU $RUN_SEP $RUN_LUXUS

            done

        done
    done
done
    
