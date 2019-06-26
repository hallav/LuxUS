

RUN_V3_SEP=0
RUN_V3=1
RUN_V1=0

declare -a muB=("mu_1_4_minus1")

for N_CYT in {1..20}
do 
    for REPS in 6
    do

        for READS in 12

        do

            for DIFF in 0 1
            do
                #Start computation cluster runs for the following script
                sbatch run_calculate_BF_luxus_N_cytosines.sh $DIFF $REPS $READS $RUN_V3_SEP $RUN_V3 $RUN_V1 $N_CYT

            done

        done
    done
done
    
