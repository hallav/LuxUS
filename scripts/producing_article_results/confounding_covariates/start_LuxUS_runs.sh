#This script starts SLURM jobs for script run_LuxUS_confounding_covariates.sh which calls srun_LuxUS_confounding_covs.sh

RUN_V3_SEP=0
RUN_V3=1


for REPS in 6 12
do

    for READS in 6 12 24
    do

        for DIFF in 0
        do

            sbatch run_LuxUS_confounding_covariates.sh $DIFF $REPS $READS $RUN_V3_SEP $RUN_V3

        done

    done
done
    
