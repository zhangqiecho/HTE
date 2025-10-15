#!/usr/bin/env bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qi.zhang2@emory.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=15
#SBATCH --time=7-00:00:00
#SBATCH --job-name=sim1000
#SBATCH --mem=30GB
#SBATCH --partition=naimi
#SBATCH --array=1-15
#SBATCH --output=/home/qzha223/scenario/output/job.%J.out

config=/home/qzha223/scenario/data/config_1000.txt
c_number=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
sample_size=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
number_sims=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)


start_time="$(date -u +%s)"
module purge
module load R/4.2.2

Rscript --no-save --no-restore --verbose run_cluster_sim.R $c_number $sample_size $number_sims  > /home/qzha223/scenario/output/sim_c${c_number}_n${sample_size}_sim${number_sims}.Rout 2>&1

end_time="$(date -u +%s)"
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the number of confounder is ${c_number}, the sample size is ${sample_size}, and the number of simulations is ${number_sims}." 
elapsed="$(($end_time-$start_time))"
echo "Total of $elapsed seconds elapsed for this job!"