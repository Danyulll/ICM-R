#!/bin/bash
#SBATCH --account=def-dsteph          # replace with your supervisor's account
#SBATCH --nodes=1                       # must be 1
#SBATCH --cpus-per-task=4               # number of processes/threads
#SBATCH --mem-per-cpu=16G             # memory per CPU (MB)
#SBATCH --time=0-00:15                  # time (DD-HH:MM)

module load r/4.5.0


# Run your script (change the filename to your .R file)
R CMD BATCH --no-save --no-restore model.R
