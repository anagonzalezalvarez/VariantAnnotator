#!/bin/bash --login

#SBATCH --time 2:00:00
#SBATCH --cpus-per-task=4  
#SBATCH --mem-per-cpu=9G
#SBATCH --partition=batch
#SBATCH --job-name=jupyter
#SBATCH -o out_ibex/%x-%j.out

# use srun to launch Jupyter server in order to reserve a port
srun --resv-ports=1 bin/launch-jupyter.srun
