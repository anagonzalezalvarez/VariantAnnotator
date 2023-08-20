#!/bin/bash --login

#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=6  
#SBATCH --mem=64G
#SBATCH --partition=batch 
#SBATCH -o out_ibex/%x-%j.out

conda activate ./env
python /ibex/scratch/gonzalac/variant_annotator/scripts/notebooks_scripts/11_Get_structure.py
 
