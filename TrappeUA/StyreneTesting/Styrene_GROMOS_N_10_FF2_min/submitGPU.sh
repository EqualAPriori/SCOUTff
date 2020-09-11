#!/bin/bash
#SBATCH -N 1 --partition=gpu --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:v100:1
#SBATCH --job-name=_N_3

cd $SLURM_SUBMIT_DIR

conda activate py3
srun --gres=gpu:1 python TIP4Pew.py
sleep 20

