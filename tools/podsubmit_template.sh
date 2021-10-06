#!/bin/bash
# ask for 16 cores on two nodes
#SBATCH --nodes=1 --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=167:00:00
#SBATCH --job-name=__label__
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=__USERNAME__@ucsb.edu

#module load intel/18

/bin/hostname
echo CUDA_VISIBLE:  $CUDA_VISIBLE_DEVICES
echo SLURM_JOB_GPU: $SLURM_JOB_GPUS
echo " "
srun --gres=gpu:1 /usr/bin/nvidia-smi


cd $SLURM_SUBMIT_DIR
#srun --gres=gpu:1 python ./run.py -ff /home/kshen/mylib/SCOUTff/tools/../TrappeUA/XML/TrappeUA_Butadiene_Gromos.xml > z.log 

