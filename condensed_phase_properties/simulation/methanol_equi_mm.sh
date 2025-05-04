#!/bin/bash
#SBATCH -p ADA
##SBATCH -p 4090
##SBATCH -p 40xx
##SBATCH --gres=gpu
##SBATCH -w n0034

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate openmm-ml-mace-nutmeg



# this block has to be edited for other simulation runs----------------------------

theory="mm" #mm/ani2x/mace-s/mace-m/mace-sl/nutemg-s
system_name="moh323" #"moh323"/"moh73"
ensemble="NPT"  #NVT/NPT
step=20000000 #1.1 ns

job_name="methanol_equi_mm"   


python equi.py ${theory} ${system_name} ${ensemble} ${step}







