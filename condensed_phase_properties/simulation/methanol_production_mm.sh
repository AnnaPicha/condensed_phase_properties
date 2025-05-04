#!/bin/bash
#SBATCH --gres=gpu
#SBATCH -p gpu
#SBATCH -p 4090
##SBATCH -p 40xx
##SBATCH --gres=gpu
##SBATCH -w n0034

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate openmm-ml-mace-nutmeg



# this block has to be edited for other simulation runs----------------------------

theory="mm" #mm/ani2x/mace-s/mace-m/mace-sl/nutemg-s
step=200000 #100 ps
system_name="moh323"
ensemble="NVT"  #NVT/NPT
int_method="NH" #""L"/"NH" choose "NH" for NVT runs

job_name="methanol_production_mm"   
last=36
# ---------------------------------------------------------------------------------



file_name=${theory}"_"${system_name}  # pass this to name .rst, .csv, .dcd files
psf="/site/raid7/anna/cp_prop/methanol/coor_moh/"${system_name}".psf"
crd="/site/raid7/anna/cp_prop/methanol/coor_moh/"${system_name}".crd"






if [ $1 -le $last ] ; then
    inp_file="production.py"
    python $inp_file -cnt $1 -p $psf -c $crd -n ${file_name} -nvt ${last} -m ${theory} -e ${ensemble} -int ${int_method} -s ${step} -sys ${system_name} 
fi

ret_val=$? # return value -> indicates wheter the command above (=running the python script) was succesful. 0 = succesful, 1 = error occured


if [ $ret_val -eq 0 -a $1 -lt $last ] ; then
    let next=$1+1
    sbatch -J ${job_name}_$next -o out/$next.out methanol_production_mm.sh $next 
elif [ ! $ret_val -eq 0 ] ; then
    echo "Error in run $1 .." >> out/error.log
fi
