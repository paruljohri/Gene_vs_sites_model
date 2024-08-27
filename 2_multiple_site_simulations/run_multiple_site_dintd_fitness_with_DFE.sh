#!/bin/bash
#SBATCH --mail-user=pjohri@unc.edu
#SBATCH --mail-type=ALL
#SBATCH -p general
#SBATCH -n 1 #number of tasks
#SBATCH --time=0-20:00
#SBATCH --mem=20m
#SBATCH -a 1-1000%500
#SBATCH -o /nas/longleaf/home/pjohri/LOGFILES/fitness_%A_rep%a.out
#SBATCH -e /nas/longleaf/home/pjohri/LOGFILES/fitness_%A_rep%a.err

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
######################################
#To be run on LONGLEAF!
#9 simulations per submitted job!
#######################################


######NOTES########
#theta = Nu = 0.01
#20 hours for 1000 sites and a DFE
#10 hours for 1000 sites
#5 hours for 100 sites
#The only difference between _LD and _LD_DFE is that in _LD_DFE we simulate a DFE, not jsut constant "s"
###################

declare -i repID=0+$SLURM_ARRAY_TASK_ID

module load slim

cd /nas/longleaf/home/pjohri/FitnessNote/programs

echo "starting simulation " $repID

num_sites=1000
theta="theta0_005"
fitness_model="additive_site" #multiplicative/additive_site/additive_gene
folder="/work/users/p/j/pjohri/FitnessNote/simulations/thousand_site/${fitness_model}/${theta}"

########################
#DFE parameters:
#gamma=2, beta=0.3: "0.610, 0.348, 0.042, 1.47e-08"
#gamma=20, beta=0.3: "0.315, 0.295, 0.348, 0.042"
#gamma=100, beta=0.3: "0.195, 0.192, 0.340, 0.273"
########################

##gamma=-2; h=0.5
#if [ ! -f "${folder}/gamma2_mean/h_0_5/output${repID}.txt" ]
#then
#	slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_mean_gamma=-2.0 -d d_beta=0.3 -d d_dom=0.5 -d "d_folder='${folder}/gamma2_mean/h_0_5'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD_DFE.slim
#fi

##gamma=-2; h=0.2
#if [ ! -f "${folder}/gamma2_mean/h_0_2/output${repID}.txt" ]
#then
#	slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_mean_gamma=-2.0 -d d_beta=0.3 -d d_dom=0.2 -d "d_folder='${folder}/gamma2_mean/h_0_2'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD_DFE.slim
#fi

##gamma=-2; h=0.0
#if [ ! -f "${folder}/gamma2_mean/h_0_0/output${repID}.txt" ]
#then
#	slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_mean_gamma=-2.0 -d d_beta=0.3 -d d_dom=0.0 -d "d_folder='${folder}/gamma2_mean/h_0_0'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD_DFE.slim
#fi

##gamma=-20; h=0.5
#if [ ! -f "${folder}/gamma20_mean/h_0_5/output${repID}.txt" ]
#then
#	slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_mean_gamma=-20.0 -d d_beta=0.3 -d d_dom=0.5 -d "d_folder='${folder}/gamma20_mean/h_0_5'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD_DFE.slim
#fi

##gamma=-20; h=0.2
#if [ ! -f "${folder}/gamma20_mean/h_0_2/output${repID}.txt" ]
#then
#	slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_mean_gamma=-20.0 -d d_beta=0.3 -d d_dom=0.2 -d "d_folder='${folder}/gamma20_mean/h_0_2'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD_DFE.slim
#fi

##gamma=-20; h=0.0
#if [ ! -f "${folder}/gamma20_mean/h_0_0/output${repID}.txt" ]
#then
#	slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_mean_gamma=-20.0 -d d_beta=0.3 -d d_dom=0.0 -d "d_folder='${folder}/gamma20_mean/h_0_0'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD_DFE.slim
#fi

##gamma=-100; h=0.5
#if [ ! -f "${folder}/gamma100_mean/h_0_5/output${repID}.txt" ]
#then
#        slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_mean_gamma=-100.0 -d d_beta=0.3 -d d_dom=0.5 -d "d_folder='${folder}/gamma100_mean/h_0_5'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD_DFE.slim
#fi

##gamma=-100; h=0.2
#if [ ! -f "${folder}/gamma100_mean/h_0_2/output${repID}.txt" ]
#then
#        slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_mean_gamma=-100.0 -d d_beta=0.3 -d d_dom=0.2 -d "d_folder='${folder}/gamma100_mean/h_0_2'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD_DFE.slim
#fi

##gamma=-100; h=0.0
#if [ ! -f "${folder}/gamma100_mean/h_0_0/output${repID}.txt" ]
#then
#        slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_mean_gamma=-100.0 -d d_beta=0.3 -d d_dom=0.0 -d "d_folder='${folder}/gamma100_mean/h_0_0'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD_DFE.slim
#fi

##gamma=-1000; h=0.5
if [ ! -f "${folder}/gamma1000_mean/h_0_5/output${repID}.txt" ]
then
        slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_mean_gamma=-1000.0 -d d_beta=0.3 -d d_dom=0.5 -d "d_folder='${folder}/gamma1000_mean/h_0_5'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD_DFE.slim
fi

##gamma=-1000; h=0.2
if [ ! -f "${folder}/gamma1000_mean/h_0_2/output${repID}.txt" ]
then
        slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_mean_gamma=-1000.0 -d d_beta=0.3 -d d_dom=0.2 -d "d_folder='${folder}/gamma1000_mean/h_0_2'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD_DFE.slim
fi

##gamma=-1000; h=0.0
if [ ! -f "${folder}/gamma1000_mean/h_0_0/output${repID}.txt" ]
then
        slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_mean_gamma=-1000.0 -d d_beta=0.3 -d d_dom=0.0 -d "d_folder='${folder}/gamma1000_mean/h_0_0'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD_DFE.slim
fi

echo "Finished simulation " $repID





