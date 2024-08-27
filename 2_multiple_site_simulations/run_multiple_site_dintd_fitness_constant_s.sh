#!/bin/bash
#SBATCH --mail-user=pjohri@unc.edu
#SBATCH --mail-type=ALL
#SBATCH -p general
#SBATCH -n 1 #number of tasks
#SBATCH --time=0-20:00
#SBATCH --mem=20m
#SBATCH -a 1-1000%1000
#SBATCH -o /nas/longleaf/home/pjohri/LOGFILES/fitness_%A_rep%a.out
#SBATCH -e /nas/longleaf/home/pjohri/LOGFILES/fitness_%A_rep%a.err

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
######################################
#To be run on LONGLEAF!
#6 simulations per submitted job!
#######################################


######NOTES########
#theta = Nu = 0.01
#10 hours for 1000 sites
#5 hours for 100 sites
#The only difference between v1 and v2 is that in v2 we also get an output of .ms files
###################

declare -i repID=0+$SLURM_ARRAY_TASK_ID

module load slim

cd /nas/longleaf/home/pjohri/FitnessNote/programs

echo "starting simulation " $repID

num_sites=1000
theta="theta0_01"
fitness_model="additive_site" #multiplicative/additive_site/additive_gene
folder="/work/users/p/j/pjohri/FitnessNote/simulations/thousand_site/${fitness_model}/${theta}/noEpistasis"

##gamma=-2; h=0.5
#if [ ! -f "${folder}/gamma2/h_0_5/output${repID}.txt" ]
#then
#	slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_gamma=-2.0 -d d_dom=0.5 -d "d_folder='${folder}/gamma2/h_0_5'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD.slim
#fi

##gamma=-2; h=0.2
#if [ ! -f "${folder}/gamma2/h_0_2/output${repID}.txt" ]
#then
#	slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_gamma=-2.0 -d d_dom=0.2 -d "d_folder='${folder}/gamma2/h_0_2'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD.slim
#fi

##gamma=-2; h=0.0
#if [ ! -f "${folder}/gamma2/h_0_0/output${repID}.txt" ]
#then
#	slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_gamma=-2.0 -d d_dom=0.0 -d "d_folder='${folder}/gamma2/h_0_0'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD.slim
#fi

##gamma=-20; h=0.5
#if [ ! -f "${folder}/gamma20/h_0_5/output${repID}.txt" ]
#then
#	slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_gamma=-20.0 -d d_dom=0.5 -d "d_folder='${folder}/gamma20/h_0_5'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD.slim
#fi

##gamma=-20; h=0.2
#if [ ! -f "${folder}/gamma20/h_0_2/output${repID}.txt" ]
#then
#	slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_gamma=-20.0 -d d_dom=0.2 -d "d_folder='${folder}/gamma20/h_0_2'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD.slim
#fi

##gamma=-20; h=0.0
#if [ ! -f "${folder}/gamma20/h_0_0/output${repID}.txt" ]
#then
#	slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_gamma=-20.0 -d d_dom=0.0 -d "d_folder='${folder}/gamma20/h_0_0'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD.slim
#fi

##gamma=-100; h=0.5
if [ ! -f "${folder}/gamma100/h_0_5/output${repID}.txt" ]
then
        slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_gamma=-100.0 -d d_dom=0.5 -d "d_folder='${folder}/gamma100/h_0_5'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD.slim
fi

##gamma=-100; h=0.2
if [ ! -f "${folder}/gamma100/h_0_2/output${repID}.txt" ]
then
        slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_gamma=-100.0 -d d_dom=0.2 -d "d_folder='${folder}/gamma100/h_0_2'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD.slim
fi

##gamma=-100; h=0.0
if [ ! -f "${folder}/gamma100/h_0_0/output${repID}.txt" ]
then
        slim -d d_seed=${repID} -d d_num_sites=${num_sites} -d d_mut_rate=1e-5 -d d_rec_rate=1e-5 -d "d_fitness_model='${fitness_model}'" -d d_gamma=-100.0 -d d_dom=0.0 -d "d_folder='${folder}/gamma100/h_0_0'" -d "d_repID='${repID}'" multiple_site_dintd_fitness_LD.slim
fi

echo "Finished simulation " $repID





