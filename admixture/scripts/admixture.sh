#!/bin/bash
#SBATCH -J admxixture
#SBATCH -o slurm_logs/admkxture_%A.out
#SBATCH -e slurm_logs/admixture_%A.err
#SBATCH -p shared,serial_requeue
#sBATCH -N 1
#SBATCH -n 6                
#SBATCH -t 6:00:00                   
#SBATCH --mem=12000

module load python
source activate admixture

K=$1
mkdir -p K${K}

for i in `seq 1 100`

do admixture -j6 -s time --cv=10 chromfix_nosingletons_Phase_3_10k_FINAL_NEUTRAL.bed $K > chromfix_nosingletons_Phase_3_10k_FINAL_NEUTRAL.bed_k${K}_rep${i}.log 2>&1
mv chromfix_nosingletons_Phase_3_10k_FINAL_NEUTRAL.bed_k${K}_rep${i}.log logs/chromfix_nosingletons_Phase_3_10k_FINAL_NEUTRAL.bed_k${K}_rep${i}.log
mv chromfix_nosingletons_Phase_3_10k_FINAL_NEUTRAL.${K}.Q K${K}/chromfix_nosingletons_Phase_3_10k_FINAL_NEUTRAL.${K}_r${i}.Q
mv chromfix_nosingletons_Phase_3_10k_FINAL_NEUTRAL.${K}.P K${K}/chromfix_nosingletons_Phase_3_10k_FINAL_NEUTRAL.${K}_r${i}.P
done 
