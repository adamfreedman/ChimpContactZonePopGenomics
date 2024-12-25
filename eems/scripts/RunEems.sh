#!/bin/bash
#SBATCH -J eems
#SBATCH -o eems_%A.out
#SBATCH -e eems_%A.err
#SBATCH -p shared,serial_requeue
#SBATCH -n 1                
#SBATCH -t 06:00:00                   
#SBATCH --mem=12000 

/PATH/TO/eems/runeems_snps/src/runeems_snps --params $1


