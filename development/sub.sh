#!/bin/bash

#SBATCH -J spatial_results
#SBATCH -o spatial_%j.txt
#SBATCH -e spatial_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=whitealj@iu.edu
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=24
#SBATCH --time=16:00:00
 
module load r/4.0.4

cd "/geode2/home/u100/whitealj/BigRed3/scripts/"
R CMD BATCH sub.R
