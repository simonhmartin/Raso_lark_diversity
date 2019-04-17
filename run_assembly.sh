#!/bin/sh

#SBATCH -p bigmem
#SBATCH -n 64
#SBATCH -N 1
#SBATCH --mem 512000
#SBATCH -t 96:00:00
#SBATCH -J assemble
#SBATCH -o assemble%j.out
#SBATCH -e assemble%j.err

RunAllPathsLG \
 PRE=/PATH/TO/PARENT/DIR\
 REFERENCE_NAME=skylark\
 DATA_SUBDIR=data\
 RUN=skylark_genome_2nd\
 HAPLOIDIFY=TRUE\
 OVERWRITE=True\
 THREADS=64\  ### The level of parallelization. For maximum performance, set this value to the number of processors available
 
