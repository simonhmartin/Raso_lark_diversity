#!/bin/sh

#SBATCH -p bigmem
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --mem 100000
#SBATCH -t 48:00:00
#SBATCH -J prepare
#SBATCH -o prepare_%j.out
#SBATCH -e prepare_%j.err


# NOTE: The option GENOME_SIZE is OPTIONAL. 
#       It is useful when combined with FRAG_COVERAGE and JUMP_COVERAGE 
#       to downsample data sets.
#       By itself it enables the computation of coverage in the data sets 
#       reported in the last table at the end of the preparation step. 

# NOTE: If your data is in BAM format you must specify the path to your 
#       picard tools bin directory with the option: 
#
#       PICARD_TOOLS_DIR=/your/picard/tools/bin

PrepareAllPathsInputs.pl\
 DATA_DIR=/PATH/TO/DATA\
 PLOIDY=2\
 IN_GROUPS_CSV=assembly_data/in_groups.csv\
 IN_LIBS_CSV=assembly_data/in_libs.csv\
 | tee prepare.out
 
