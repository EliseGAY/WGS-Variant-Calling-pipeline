#!/bin/sh
#SBATCH --clusters=mesopsl1
#SBATCH --account=gay
#SBATCH --partition=def
#SBATCH --qos=mesopsl1_def_long
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --job-name=Index_BCFtools
#SBATCH --time=20:00:00
#SBATCH -o index_BCFtools.o
#SBATCH -e index_BCFtools.e

# IMPORT MODULE
module load gcc/9.2.0
module load samtools/1.10
module load gatk/4.2.0.0

# Run GATK indexation
bcftools index 57st_bcftools.vcf.gz