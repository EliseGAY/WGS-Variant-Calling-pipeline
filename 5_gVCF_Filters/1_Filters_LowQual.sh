#!/usr/bin/bash
#SBATCH -A nuclear_genic_sequences_reconstruction
#SBATCH -o Filter.o
#SBATCH -e Filter.e
#SBATCH --job-name=Filter
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=01-00:00:00
#SBATCH --partition=long
#SBATCH --mem-per-cpu=100G

Input="/your/path/XXX.TAG.vcf.gz"

#------------------------------------------------------------#
# Filtering: Remove low-quality positions and unwanted variants
# - Remove LowQual, MQFILTER, Repeat, and Indels
#------------------------------------------------------------#

bcftools view ${Input} -e 'FILTER~"MQFILTER|LowQual|Repeat" || TYPE="indel"' --threads 8 -Oz -o VCF.TAG.Flowqual_Noindels_Norepeat.vcf.gz

# Example of filters on SNP (to keep snp only), missing data, and MAF and. To filter the DP and missing data along with custom vizualisation see the DP_Na_Filter R script.

##TO DO
