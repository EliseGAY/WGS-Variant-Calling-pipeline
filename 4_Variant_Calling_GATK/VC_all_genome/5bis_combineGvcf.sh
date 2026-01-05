#!/usr/bin/bash
#SBATCH -p std
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=gatk
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=200G
#SBATCH -V

# IMPORT MODULE
# module load singularity/4.2.2 / if gatk not installed

# INPUTS : 

# File with one scaffold or chromosome per lines
chr_list=$(cat chr.list)

# List of vcf file in the following format : --variant /PATH_TO_SAMPLE1/SAMPLE1_gatk.vcf.gz --variant  /PATH_TO_SAMPLE1/SAMPLE1_gatk.vcf.gz
VCFs_File="/path_to/VCF_list.txt"
Genome="PATH_TO/Genome.fasta" # index .fai and .dict have to be present in the same directory (see samtools faidx and gatk CreateSequenceDictionary functions)

# Run GATK
Gatk --java -jar CombineGVCFs \
-R ${Genome} \
--variant ${VCFs_File} \
# -L $CHR \ : if you want to do that by chromosome or on a subset
-O All.vcf.gz
