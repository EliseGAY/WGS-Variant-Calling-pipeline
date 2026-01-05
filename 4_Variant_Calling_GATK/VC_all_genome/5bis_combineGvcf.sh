#!/usr/bin/bash
#SBATCH -p std
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --job-name=gatk
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH -V

# IMPORT MODULE
# module load singularity/4.2.2 / if gatk not installed

# INPUTS : 

# File with one scaffold or chromosome per lines
chr_list=$(cat chr.list)

# List of vcf file in the following format : --variant /PATH_TO_SAMPLE1/SAMPLE1_gatk.vcf.gz --variant  /PATH_TO_SAMPLE1/SAMPLE1_gatk.vcf.gz
VCFs_Path="/path_to_vcfs/"
Samples="sample1 sample2 sample3"
Genome="/scratch/lasojada/Fourmis/ref_seq/flye_polished_yahs_curation_round2.1.primary.curated.fasta"
DIRECTORY=$PWD

# create variant vcf list formated for the Combine function
VARIANTS=() # initialize the variant list
for i in $Samples; do
    file="$VCFs_Path/${i}.vcf.gz" # create full vcf path
    [[ -f "$file" ]] && VARIANTS+=( "--variant" "$file" )
done
echo "${VARIANTS[@]}" # print variant list 

mkdir All
cd All

# Run GATK
Gatk --java -jar CombineGVCFs \
-R ${Genome} \
${ListOfFiles} \
# -L $CHR \ : if you want to do that by chromosome
-O Cata_All.vcf.gz
