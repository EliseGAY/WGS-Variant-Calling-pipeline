#!/bin/bash/

#==================#
# Load Directories
#==================#

interval_list="/travail/egay/Whole_Genome_analysis_GWS/GATK/GATK_23samples/VCF_Final/VCF_Filters/interval.list"

#=======================#
# loop on chr
#=======================#
while read chr;
 do

#==================#
# Create variables
#==================#

Input="/travail/egay/Whole_Genome_analysis_GWS/GATK/GATK_23samples/VCF_Final/${chr}_GATK_TAG.vcf.gz"
Basename="WGS_23_${chr}"

#========================================#
# Write on script for each chr.vcf files
#========================================#
        cat > ${chr}_Filters_pos.sh << EOF
#!/bin/sh

# SLURM job configurations
#SBATCH --clusters=mesopsl1              # Cluster to use
#SBATCH --account=gay                    # Account for billing
#SBATCH --partition=def                  # Partition to run the job on
#SBATCH --qos=mesopsl1_def_long          # Quality of service
#SBATCH --nodes=1                        # Number of nodes to use
#SBATCH --ntasks-per-node=10             # Number of tasks per node
#SBATCH --job-name=${chr}_pos            # Job name (using chromosome as part of the name)
#SBATCH --time=48:00:00                  # Maximum runtime (48 hours)
#SBATCH -o ${chr}_filters.o             # Standard output file
#SBATCH -e ${chr}_filters.e             # Standard error file

# Load necessary modules
module load gcc/9.2.0                   # Load GCC compiler
module load samtools/1.10               # Load Samtools for working with BAM/VCF files
module load gatk/4.2.0.0                # Load GATK for variant calling

#------------------------------------------------------------#
# Filtering: Remove low-quality positions and unwanted variants
# - Remove LowQual, MQFILTER, Repeat, and Indels
# - Keep only SNPs in the final output
#------------------------------------------------------------#

# Step 1: Extract the header from the input VCF and store it in a separate file
bcftools view ${Input} | grep "#" >> ${Basename}_GATK.header

# Step 2: Filter out low-quality positions, indels, and repeats:
# - Remove positions marked with MQFILTER, LowQual, or Repeat flags
# - Exclude positions with "*" in the REF/ALT columns and clean up the file
bcftools view -H ${Input} | grep -v "MQFILTER\|LowQual\|Repeat" |  awk '!(\$5 ~ "^*")' | awk '!(\$5 ~ "^[ACGT*].")' | awk '!(\$4 ~ "^[ACGT*].")' >> ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat_Temp.vcf

# Step 3: Re-add the header to the filtered VCF data
cat ${Basename}_GATK.header ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat_Temp.vcf >> ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat.vcf

# Step 4: Compress the filtered VCF file using bgzip
bcftools view ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat.vcf -o ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat.vcf.gz -Oz

# Step 5: Extract positions (no missing data) from the VCF
# - Filter the VCF to remove any sites with missing genotypes (max-missing 1)
# - Save the position (chromosome and position) to a separate file
bcftools view ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat.vcf.gz | vcftools --vcf '-' --max-missing 1 --stdout --recode | grep -v "#" | cut -f1,2 >> ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat_NoNa.positions

# Step 6: Generate a VCF file with no missing data
# - Remove positions with missing genotypes and save the cleaned VCF as a gzipped file
bcftools view ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat.vcf.gz | vcftools --vcf '-' --max-missing 1 --stdout --recode | bcftools view -o ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat_NoNa.vcf.gz -Oz

# Step 7: Clean up temporary files to save space
rm ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat_Temp.vcf ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat.vcf

# Step 8: Extract only SNPs (Single Nucleotide Polymorphisms) from the VCF
# - Remove indels and retain only SNPs in the final VCF file
bcftools view -v snps ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat.vcf.gz -o ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat_SNP.vcf.gz -Oz

# Step 9: Index the VCF files for efficient querying
bcftools index ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat.vcf.gz
bcftools index ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat_SNP.vcf.gz
bcftools index ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat_NoNa.vcf.gz

# Submit the SLURM job for the next chromosome or interval in the list
EOF
sbatch ${chr}_Filters_pos.sh
done < ${interval_list}
