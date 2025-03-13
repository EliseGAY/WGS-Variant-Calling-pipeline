## 📂 Inputs

1. **Input VCF file**: The VCF file from the previous step in the pipeline (e.g., from GATK's HaplotypeCaller).
2. **Interval list**: A list of intervals (chromosomes or scaffolds) for which to perform variant calling.
3. **Basename**: The base name used for output files (typically the sample name).
4. **Chromosome (chr)**: The current chromosome or scaffold being processed (set dynamically during job execution).

---

## 📝 Script Overview

The script performs the following steps:

1. **Header Extraction**:
   - Extracts the header from the input VCF file and saves it in a separate file.

2. **Variant Filtering**:
   - Removes low-quality variants (using `MQFILTER`, `LowQual`, `Repeat` flags).
   - Discards indels and other unwanted variants.
   - Uses `awk` commands to remove specific unwanted patterns, such as `"*"` or certain nucleotide patterns.

3. **VCF Compression**:
   - Compresses the filtered VCF file using `bcftools view` and `bgzip` to generate `.vcf.gz` files.

4. **No Missing Data**:
   - Filters out variants with missing genotypes using `vcftools` (removes positions with missing data).
   - Generates a position list and a VCF file without missing data.

5. **SNP Extraction**:
   - Extracts only SNPs (Single Nucleotide Polymorphisms) from the VCF, removing indels.

6. **VCF Indexing**:
   - Indexes the resulting VCF files for efficient access using `bcftools index`.

7. **Job Submission**:
   - Submits the next SLURM job to process the next interval/chromosome in the list using the `sbatch` command.

---

## 🖥️ How to Run

### Step 1: Prepare Input Files
Ensure that you have the following files prepared before running the script:

- A list of samples (VCF files) generated from previous steps.
- An interval list containing the chromosomes/scaffolds to be processed.
- Set the `Basename` and `chr` (chromosome) variables in the script as appropriate.

### Step 2: Run the Script

To execute the script on the cluster, use the following command:

`bash
sbatch filters_LowQual.sh
`
#======================#
# 02/2022
# Elise GAY
# Run Filters on gVCF
# please inform the authors before sharing
#======================#

# Aim : 
#------#
Run filters on gVCF with BCFTools

# Input :
#----------#
gVCF tagged

# Methods :
#----------#
run one job by chromosome : get each gVCF by loopoing on chromosome list
Get several vcf by filtering fo quality, SNP etc.. 

# output :
#----------#
# VCF with all clean position : ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat.vcf.gz
# VCF with SNP clean position : ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat_SNP.vcf.gz
# VCF with SNP clean position with no NA : ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat_NoNa.vcf.gz
