## 📂 Inputs

1. **Input VCF file**: The VCF file from the previous step in the pipeline (e.g., from GATK's HaplotypeCaller).
2. **Interval list**: A list of intervals (chromosomes or scaffolds) for which to perform variant calling.
3. **Basename**: The base name used for output files (typically the sample name).
4. **Chromosome (chr)**: The current chromosome or scaffold being processed (set dynamically during job execution).

---

## 📝 Script Overview

The script performs the following steps:

**Variant Filtering**:
   - Removes low-quality variants (using `MQFILTER`, `LowQual`, `Repeat` flags).
   - Discards indels and other unwanted variants.

**SNP Extraction**:
   - Extracts only SNPs (Single Nucleotide Polymorphisms) from the VCF, removing indels.

**Filter for Missing Data**: TODO
   - Filters out variants with missing genotypes using `vcftools` (removes positions with missing data).
   - Generates a position list and a VCF file without missing data.

---

## 🖥️ How to Run

To execute the script on the cluster, use the following command:

```
sbatch filters_LowQual.sh
```

## 📤 example output :

VCF with all clean position : ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat.vcf.gz

VCF with SNP clean position : ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat_SNP.vcf.gz

VCF with SNP clean position with no NA : ${Basename}_GATK_TAG_Flowqual_Noindels_Norepeat_NoNa.vcf.gz

