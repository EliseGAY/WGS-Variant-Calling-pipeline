# Variant Calling Filtering Pipeline

## Author: Elise GAY  
## Date: 02/2022

This pipeline filters **gVCF** files based on **DP (depth) range** and **genotype frequencies**. It is designed to work with gVCF files generated using **GATK**. An example of its application on a rat dataset (Y haploid chromosome) is available in the corresponding HTML file.

---

## 🛠️ Prerequisites

This pipeline relies on **R** functions located in the "R_functions_ploidy1" or "R_functions_ploidy2" folders. The primary difference between these two is the genotype search handling, where `"."` or `"./."` are treated differently.

You will need an R environment with the necessary packages installed to run the functions for filtering.

---

## 📂 Inputs

1. **gVCF**: The gVCF file generated from previous variant calling steps (e.g., using GATK).
2. **R functions**:
   - Functions located in either the "R_functions_ploidy1" or "R_functions_ploidy2" folder.
   - These functions perform the necessary genotype searches for filtering the VCF files.

---

## 📝 Methods

### 1. **DP Distribution & Genotype Frequencies**
   - The pipeline first describes the **DP (depth) distribution** and **genotype frequencies** in the provided gVCF file.
   - It then applies filters based on **DP range** and **genotype frequency thresholds** using the provided R functions.

### 2. **Filtering with R Functions**
   - The functions filter the gVCF file based on the thresholds defined for DP range and genotype frequencies.
   - The functions used for filtering depend on whether the gVCF is haploid or diploid.

### Script Workflow:
1. Review the results from the example dataset:
   - `.rmd` and `.html` files contain the results of running the pipeline on the rat Y chromosome dataset.
   
2. To run the pipeline on your own dataset, you should adapt the provided `.r` script to suit your needs.

---

## 📥 How to Run

1. **Prepare Input Files**:
   - Ensure that you have the gVCF files ready for processing.
   - Adjust the input files as required and make sure the correct R function folder ("R_functions_ploidy1" or "R_functions_ploidy2") is used based on your gVCF file type.

2. **Run the R Script**:
   - Adapt the `.r` script for your dataset, setting the desired thresholds for DP and genotype frequencies.
   - Execute the R script to filter the gVCF file.

---

## 📁 Output

The output consists of the filtered gVCF file based on the specified DP and genotype frequency thresholds, which will be saved in the **"data"** folder.

The following files will be generated:
- **Filtered VCF**: The gVCF file after applying the DP range and genotype frequency filters.

---

## 📚 Example

For a real-world example, check the `.rmd` and `.html` files, which show the pipeline results for the rat dataset on the Y haploid chromosome.
