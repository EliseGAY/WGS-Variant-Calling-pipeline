#####===============================#
##### WGS-Variant-Calling-pipeline
##### Elise GAY  - MNHN - EPHE
##### Romuald Laso-Jadart - MNHN - EPHE
#####================================#

This pipeline aim to run a variant calling on Whole Genome Sequencing data from the quality check of Fastq file to the VCF filtering.

All scripts are adapted to be run on supercalculator with SLURM scheduler.

Each folder contains generical scripts and a readme with information to run the step. 

- 1_Quality_Control : FASTQ quality check with fastqc
- 2_Trimming_Fastq : Read trimming with trimmomatic
- 3_Mapping : mapping and pre-process bam file with BWA, Picard and Samtools
- 4_Variant_Calling_GATK : Variant Calling with GATK
- 5_gVCF_Filters : Filtering of gVCF with GATK
- 6_DP_NA_Filters : Vizualisatin and filters for sequencing depth and genotype frequencies on R.
  
Requierment : 

FastQC / Trimmomatic / Samtools / BWA / Picard / GATK / BCFTools / R
