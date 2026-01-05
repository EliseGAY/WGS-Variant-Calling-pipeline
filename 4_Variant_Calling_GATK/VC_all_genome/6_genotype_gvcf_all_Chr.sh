#!/usr/bin/bash
#SBATCH --clusters=mesopsl1
#SBATCH --account=gay
#SBATCH --partition=def
#SBATCH --qos=mesopsl1_def_long
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --job-name=genotyping
#SBATCH --time=56:00:00
#SBATCH -o genotyping.o
#SBATCH -e genotyping.e

# IMPORT MODULE
module load gcc/9.2.0
module load samtools/1.10
module load gatk/4.2.0.0

#==================#
# Load Directories
#==================#
# don't know why but I can't manage to create variable with path and input / output. I have to write the path directly in the GATK command line

# Just as a reminder, those are paths to load in the GATK command #
# folder of DBImport/DB : (with three /// at the begining) 
# Do not use if you use combinegVCF
DB="///travail/egay/capture_analysis_GWS/Variant_Calling/GATK/DB_Import/DB/"

# Temporary path (optional if you have a OOM kill event)
Temp_path="/travail/egay/capture_analysis_GWS/Variant_Calling/GATK/VCF_214samples/Temp/"

# PATH to the reference genome fasta file (masked or not masked)
Genome="/travail/egay/Genome_Reference/CarCar2.pri.cur.20210205.fasta"

# interval list (to be mentioned even if you do your calling on all your chr) : if you need to subset your vcf 
interval_list="/travail/egay/capture_analysis_GWS/Variant_Calling/GATK/VCF_214samples/interval.list" 

Output="/travail/egay/capture_analysis_GWS/Variant_Calling/GATK/VCF_214samples/VCF_214_samples_interval.vcf.gz"

# create VCF folder outputs
mkdir "/travail/egay/capture_analysis_GWS/Variant_Calling/GATK/VCF_214samples/Temp/"

#==================#
# Run genotyping
#==================#

gatk --java-options "-Xmx100g" GenotypeGVCFs \
-R "/travail/egay/Genome_Reference/CarCar2.pri.cur.20210205.fasta" \
-V gendb:///travail/egay/capture_analysis_GWS/Variant_Calling/GATK/DB_Import/DB/ \ # if you used CombinegVCF just give the absolute path of the all.VCF.gz
-O "/travail/egay/capture_analysis_GWS/Variant_Calling/GATK/VCF_214samples/VCF_214_samples_interval.vcf.gz" \
# --tmp-dir "/travail/egay/capture_analysis_GWS/Variant_Calling/GATK/VCF_214samples/Temp/" \ # if needed
--include-non-variant-sites true \
--sample-ploidy 1 \
# -L "/travail/egay/capture_analysis_GWS/Variant_Calling/GATK/VCF_214samples/interval.list" # if you need to subset your vcf


