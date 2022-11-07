#!/bin/sh
#SBATCH --clusters=mesopsl1
#SBATCH --account=gay
#SBATCH --partition=def
#SBATCH --qos=mesopsl1_def_long
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --job-name=combine_${name}
#SBATCH --time=48:00:00
#SBATCH -o combine.o
#SBATCH -e combine.e

# IMPORT MODULE
module load gcc/9.2.0
module load samtools/1.10
module load gatk/4.2.0.0


#==================#
# Load Directories
#==================#
    # Input of intervall list
    Intervall_list="/travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/interval.list"
        #'- Interval.list input :

        #SUPER_1
        #SUPER_2
        #SUPER_3
        #...
        #'

    # PATH to the reference genome fasta file (masked or not masked) : java -jar /travail/egay/software/picard.jar CreateSequenceDictionary R=CarCar2.pri.cur.20210205.fasta
    Genome="/travail/egay/Genome_Reference/CarCar2.pri.cur.20210205.fasta"
    # Output path and name
    Output="/travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/Combine_GVCF/Combine57sp_MapWG_GATK_GWS.vcf.gz"

#==================#
# run CombineGVCF
#==================#
gatk --java-options "-Xmx10g" CombineGVCFs \
 -R ${Genome} \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN13308_step1_variantcalling/GN13308_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN13336_step1_variantcalling/GN13336_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN13337_step1_variantcalling/GN13337_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN13338_step1_variantcalling/GN13338_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN1426_step1_variantcalling/GN1426_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN1429_step1_variantcalling/GN1429_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN17544_step1_variantcalling/GN17544_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN17559_step1_variantcalling/GN17559_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN17588_step1_variantcalling/GN17588_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN17767_step1_variantcalling/GN17767_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18328_step1_variantcalling/GN18328_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18330_step1_variantcalling/GN18330_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18396_step1_variantcalling/GN18396_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18397_step1_variantcalling/GN18397_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18402_step1_variantcalling/GN18402_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18405_step1_variantcalling/GN18405_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18409_step1_variantcalling/GN18409_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18410_step1_variantcalling/GN18410_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18411_step1_variantcalling/GN18411_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18412_step1_variantcalling/GN18412_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18419_step1_variantcalling/GN18419_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18422_step1_variantcalling/GN18422_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18424_step1_variantcalling/GN18424_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18425_step1_variantcalling/GN18425_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18429_step1_variantcalling/GN18429_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18539_step1_variantcalling/GN18539_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18541_step1_variantcalling/GN18541_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18546_step1_variantcalling/GN18546_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18937_step1_variantcalling/GN18937_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18940_step1_variantcalling/GN18940_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18941_step1_variantcalling/GN18941_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18942_step1_variantcalling/GN18942_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18943_step1_variantcalling/GN18943_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18944_step1_variantcalling/GN18944_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18945_step1_variantcalling/GN18945_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18946_step1_variantcalling/GN18946_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18947_step1_variantcalling/GN18947_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18954_step1_variantcalling/GN18954_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18966_step1_variantcalling/GN18966_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18967_step1_variantcalling/GN18967_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18968_step1_variantcalling/GN18968_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18969_step1_variantcalling/GN18969_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18971_step1_variantcalling/GN18971_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18972_step1_variantcalling/GN18972_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18973_step1_variantcalling/GN18973_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18974_step1_variantcalling/GN18974_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18976_step1_variantcalling/GN18976_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18978_step1_variantcalling/GN18978_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18981_step1_variantcalling/GN18981_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18982_step1_variantcalling/GN18982_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18983_step1_variantcalling/GN18983_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18985_step1_variantcalling/GN18985_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18986_step1_variantcalling/GN18986_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18989_step1_variantcalling/GN18989_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN18990_step1_variantcalling/GN18990_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN2425_step1_variantcalling/GN2425_gatk.vcf.gz \
 --variant /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/GN2_step1_variantcalling/GN2_gatk.vcf.gz \
 -O ${Output} \
 -L ${Intervall_list}
