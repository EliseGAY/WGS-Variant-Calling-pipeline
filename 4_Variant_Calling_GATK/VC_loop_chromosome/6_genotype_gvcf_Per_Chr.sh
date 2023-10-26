#!/usr/bin/bash

#==================#
# Load Directories
#==================#
# Intervall list
Intervall_list="/travail/egay/capture_analysis_GWS/Variant_Calling/GATK/VCF_214samples/interval.list"

# loop the script by chromosome
while read chr;
 do
        cat > step3_${chr}_genotype.sh << EOF
#!/bin/sh
#SBATCH --clusters=mesopsl1
#SBATCH --account=gay
#SBATCH --partition=def
#SBATCH --qos=mesopsl1_def_long
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --job-name=genotyping${chr}
#SBATCH --time=56:00:00
#SBATCH -o genotyping${chr}.o
#SBATCH -e genotyping${chr}.e

# IMPORT MODULE
module load gcc/9.2.0
module load samtools/1.10
module load gatk/4.2.0.0


#==================#
# Run genotyping
#==================#

mkdir /travail/egay/Whole_Genome_analysis_GWS/GATK/GATK_21samples/VCF_Final/${chr}_temp/

gatk --java-options "-Xmx100g" GenotypeGVCFs \
-R /travail/egay/Genome_Reference/CarCar2.pri.cur.20210205.fasta \
-V gendb:///travail/egay/Whole_Genome_analysis_GWS/GATK/GATK_21samples/DB_VCF/${chr}_DB/ \
-O /travail/egay/Whole_Genome_analysis_GWS/GATK/GATK_21samples/VCF_Final/${chr}_GATK.vcf.gz \
--tmp-dir /travail/egay/Whole_Genome_analysis_GWS/GATK/GATK_21samples/VCF_Final/${chr}_temp/ \
--include-non-variant-sites true \
-L ${chr}
EOF
        sbatch step3_${chr}_genotype.sh
done < ${Intervall_list}
