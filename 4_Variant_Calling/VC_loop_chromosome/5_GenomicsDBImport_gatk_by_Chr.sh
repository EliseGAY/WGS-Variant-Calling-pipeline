#!/usr/bin/bash

#==================#
# INFOS :
#==================#
#'''
#GenomicsDBImport create a database from all the ouputs given by haplotypeCaller and with specified scafold/chr

#
#GenomicsDBImport create a database from all the ouputs given by haplotypeCaller and with specified scafold/chr
#INPUT : GenomicsDBImport merge GVCFs from multiple samples
#------------------------------------------------------------
#- Gvcf inputs :
#       One or more GVCFs produced by in HaplotypeCaller with the `-ERC GVCF` or `-ERC BP_RESOLUTION` settings, containing the samples to joint-genotype.
#
#- sample map file input :
#
#       sample1      sample1.vcf.gz
#       sample2      sample2.vcf.gz
#       sample3      sample3.vcf.gz
#
#- Interval.list input :
#
#       SUPER_1
#       SUPER_2
#       SUPER_3
#       ...
#
#- Create a directory "BD" during the run
#- "Temp directory" has to be prealably created
#
#Ouput :
#-------
#A GenomicsDB workspace
#'''


#==================#
# Load Directories
#==================#
# Interval list
Intervall_list="/travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/DB_VCF/interval.list"
# INPUT SAMPLES MAP FILE
SAMPLES_MAP="/travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/DB_VCF/samples_lists"
# Temporary path
Temp_path="/travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/DB_VCF/DB-TEMP"

#==============================================#
# Run DBImport for each Chr in 'interval.list'
#==============================================#
while read chr;
        do

        cat > ${chr}_genomicdb.sh << EOF
#!/bin/sh
#SBATCH --clusters=mesopsl1
#SBATCH --account=gay
#SBATCH --partition=def
#SBATCH --qos=mesopsl1_def_long
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --job-name=DBImport
#SBATCH --time=96:00:00
#SBATCH -o DBImport${chr}.o
#SBATCH -e DBImport${chr}.e

# IMPORT MODULE
module load gcc/9.2.0
module load samtools/1.10
module load gatk/4.2.0.0

# Run GenomicsDBImport

gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
--sample-name-map ${SAMPLES_MAP} \
--reader-threads 10 \
--genomicsdb-workspace-path /travail/egay/capture_analysis_GWS/Variant_Calling/subset_best_samples/GATK/DB_VCF/${chr}_DB/ \
--tmp-dir ${Temp_path} \
-L ${chr}
EOF
        sbatch ${chr}_genomicdb.sh
done < ${Intervall_list}

