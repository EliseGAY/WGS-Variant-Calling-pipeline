#!/bin/bash

# sample basenames present in bam file names and bam folders
Basename_samples="GN1060_S4_L001 GN1061_S5_L001 GN10687_S7_L001 GN10688_S8_L001 GN13308_S9_L001 GN13309_S10_L001 GN13336_S11_L001 GN13337_S15_L001"

for name in $Basename_samples
do
        #==================#
        # Load Directories
        #==================#

        # Absolute Path of your local directory (where the script is launch)
        Local_PATH="/travail/egay/capture_analysis_GWS/Variant_Calling/GATK/"

        # Input bam file from mapping step
        BAM="/travail/egay/capture_analysis_GWS/MAPPING/mapping_BAITS/${name}_mapping/${name}.sorted.duplicates_BAITS.bam.gz"

        # Interval.list file
        # Interval list file contains all contig/chr/scaffold you want to do the calling on. Put all your contig name if you want to do the calling in all the genome
        Intervall_list="/travail/egay/capture_analysis_GWS/Variant_Calling/GATK/interval.list"

        #'- Interval.list input example:

        #SUPER_1
        #SUPER_2
        #SUPER_3
        #'

        # PATH to the reference genome fasta file : 
		# Note : If genome index is not present along with the fasta file run the following command line :
			# java -jar /travail/egay/software/picard.jar CreateSequenceDictionary R=CarCar2.pri.cur.20210205.fasta
        Genome="/travail/egay/Genome_Reference/CarCar2.pri.cur.20210205.fasta"

        # Output path and name (folder name= same as the mkdir command line below)
        Output="/travail/egay/capture_analysis_GWS/Variant_Calling/GATK/${name}_step1_variantcalling/${name}_gatk.vcf.gz"

    # Create and go in the output gatk directory for the sample x
    mkdir ${name}_step1_variantcalling
    cd ${name}_step1_variantcalling
    cat > ${name}_haplotypecaller.sh << EOF
#!/bin/sh
#SBATCH --clusters=mesopsl1
#SBATCH --account=gay
#SBATCH --partition=def
#SBATCH --qos=mesopsl1_def_long
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --job-name=gatk_${name}
#SBATCH --time=24:00:00
#SBATCH -o gatk_${name}.o
#SBATCH -e gatk_${name}.e

# IMPORT MODULE
module load gcc/9.2.0
module load samtools/1.10
module load gatk/4.2.0.0

# Run HaplotypeCaller / BP_RESOLUTION = option which print all the bases, not only variant bases

gatk --java-options "-Xmx10g" HaplotypeCaller -ERC BP_RESOLUTION \
-R ${Genome} \
-I ${BAM} \
-O ${Output} \
--intervals ${Intervall_list} \
--do-not-run-physical-phasing true
EOF
    sbatch ${name}_haplotypecaller.sh
    cd ${Local_PATH}
done
