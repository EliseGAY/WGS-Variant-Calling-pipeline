#!/usr/bin/bash

# Write your samples prefix in the list
Basename_samples="GN17526_S5_L004 GN18410_S9_L004  GN18417_S6_L004 GN18424_S7_L004 GN18426_S8_L004 GN18972_S1_L004 GN18973_S2_L004 GN18983_S3_L004 GN18986_S4_L004"

#Basename_samples="GN17525_S1_L002 GN17768_S2_L002"
#Basename_samples="GN1434_S33_L002 GN1434_S45_L003 GN1436_S34_L002 GN1436_S46_L003 GN17774_S42_L002 GN17774_S54_L003 GN17775_S43_L002 GN17775_S55_L003 GN17525_S35_L002 GN17525_S47_L003 GN17608_S40_L002 GN17608_S52_L003 GN17582_S38_L002 GN17582_S50_L003 GN17597_S39_L002 GN17597_S51_L003 GN17768_S41_L002 GN17768_S53_L003 GN18404_S36_L002 GN18404_S48_L003 GN18411_S37_L002 GN18411_S49_L003 GN17525_S47_L003 GN17525_S35_L002"

for name in $Basename_samples
do
	#==================#
	# Load Directories
	#==================#


	
	# Samples names for BAM file
	sample_name=${name::-5}
	# Absolute Path of your local directory (where the script is launch)
	Local_PATH="/work/egay/white_shark_project/samples_mapping/MAPPING/"
	
	# Picard PATH ".jar" file
	Picard="/usr/local/bioinfo/src/picard-tools/picard-2.20.7/picard.jar"
	
	# Temp folder to stock temp file in markduplicates
	   Temp_duplicates_folder=${Local_PATH}${name}_mapping/TMP

	# INPUT FASTQ TRIMMED DIRECTORY
	DIR_samples="/travail/egay/Whole_Genome_analysis_GWS/TRIMMING/"
	# Absolute path of fastq_trimmed_R1/R2.fastq.gz files
    fastq_R1=${DIR_samples}${name}_R1.trim.paired.fastq.gz
    fastq_R2=${DIR_samples}${name}_R2.trim.paired.fastq.gz
	
	# PATH to the reference genome fasta file
	Genome="/work/egay/white_shark_project/samples_mapping/MAPPING/index_genome_male/CarCar2.pri.cur.20210205.fasta"
	
	# Create and go in the output mapping directory for the sample x
    mkdir ${name}_mapping
    cd ${name}_mapping
    mkdir ${Temp_duplicates_folder}
	
	#======================================================#
	# Mapping : write and launch a script for each sample
	#======================================================#
    cat > ${name}_mapping.sh << EOF
#!/bin/sh
#SBATCH --clusters=mesopsl1
#SBATCH --account=gay
#SBATCH --partition=def
#SBATCH --qos=mesopsl1_def_long
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --job-name=mapp_${name}
#SBATCH --time=20:00:00

# IMPORT MODULE

module load gcc/9.2.0
module load samtools/1.10

# MAPPING STEP
/travail/egay/software/bwa-0.7.17/bwa mem -t 30 -M -R "@RG\tID:${name}_1\tSM:${sample_name}" "/travail/egay/Genome_Reference/CarCar2.pri.cur.20210205.fasta" ${fastq_R1} ${fastq_R2} | samtools view -bS -1 -h > ${name}.bam.gz

echo "mapping done"

# SORTING STEP
samtools sort -l 6 -o ${name}.sorted.bam.gz -O bam -@ 20 ${name}.bam.gz
echo "sorting done"

# MARK DUPLICATES
java -Djava.io.tmpdir=${Temp_duplicates_folder} -jar /travail/egay/software/picard.jar MarkDuplicates \
I=${name}.sorted.bam.gz \
M=metrics_duplicates.txt \
O=${name}.sorted.duplicates.bam.gz \
COMPRESSION_LEVEL=5
echo "mark duplicates done"

# GEt statistics
samtools index ${name}.sorted.duplicates.bam.gz
samtools stats -@10 ${name}.sorted.duplicates.bam.gz >> ${name}.bam.stats
samtools coverage -o ${name}_metrics_coverage.txt ${name}.sorted.duplicates.bam.gz

EOF
    sbatch ${name}_mapping.sh # Launch the script to map the sample x
    cd ${Local_PATH} # return to the location
done
