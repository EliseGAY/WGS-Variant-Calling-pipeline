#!/usr/bin/bash

# Write your samples prefix int the list
Basename_samples="GN17608_S40_L002 GN17608_S52_L003 GN17768_S41_L002 GN1434_S33_L002 GN17768_S53_L003 GN1434_S45_L003 GN17774_S42_L002 GN1436_S34_L002 GN17774_S54_L003 GN1436_S46_L003 GN17775_S43_L002 GN17525_S35_L002 GN17775_S55_L003 GN17525_S47_L003 GN18404_S36_L002 GN17582_S38_L002 GN17582_S50_L003 GN18404_S48_L003 GN18411_S37_L002 GN17597_S39_L002 GN18411_S49_L003"

for name in $Basename_samples
do

	#==================#
	# Load Directories
	#==================#
	
	# Absolute Path of your local directory (where the script is launch)
	Local_PATH="/work/egay/white_shark_project/samples_mapping/MAPPING/"
	
	# Picard PATH ".jar" file
	Picard="/usr/local/bioinfo/src/picard-tools/picard-2.20.7/picard.jar"
	
	# INPUT BAM
	BAM="/save/egay/white_shark_project/MAPPING_save/${name}_mapping/${name}.sorted.duplicates.bam.gz"

	# Output BAM
	BAM_output="/save/egay/white_shark_project/MAPPING_save/${name}_mapping/${name}.sorted.duplicates.GROUP.bam.gz"
	

	# PATH to the reference genome fasta file
	Genome="/work/egay/white_shark_project/samples_mapping/MAPPING/index_genome_male/CarCar2.pri.cur.20210205.fasta"
	
	# Create and go in the output mapping directory for the sample x
	# mkdir ${name}_mapping # folder already created
    cd ${name}_mapping

	
	#======================================================#
	# Mapping : write and launch a script for each sample
	#======================================================#
    cat > ${name}_mapping.sh << EOF
#!/bin/bash
#SBATCH -V
#SBATCH -J readgroup
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=10G
#SBATCH --tasks-per-node=2
#SBATCH --nodes=1

# IMPORT MODULE
module load bioinfo/picard-2.20.7
module load bioinfo/gatk-4.2.0.0
module load bioinfo/samtools-1.9

# Add samples ID in BAM file
java -jar ${Picard} AddOrReplaceReadGroups \
I= ${BAM} \
O= ${BAM_output} \
RGID=${name} \
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=unit1 \
RGSM=${name}

# Index Bam file (needed for Variant Calling)
samtools index ${BAM_output}

EOF
    sbatch ${name}_mapping.sh # Launch the script to map the sample x
    cd ${Local_PATH} # return to the location
done
