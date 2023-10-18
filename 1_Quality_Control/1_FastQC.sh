#!/usr/bin/bash

#=============#
# Directories
#=============#

# INPUT FASTQ DIRECTORY
# Put all the fastq files you want to analyze in the 'RAWDATA/' directory
Dir="/work/egay/RAWDATA/"

#=============#
# fastqc
#=============#

# GET FASTQ FILES LIST
fastq_files=$(ls ${Dir} | grep "fastq.gz")

# PRINT FASTQ FILE LIST TO CHECK
echo $fastq_files

# LOOP ON FASTQ FILE NAME + CREATION OF *_fastqc.sh SCRIPT FOR EACH FASTQ FILE + CREATE A JOB FOR EACH *_fastqc.sh ON THE CLUSTER
for file in $fastq_files
do
    cat > ${file}_fastqc.sh << EOF
#!/bin/bash
#SBATCH -V
#SBATCH --nodes=2
#SBATCH --time=01:30:00
#SBATCH --tasks-per-node=1
#SBATCH --mem-per-cpu=5G
#SBATCH -o fastqc.out
#SBATCH -e fastqc.err
#SBATCH -J fastqc

# IMPORT MODULE
module load bioinfo/FastQC_v0.11.7

fastqc -o fastqc_trim -t 16 $Dir$file

EOF
    sbatch ${file}_fastqc.sh
done
