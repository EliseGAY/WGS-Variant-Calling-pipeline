# Run TRIMMOMATIC  

## Author  
**Elise GAY**  
📅 *02/2022*  

---  

## 📌 Aim  
Run BWA

## 📂 Input  
- **R1 and R2 Trimed.FASTQ.gz files**  
- **One reference genome in fasta format**
- **Basename_samples** : list of basename of your samples

## 🛠 Methods

1. Genome indexation with bwa (needed before the mapping)
`sbatch index.sh `

2. Mapping : 

`sh mapping.sh`

Fill the several variables at the begining of the script (adapt variable declaration if necessary):
`sample_name` = 
`Local_PATH`
`Picard`
`Temp_duplicates_folder`
`DIR_samples`
`fastq_R1`
`fastq_R2`
`Genome`

	The mapping.sh script run a loop on sample list and launch a job in the cluster for every iteration (every sample)
	
	- Mapping with BWA
	- sorting by coordinates with samtools
	- mark duplicated read with PICARD
	- get several statistic with samtools 

(optional) 3. Add_sample_ID.sh
only if : -R (add tag "RGID" on reads : needed for GATK) was not used in the mapping.sh script

```sh trim.sh```

Note: The command
```sbatch script.sh```

is included inside trim.sh itself (line 53) to run one SLURM script per sample.

## 📤 Output
Four FASTQ.gz files:

- paired_R1.fastq.gz
- paired_R2.fastq.gz
- unpaired_R1.fastq.gz
- unpaired_R2.fastq.gz

One trimming report per sample




#==================#
# 02/2022
# Elise GAY
# Run BWA
# please inform the authors before sharing
#==================#

# Aim : 
#------#
Run mapping of illumina paired-end seuqences on whole reference genome

# Input :
#----------#
One reference genome in fasta format
R1 and R2 trimmed fastq.gz files
Basename_samples : list of basename of your samples

# Methods :
#----------#
Script adapted to run on the PSL cluster 

1. index.sh 
	genome indexation with bwa (needed before the mapping)

2. mapping.sh
	Fill the several variables at the begining of the script (sample_name, Local_PATH, Picard, Temp_duplicates_folder, DIR_samples, fastq_R1, fastq_R2, Genome)
	The mapping.sh script run a loop on sample list and launch a job in the cluster for every iteration (every sample)
	
	- Mapping with BWA
	- sorting by coordinates with samtools
	- mark duplicated read with PICARD
	- get several statistic with samtools 

(optional) 3. Add_sample_ID.sh
only if : -R (add tag "RGID" on reads : needed for GATK) was not used in the mapping.sh script

# output :
#----------#
# indexed final bam files : "${name}.sorted.duplicates.bam.gz"
# Statistics on final bam : ${name}.bam.stats and ${name}_metrics_coverage.txt
