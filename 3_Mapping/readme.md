# Run TRIMMOMATIC  

## Author  
**Elise GAY**  
📅 *02/2022*  

---  

## 📌 Aim  
Run BWA

`1_index.sh`
`1_mapping.sh`
(optional) `3_Add_samples_ID.sh`

## 📂 Input  
- **R1 and R2 Trimed.FASTQ.gz files**  
- **One reference genome in fasta format**
- **Basename_samples** : list of basename of your samples

## 🛠 Methods

*️⃣ Genome indexation with bwa (needed before the mapping)

run : 

`sbatch index.sh `

*️⃣ Mapping : 

Fill the several variables at the begining of the script (adapt variable declaration if necessary):

`sample_name` = your list of sample
`Local_PATH` = Root of working dir
`Picard` = Path to Jar file
`Temp_duplicates_folder` = Path to temp folder 
`DIR_samples` = Path to bam files
`fastq_R1` = R1 name
`fastq_R2` = R2 name
`Genome`  = Path to the folder containing the indexed genome

run :

`sh mapping.sh`

*️⃣ Steps by steps : 

`sh mapping.sh` : Loop over the sample list and launch as much job as samples number

- Mapping with BWA

- sorting by coordinates with samtools

- mark duplicated read with PICARD

- get several statistic with samtools 



## 📤 Output

Indexed genome file

Indexed final bam files : `${name}.sorted.duplicates.bam.gz` `${name}.sorted.duplicates.bam.bai` 

Statistics on final bam : `${name}.bam.stats` and `${name}_metrics_coverage.txt`




