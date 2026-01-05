## Variant Calling

## Author  
**Elise GAY**  
📅 *02/2022*  

---  

### 📌 Aim  

Create a gVCF file from all chr at once

Two options with GATK : `CombineGVCF` or `DBimport`

### STEP 1 option 1 : DBIMPORT 

📂 Input (see details format in the sh script)

Fill the variables in the corresponding script : 

`samples_step1_variantcalling` : Folder created the step before

`sample map file` tab delimited column , no header :

```
sample1      /your_absolute_path/sample1.vcf.gz

sample2      /your_absolute_path/sample2.vcf.gz
  
sample3      /your_absolute_path/sample3.vcf.gz
   
```

`Interval list` : file with chr name

`DB` : Path to the DB database

```
SUPER_1

SUPER_2

SUPER_3
```

`Temp directory` = Path to temp directory


🛠 Methods

*️⃣ Run the GenomicDBImport function of GATK : create a database needed for the genotyping step

`sh GenomicsDBImport_gatk_all_Chr.sh`

📤 Outputs :

- A DB folder, no need to go get the insight of the folders

### STEP 1 option 2 : CombineGVCF

📂 Input :

list of chr in `chr.list` file
All vcf in one folder `vcf_input` using to search for vcf `vcf_input/*.vcf.gz`
Genome_path `Genome` 
current dir `PWD`

🛠 Methods

use the `Combinegvcf` function of GATK

📤 Outputs :
 All vcf

### STEP 2 : Generate the gvcf

`sh Genotype_gvcf_all_Chr.sh`

📤 Output :

- One gVCF with genotypes

### STEP 3 : Tag the position according to quality criteria

📂 Input :

Your gVCF path

`repeat` = path to repeat file in bed format

📓Note :  

To obtain the repeat coordinate with RepeatMasker run the script `Repeats_detection.sh` (check for corresponding repo)

Briefly, to get the repeat bed file from repeatmakser output :

`cut -f5,6,7 'genome.fasta.out' >> Repeats.bed`

`bedtools sort -i Repeats.bed >> Repeats_sorted.bed`

`index bed file with GATK IndexFeatureFile -I repeats.bed`

🛠 Methods :

`sh VCF_tagging_all_chr.sh`

📤 Outputs :

- One gVCF with genotypes and quality tags 






