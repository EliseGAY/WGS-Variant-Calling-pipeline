#======================#
# 02/2022
# Elise GAY
# Run Variant Calling
# please inform the authors before sharing
#======================#

# Aim : 
#------#
Continue the variant calling on the whole genome. Create a gVCF containing all sample and all chromosome 
Not suitable for large genome 

# Input :
#----------#
- A folder created in the previous step "samples_step1_variantcalling" with the haplotype_caller.sh script

- sample map file input (tab delimited column , no header) :

       sample1      /your_absolute_path/sample1.vcf.gz
       sample2      /your_absolute_path/sample2.vcf.gz
       sample3      /your_absolute_path/sample3.vcf.gz

- Interval.list input :
       SUPER_1
       SUPER_2
       SUPER_3
       ...

- Create a directory "BD" during the run
- "Temp directory" has to be prealably created

# Methods :
#----------#

5. GenomicsDBImport_gatk_all_Chr.sh
	Run the GenomicDBImport function of GATK : Allow to create a database needed for the genotyping step
	
6. Genotype_gvcf_all_Chr.sh
	Get the final gVCF with all chromosome and all samples combined
	Needed some variable (path to temp directory, interval list, genome, the DB database) to fill in directly in the GATK command. 
	Please read the comments in the script

7. VCF_tagging_all_chr.sh
	Tag position with a custom filter : MQ, quality, DP threshold OR mask some region : repeat.
	

# output :
#----------#
from (5.) : A DB folder , no need to go get the insight of the folder
from (6.) : A final GVCF 
from (7.) : A Tagged GVCF with position tagged for quality filters
