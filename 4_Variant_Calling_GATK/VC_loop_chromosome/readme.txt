#======================#
# 02/2022
# Elise GAY
# Run Variant Calling
#======================#

# Aim : 
#------#
Continue the variant calling by chromosome 
Suitable for large genome
The scripts have to be run one by one independantly

# Input :
#----------#
- sample map file input (tab delimited column , no header) :

       sample1      /your_absolute_path/sample1.vcf.gz
       sample2      /your_absolute_path/sample2.vcf.gz
       sample3      /your_absolute_path/sample3.vcf.gz

- Interval.list input :
       SUPER_1
       SUPER_2
       SUPER_3
       ...

# Methods :
#----------#
5. GenomicsDBImport_gatk_by_Chr.sh
	Run the GenomicDBImport function of GATK : Allow to create a database needed for the genotyping step
	The loop "while" iterate among the chromosome listed in the interval list file 

	Create a directory "BD" during the run
	"Temp directory" has to be prealably created for
	
6. Genotype_gvcf_per_Chr.sh
	Get on final gVCF by chr with all samples combined in each
	Needed some variable (path to temp directory, interval list, genome, the DB database) to fill in directly in the GATK command. 
	Please read the comments in the script

7. VCF_tagging_all_chr.sh
	Tag position with a custom filter : MQ, quality, DP threshold OR mask some region : repeat.
	
# output :
#----------#
from (5.) : A DB folder for each chromosome, no need to get the insight of the folder.
from (6.) : One gVCF for each Scaffold / chromosome provided in the interval_list
from (7.) : One tagged gVCF for each Scaffold / chromosome provided in the interval_list
