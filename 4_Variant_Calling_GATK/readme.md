Elise GAY
📅 02/2022

    ⚠️ Please inform the authors before sharing.

📌 Aim

Run variant calling with GATK on a list of samples.
📂 Input (see details format in the sh script)

    List of samples
    BAM files obtained from the mapping step
    Interval_list: File with chromosome (or loci, or scaffold) names
    Reference genome path

🛠 Methods

    haplotypecaller_gatk.sh
    Runs variant calling for each sample provided in the list. The loop launches a job on the cluster for each sample.

    VC_allgenome or VC_loop_chromosome
    Choose to run the script either by chromosome (useful on large genomes) or on the complete genome at once.
    Look at the README in each folder to see how to run the pipeline.

📤 Output

    A folder is created for each sample: samplesX_step1_variantcalling
    A file ${name}_gatk.vcf.gz is created in each folder

#======================#
# 02/2022
# Elise GAY
# Run Variant Calling
# please inform the authors before sharing
#======================#

# Aim : 
#------#
Run variant calling with GATK on list of samples

# Input (see details format in the sh script) :
#-----------------------------------------------#
list of sample
BAM files obtained from mapping step
Intervall_list : file with chromosome (or loci, or scaffold..) names
Reference genome path

# Methods :
#----------#
1. haplotypecaller_gatk.sh
Run the calling for each sample provided in the list. The loop launch a job on the cluster for each samples.

2. VC_allgenome or VC_loop_chromosome : Choose to run script either by chromosome (usefull on large genome) or on the complete genome at once
	Look at the readme in each folder to see how to run the pipeline

# output :
#----------#

A folder is created for each samples "samplesX_step1_variantcalling" 
A "${name}_gatk.vcf.gz" is created in each folder
