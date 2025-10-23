# credit to Romuald Laso-Jadart
# To do 

#!/usr/bin/bash
#SBATCH -p std
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --job-name=gatk
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH -V

# IMPORT MODULE
#module load singularity/4.2.2 / if gatk not installed

chr_list=$(cat chr.list)
ListOfFiles=$(realpath *_step1_haplotypecaller/*.vcf.gz  | sed 's/^/--variant /g' | tr '\n' ' ' )
Genome="/scratch/lasojada/Fourmis/ref_seq/flye_polished_yahs_curation_round2.1.primary.curated.fasta"
DIRECTORY=$PWD

mkdir All
cd All

#### singularity run -B /scratch:/scratch $HOME/gatk_latest.sif gatk CombineGVCFs \ : / if gatk not installed

Gatk --java -jar CombineGVCFs \
-R ${Genome} \
${ListOfFiles} \
# -L $CHR \ : if you want to do that by chromosome
-O Cata_All.vcf.gz

EOF

sbatch combine_$CHR.sh
cd $DIRECTORY
done
