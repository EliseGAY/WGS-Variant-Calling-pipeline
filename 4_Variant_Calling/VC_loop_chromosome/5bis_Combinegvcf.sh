List_Chr=$(cat intervals.list)

for chr in ${List_Chr}
	do 
	cd /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/All
	cat > ${chr}_combine.sh  << EOF
#!/bin/bash
#SBATCH --partition=type_2
#SBATCH -V
#SBATCH -J combinegvcfs
#SBATCH --time=50:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G

# IMPORT MODULE
module load biology
module load userspace/tr17.10
module load java-JDK-OpenJDK/11.0.9
module load gatk/4.2.0.0

 gatk CombineGVCFs \
   -R "/mnt/beegfs/rlasojadart/Rats/ref_seq/mRatBN7_20st_XYmito_masked_linear.fasta" \
   --variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/ANR0085_S11_variantcalling/ANR0085_S11_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/ANR0086_S14_variantcalling/ANR0086_S14_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/ANR0087_S17_variantcalling/ANR0087_S17_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/MEN307_S5_variantcalling/MEN307_S5_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PAR005_S10_variantcalling/PAR005_S10_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PAR009_S22_variantcalling/PAR009_S22_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PAR011_S21_variantcalling/PAR011_S21_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PAR057_S2_variantcalling/PAR057_S2_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PAR068_S3_variantcalling/PAR068_S3_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PAR083_S23_variantcalling/PAR083_S23_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PAR090_S20_variantcalling/PAR090_S20_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PAR1165_S9_variantcalling/PAR1165_S9_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PAR1920_S13_variantcalling/PAR1920_S13_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PAR1934_S8_variantcalling/PAR1934_S8_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PAR2148_S18_variantcalling/PAR2148_S18_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PAR320_S1_variantcalling/PAR320_S1_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PAR341_S6_variantcalling/PAR341_S6_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PAR379_S16_variantcalling/PAR379_S16_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PAR409_S4_variantcalling/PAR409_S4_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PAR412_S7_variantcalling/PAR412_S7_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PAR431_S15_variantcalling/PAR431_S15_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PAR993_S19_variantcalling/PAR993_S19_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/PZP033_S12_variantcalling/PZP033_S12_gatk.vcf.gz \
--variant /mnt/beegfs/rlasojadart/Rats/Paris/Variant_calling/RBS_S24_variantcalling/RBS_S24_gatk.vcf.gz \
   -O Paris_rats_${chr}.vcf.gz \
   -L ${chr}	

EOF
	sbatch ${chr}_combine.sh
done 
