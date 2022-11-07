#!/usr/bin/bash

#==================#
# INFOS :
#==================#
'''
INPUT : VCF with all tag added from the 7_VCF_filter.sh script
---------------------------------------------------------------

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  GN1430  GN1434  GN1436  GN17525 GN17582 GN17597 GN17608 GN17768 GN17774 GN17775 GN18404 GN18411
SUPER_9 7205    .       T       .       .       Repeat  DP=2    GT:AD:DP:RGQ    ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       ./.:0:1:0       ./.:0:0:0       ./.:0:0:0       ./.:0:1:0
SUPER_9 7206    .       A       .       .       Repeat  DP=2    GT:AD:DP:FT:RGQ ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:1:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:1:1:DPFILTER:3
SUPER_9 7207    .       A       .       .       Repeat  DP=2    GT:AD:DP:RGQ    ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       ./.:0:1:0       ./.:0:0:0       ./.:0:0:0       ./.:0:1:0
SUPER_9 7208    .       C       .       .       Repeat  DP=2    GT:AD:DP:FT:RGQ ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:1:1:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:1:1:DPFILTER:3
SUPER_9 7209    .       C       .       .       Repeat  DP=2    GT:AD:DP:FT:RGQ ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:1:1:PASS:0  ./.:0:0:PASS:0  ./.:0:0:PASS:0  ./.:1:1:DPFILTER:3
SUPER_9 7210    .       C       .       .       Repeat  DP=2    GT:AD:DP:RGQ    ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       ./.:1:1:0       ./.:0:0:0       ./.:0:0:0       ./.:0:1:0

SUPER_2 10669   .       ATG     A       125.87  MQFILTER        DP=16;ExcessHet=0.2482;FS=0;InbreedingCoeff=0.3763;MLEAC=3;MLEAF=0.3;MQ=24;QD=27.24;SOR=2.833   GT:AD:DP:FT:GQ:PL       ./.:1,0:1:PASS:.:0,0,0  ./.:1,0:1:PASS:.:0,0,0  ./.:1,0:1:PASS:.:0,0,0  ./.:0,3:3:DPFILTER:9:135,9,0    ./.:1,0:1:PASS:.:0,0,0  ./.:2,0:2:DPFILTER:3:0,3,45     ./.:0,0:0:PASS:.:0,0,0  ./.:0,0:0:PASS:.:0,0,0  ./.:2,0:2:DPFILTER:6:0,6,84     ./.:2,0:2:DPFILTER:3:0,3,45     ./.:2,0:2:PASS:.:0,0,0  ./.:1,0:1:DPFILTER:3:0,3,42
SUPER_2 29615   .       A       .       15.64   LowQual;MQFILTER        BaseQRankSum=-1.282;DP=65;ExcessHet=3.01;MLEAC=.;MLEAF=.;MQ=27;MQRankSum=0;ReadPosRankSum=0.431 GT:DP:FT:RGQ    ./.:1:DPFILTER:3        ./.:4:DPFILTER:12       ./.:1:DPFILTER:3        ./.:5:DPFILTER:12       ./.:8:DPFILTER:21       ./.:3:DPFILTER:0        ./.:5:DPFILTER:0        ./.:9:DPFILTER:0        ./.:6:DPFILTER:0        ./.:8:DPFILTER:0        ./.:7:DPFILTER:1        ./.:5:DPFILTER:15

SUPER_25        1919    .       A       .       .       Repeat  AN=6;DP=90      GT:AD:DP:FT:RGQ 0/0:14:14:PASS:42       0/0:11:12:PASS:6        ./.:4:4:DPFILTER:12     ./.:9:9:DPFILTER:27     ./.:10:11:PASS:0        ./.:4:4:DPFILTER:12     0/0:10:10:PASS:30       ./.:7:7:DPFILTER:21     ./.:1
SUPER_25        1920    .       A       .       .       Repeat  AN=8;DP=94      GT:AD:DP:FT:RGQ 0/0:14:14:PASS:36       0/0:12:12:PASS:33       ./.:6:6:DPFILTER:18     ./.:9:9:DPFILTER:27     0/0:10:11:PASS:3        ./.:6:6:DPFILTER:12     0/0:10:10:PASS:30       ./.:7:7:DPFILTER:21     ./.:1
SUPER_25        1921    .       C       .       .       Repeat  AN=6;DP=89      GT:AD:DP:FT:RGQ 0/0:14:14:PASS:42       0/0:12:12:PASS:33       ./.:5:5:DPFILTER:15     ./.:9:9:DPFILTER:24     ./.:8:9:PASS:0  ./.:4:4:DPFILTER:12     0/0:10:10:PASS:30       ./.:7:7:DPFILTER:21     ./.:0:1:PASS:
SUPER_25        1922    .       C       .       .       Repeat  AN=6;DP=88      GT:AD:DP:FT:RGQ 0/0:14:14:PASS:42       0/0:12:12:PASS:33       ./.:4:4:DPFILTER:12     ./.:9:9:DPFILTER:24     ./.:8:8:DPFILTER:24     ./.:4:5:PASS:0  0/0:10:10:PASS:30       ./.:7:7:DPFILTER:21     ./.:1:1:DPFIL

SUPER_25        11343   .       C       T       3933.04 PASS    AC=12;AF=0.545;AN=22;BaseQRankSum=-1.304;DP=214;ExcessHet=4.5806;FS=3.163;InbreedingCoeff=-0.1748;MLEAC=13;MLEAF=0.542;MQ=60;MQRankSum=0;QD=24.58;ReadPosRankSum=0.411;SOR=0.809        GT:AD:DP:FT:GQ:PL       0/1:8,9:17:PASS:99:293,0,2831/1:0,23:23:PASS:69:794,69,0     ./.:7,2:9:DPFILTER:56:56,0,220  0/1:7,6:13:PASS:99:177,0,253    0/1:5,13:18:PASS:99:386,0,155   0/0:19,0:19:PASS:57:0,57,657    0/1:7,12:19:PASS:99:376,0,235   1/1:0,20:20:PASS:60:767,60,0    1/1:0,18:18:PASS:54:684,54,0    0/1:6,4:10:PASS:99:137,0,221    0/0:32,0:32:PASS:90:0,90,1350        0/1:4,9:13:PASS:99:269,0,128
SUPER_25        14199   .       G       T       165.87  PASS    AC=1;AF=0.045;AN=22;BaseQRankSum=0.707;DP=242;ExcessHet=3.0103;FS=0;InbreedingCoeff=-0.0438;MLEAC=1;MLEAF=0.042;MQ=60;MQRankSum=0;QD=8.29;ReadPosRankSum=0.454;SOR=0.473        GT:AD:DP:FT:GQ:PL       0/0:27,0:27:PASS:78:0,78,1170   0/0:19,0:19:PASS:51:0,51,765 ./.:9,0:9:DPFILTER:21:0,21,315  0/1:14,6:20:PASS:99:178,0,462   0/0:16,0:16:PASS:48:0,48,572    0/0:21,0:21:PASS:63:0,63,805    0/0:23,0:23:PASS:69:0,69,927    0/0:20,0:20:PASS:54:0,54,810    0/0:16,0:16:PASS:48:0,48,564    0/0:21,0:21:PASS:60:0,60,900    0/0:25,0:25:PASS:72:0,72,1080        0/0:24,0:24:PASS:63:0,63,945

SUPER_25        11196   .       G       .       .       PASS    DP=252  GT:AD:DP:RGQ    0/0:20:20:57    0/0:23:23:66    0/0:22:22:66    0/0:19:19:57    0/0:19:19:51    0/0:23:23:60    0/0:29:29:84    0/0:10:10:27    0/0:24:24:66    0/0:13:13:33    0/0:29:29:84    0/0:21:21:51
SUPER_25        11197   .       C       .       .       PASS    DP=255  GT:AD:DP:RGQ    0/0:20:20:54    0/0:23:23:66    0/0:22:22:63    0/0:20:20:60    0/0:19:19:48    0/0:23:23:60    0/0:29:29:84    0/0:11:11:27    0/0:24:24:60    0/0:14:14:36    0/0:29:29:81    0/0:21:21:51

	...
	
Output different vcf files depending on the filters chosen
-----------------------------------------------------------
'''

for chr in $Contigs
        do
        cat > ${chr}_extract_vcf.sh << EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --job-name=${chr}_vcf
#SBATCH --time=24:00:00
#SBATCH -V

module load bioinfo/tabix-0.2.5

#============================#
# Initiate PATH and variable
#============================#
initial_VCF_FilterTag="/travail/egay/capture_analysis/GATK/VCF/${chr}_GATK_filtered.vcf.gz"

# Chose Tag on which to be discarded (change tags if needed)
Filters='MQFILTER\|LowQual\|Repeat'

#=============#
# Run filters
#=============#
# create header
#zgrep "#" /work/rlaso/Sharks/Variant*/All*/VCFs/Filter/${chr}_GATK_filtered.vcf.gz > inizio.txt

# From NoDP.filter.vcf discard lines with "MQFILTER\|LowQual\|Repeat" TAGs. NOINDELS_NOREPEATS
vcftools --gzvcf /travail/egay/GATK/${name}_GATK/${name}_SUPER_Y_gatk_DP10_DP50_Q30.vcf.gz --remove-indels --recode --recode-INFO-all --out

zgrep -v "${Filters} ${initial_VCF_FilterTag} | vcftools --remove-indels --recode --recode-INFO-all

zgrep -v "${Filters} ${initial_VCF_FilterTag} | awk '!(\$5 ~ "^*")' | awk '!(\$5 ~ "^[ACGT*].")' | awk '!(\$4 ~ "^[ACGT*].")' | sed 's/##source=VariantFiltration/##source=VariantFiltration\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGN1430\tGN1434\tGN1436\tGN17525\tGN17582\tGN17597\tGN17608\tGN17768\tGN17774\tGN17775\tGN18404\tGN18411/' | bgzip > ${chr}_GATK_NOINDELS_NOREPEATS.vcf.gz

# For SFS.r
# From filter.vcf discard lines. Discard line with ./. ; Discard lines with "MQFILTER" / "LowQual" AND "Repeat" Tags. Discard Indels. GATK_NOMISSING_NOINDELS.vcf.gz
#zgrep -v "\.\/\." /work/rlaso/Sharks/Variant*/All*/VCFs/Filter/${chr}_GATK_filtered.vcf.gz | grep -v "MQFILTER" | grep -v "LowQual" | grep -v "Repeat" | awk '!(\$5 ~ "^*")' | awk '!(\$5 ~ "^[ACGT*].")' | awk '!(\$4 ~ "^[ACGT*].")' | sed 's/##source=VariantFiltration/##source=VariantFiltration\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGN1430\tGN1434\tGN1436\tGN17525\tGN17582\tGN17597\tGN17608\tGN17768\tGN17774\tGN17775\tGN18404\tGN18411/' | bgzip > ${chr}_GATK_NOMISSING_NOINDELS.vcf.gz ### stampo linee senza missing e senza indels/triallelic

# For sNMF (PCA)
# From GATK_NOMISSING_NOINDELS.vcf.gz discard line with "." in ALT allele. NOMISSING_NOINDELS_NOREF
#zcat ${chr}_GATK_NOMISSING_NOINDELS.vcf.gz | awk '!(\$5 ~ "\\\.")' | bgzip > ${chr}_GATK_NOMISSING_NOINDELS_NOREF.vcf.gz

# For clean_vcf_FST
# From filter.vcf. discard MQFILTER; LowQual ; Repeat. Discard line with indels and "." (ref) in alt allele. NOREF_NOINDELS
#zgrep -v "LowQual" /work/rlaso/Sharks/Variant*/All*/VCFs/Filter/${chr}_GATK_filtered.vcf.gz | grep -v "MQFILTER" | grep -v "Repeat" |  awk '!(\$5 ~ "^*")'  | awk '!(\$5 ~ "\\\.")'  | awk '!(\$5 ~ "^[ACGT*].")' | awk '!(\$4 ~ "^[ACGT*].")' | sed 's/##source=VariantFiltration/##source=VariantFiltration\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGN1430\tGN1434\tGN1436\tGN17525\tGN17582\tGN17597\tGN17608\tGN17768\tGN17774\tGN17775\tGN18404\tGN18411/' | bgzip > ${chr}_GATK_NOREF_NOINDELS.vcf.gz ### elimina REF==ALT

EOF

sbatch ${chr}_extract_vcf.sh

done

'''
bcftools view -H SUPER_1_GATK_NOINDELS_NOREPEATS.vcf.gz | head
SUPER_1 1       .       A       .       .       PASS    DP=7    GT:AD:DP:RGQ    0/0:1:1:3       ./.:1:2:0       ./.:0:0:0       ./.:0:0:0       0/0:1:1:3       ./.:0:0:0       0/0:1:1:3       ./.:0:0:0       ./.:0:0:0       ./.:0:0:0       0/0:1:1:3       0/0:1:1:3
SUPER_1 44731   .       T       .       .       PASS    DP=97   GT:AD:DP:RGQ    0/0:10:10:30    0/0:13:13:39    0/0:4:4:12      0/0:13:13:39    0/0:10:10:30    0/0:7:7:21      0/0:6:7:5       0/0:16:17:35    0/0:3:3:9       0/0:3:3:9       0/0:9:9:27      0/0:1:1:3
SUPER_1 44732   .       G       .       .       PASS    DP=98   GT:AD:DP:RGQ    0/0:10:10:30    0/0:13:13:39    0/0:4:4:12      0/0:13:13:33    0/0:10:10:30    0/0:7:7:21      0/0:7:7:18      0/0:17:17:51    0/0:3:3:9       0/0:3:3:9       0/0:10:10:30    0/0:1:1:3
SUPER_1 44733   .       G       .       .       PASS    DP=100  GT:AD:DP:RGQ    0/0:10:10:30    0/0:14:14:41    0/0:5:5:15      0/0:13:13:33    0/0:9:10:14     0/0:7:7:21      0/0:7:8:8       0/0:14:15:29    0/0:3:3:9       0/0:3:3:9       0/0:11:11:33    0/0:1:1:3
SUPER_1 44734   .       G       .       .       PASS    DP=103  GT:AD:DP:RGQ    0/0:10:10:30    0/0:14:14:42    0/0:6:6:18      0/0:13:13:33    0/0:10:10:30    0/0:7:7:21      0/0:8:8:21      0/0:17:17:51    0/0:3:3:9       0/0:3:3:9       0/0:11:11:33    0/0:1:1:3
SUPER_1 44735   .       A       .       .       PASS    DP=103  GT:AD:DP:RGQ    0/0:10:10:30    0/0:14:14:42    0/0:6:6:18      0/0:13:13:33    0/0:10:10:30    0/0:7:7:[W::vcf_parse_info] INFO 'AF' is not defined in the header, assuming Type=String
21      0/0:8:8:21      0/0:17:17:51    0/0:3:3:9       0/0:3:3:9       0/0:11:11:33    0/0:1:1:3
SUPER_1 44736   .       C       .       .       PASS    DP=98   GT:AD:DP:RGQ    0/0:8:8:24      0/0:14:14:42    0/0:6:6:18      0/0:13:13:33    0/0:9:10:14     0/0:7:7:21      0/0:5:6:3       0/0:15:15:45    0/0:3:3:9       0/0:3:3:9       0/0:11:12:21    0/0:1:1:3
SUPER_1 44737   .       C       .       .       PASS    DP=107  GT:AD:DP:RGQ    0/0:10:10:30    0/0:13:14:26    0/0:5:6:2       0/0:14:14:36    0/0:9:10:14     0/0:7:7:21      0/0:9:10:14     0/0:16:17:35    0/0:3:3:9       0/0:3:3:9       0/0:12:12:33    0/0:1:1:3
SUPER_1 46168   .       T       .       .       PASS    DP=243  GT:AD:DP:RGQ    0/0:21:21:60    0/0:39:39:99    0/0:16:16:48    0/0:18:18:54    0/0:9:9:27      0/0:22:22:66    0/0:16:16:48    0/0:21:21:63    0/0:14:14:39    0/0:14:14:39    0/0:35:35:99    0/0:18:18:45
SUPER_1 46169   .       C       .       .       PASS    DP=244  GT:AD:DP:RGQ    0/0:21:21:60    0/0:39:39:99    0/0:16:16:48    0/0:18:18:54    0/0:9:9:27      0/0:22:22:66    0/0:16:16:48    0/0:21:21:63    0/0:14:14:39    0/0:14:14:39    0/0:35:35:99    0/0:19:19:45


bcftools view -H SUPER_1_GATK_NOMISSING_NOINDELS_NOREFNOMISSING_NOINDELS_NOREF.vcf.gz | head
[W::vcf_parse_info] INFO 'AF' is not defined in the header, assuming Type=String
SUPER_1 65437   .       G       A       478.87  PASS    AC=1;AF=0.042;AN=24;BaseQRankSum=-2.027;DP=235;ExcessHet=3.0103;FS=1.603;InbreedingCoeff=-0.0437;MLEAC=1;MLEAF=0.042;MQ=60;MQRankSum=0;QD=19.15;ReadPosRankSum=-1.728;SOR=0.412 GT:AD:DP:GQ:PL  0/0:21,0:21:60:0,60,900 0/0:12,1:13:24:0,24,418 0/0:14,0:14:36:0,36,540      0/0:27,0:27:72:0,72,1080        0/0:21,0:21:60:0,60,900 0/0:19,0:19:54:0,54,810 0/0:20,0:20:54:0,54,810 0/0:23,0:23:66:0,66,990 0/0:19,0:19:57:0,57,760 0/0:15,0:15:45:0,45,592 0/0:18,0:18:51:0,51,765 0/1:9,16:25:99:491,0,301
SUPER_1 91128   .       G       A       2041.05 PASS    AC=6;AF=0.250;AN=24;BaseQRankSum=1.32;DP=230;ExcessHet=0.2035;FS=1.059;InbreedingCoeff=0.5502;MLEAC=6;MLEAF=0.25;MQ=60;MQRankSum=0;QD=27.58;ReadPosRankSum=-0.144;SOR=0.5       GT:AD:DP:GQ:PL  1/1:0,16:16:48:668,48,0 0/0:16,0:16:42:0,42,630 0/1:14,11:25:99:350,0,407    0/0:31,0:31:81:0,81,1215        0/0:15,0:15:42:0,42,630 1/1:0,18:18:54:715,54,0 0/0:17,0:17:48:0,48,720 0/0:18,1:19:15:0,15,630 0/0:21,0:21:63:0,63,809 0/0:14,0:14:42:0,42,512 0/0:23,0:23:60:0,60,900 0/1:3,12:15:80:342,0,80
SUPER_1 92206   .       G       C       9508.98 PASS    AC=22;AF=0.917;AN=24;BaseQRankSum=1.97;DP=276;ExcessHet=3.2034;FS=1.916;InbreedingCoeff=-0.0909;MLEAC=22;MLEAF=0.917;MQ=60;MQRankSum=0;QD=30.05;ReadPosRankSum=0.687;SOR=0.645  GT:AD:DP:GQ:PL  1/1:0,34:34:99:1254,102,0       0/1:8,4:12:99:130,0,264      1/1:0,19:19:57:722,57,0 1/1:0,25:25:75:997,75,0 1/1:0,14:14:42:559,42,0 1/1:0,20:20:60:766,60,0 1/1:0,18:18:54:674,54,0 1/1:0,29:29:87:1143,87,0        1/1:0,24:24:72:967,72,0 0/1:12,5:17:99:158,0,343        1/1:0,39:39:99:1433,117,0       1/1:0,17:17:51:642,51,0
SUPER_1 93119   .       T       A       442.87  PASS    AC=1;AF=0.042;AN=24;BaseQRankSum=-0.429;DP=224;ExcessHet=3.0103;FS=1.685;InbreedingCoeff=-0.0436;MLEAC=1;MLEAF=0.042;MQ=60;MQRankSum=0;QD=18.45;ReadPosRankSum=2.17;SOR=0.346   GT:AD:DP:GQ:PL  0/0:24,0:24:72:0,72,830 0/0:16,0:16:45:0,45,675 0/0:14,0:14:42:0,42,517      0/0:11,0:11:30:0,30,450 0/0:11,0:11:33:0,33,399 0/0:19,0:19:54:0,54,810 0/0:17,0:17:48:0,48,720 0/0:30,0:30:84:0,84,1260        0/0:14,0:14:39:0,39,585 0/0:20,0:20:57:0,57,855 0/0:21,0:21:63:0,63,729 0/1:10,14:24:99:455,0,374
SUPER_1 93120   .       T       A       478.87  PASS    AC=1;AF=0.042;AN=24;BaseQRankSum=-0.74;DP=225;ExcessHet=3.0103;FS=3.982;InbreedingCoeff=-0.0436;MLEAC=1;MLEAF=0.042;MQ=60;MQRankSum=0;QD=20.82;ReadPosRankSum=2.02;SOR=0.209    GT:AD:DP:GQ:PL  0/0:24,0:24:66:0,66,990 0/0:16,0:16:45:0,45,675 0/0:14,0:14:42:0,42,517      0/0:11,0:11:30:0,30,450 0/0:11,0:11:33:0,33,413 0/0:19,0:19:54:0,54,810 0/0:17,0:17:48:0,48,720 0/0:30,0:30:78:0,78,1170        0/0:14,0:14:39:0,39,585 0/0:20,0:20:57:0,57,855 0/0:22,0:22:63:0,63,945 0/1:9,14:23:99:491,0,371
SUPER_1 96422   .       T       A       643.01  PASS    AC=2;AF=0.083;AN=24;BaseQRankSum=2.38;DP=252;ExcessHet=3.2034;FS=0;InbreedingCoeff=-0.092;MLEAC=2;MLEAF=0.083;MQ=60;MQRankSum=0;QD=14.61;ReadPosRankSum=0.561;SOR=0.715 GT:AD:DP:GQ:PL  0/1:13,9:22:99:309,0,375        0/0:23,1:24:57:0,57,775 0/0:20,2:22:17:0,17,684      0/0:25,0:25:72:0,72,1080        0/0:11,0:11:33:0,33,425 0/0:22,0:22:63:0,63,945 0/0:10,0:10:27:0,27,405 0/0:30,0:30:75:0,75,1125        0/0:15,0:15:42:0,42,630 0/0:25,0:25:75:0,75,976 0/0:23,0:23:69:0,69,824 0/1:10,12:22:99:352,0,352
SUPER_1 99000   .       A       G       151.87  PASS    AC=1;AF=0.042;AN=24;BaseQRankSum=1.68;DP=275;ExcessHet=3.0103;FS=5.497;InbreedingCoeff=-0.0435;MLEAC=1;MLEAF=0.042;MQ=60;MQRankSum=0;QD=10.12;ReadPosRankSum=1.9;SOR=1.911      GT:AD:DP:GQ:PL  0/0:32,0:32:93:0,93,1395        0/0:15,0:15:45:0,45,592      0/0:28,0:28:75:0,75,1125        0/1:10,5:15:99:164,0,303        0/0:22,0:22:63:0,63,945 0/0:27,0:27:78:0,78,1170        0/0:20,1:21:47:0,47,755 0/0:32,0:32:78:0,78,1170        0/0:21,0:21:54:0,54,810 0/0:22,0:22:54:0,54,810 0/0:22,0:22:63:0,63,945 0/0:18,0:18:48:0,48,720
SUPER_1 107083  .       C       T       4033.62 PASS    AC=11;AF=0.458;AN=24;BaseQRankSum=-0.251;DP=240;ExcessHet=0.2176;FS=0;InbreedingCoeff=0.4909;MLEAC=11;MLEAF=0.458;MQ=60;MQRankSum=0;QD=29.66;ReadPosRankSum=0.512;SOR=0.772     GT:AD:DP:GQ:PL  0/0:29,0:29:78:0,78,1170        0/1:7,6:13:99:212,0,185      0/1:8,8:16:99:286,0,252 1/1:0,33:33:98:1184,98,0        1/1:0,16:16:48:550,48,0 0/0:22,0:22:57:0,57,855 1/1:1,23:24:30:872,30,0 0/1:5,10:15:99:326,0,164        0/0:15,0:15:42:0,42,630 0/0:18,0:18:51:0,51,765 1/1:1,18:19:15:630,15,0 0/0:17,0:17:51:0,51,676
SUPER_1 107241  .       T       C       352.87  PASS    AC=1;AF=0.042;AN=24;BaseQRankSum=0.274;DP=253;ExcessHet=3.0103;FS=0;InbreedingCoeff=-0.0435;MLEAC=1;MLEAF=0.042;MQ=60;MQRankSum=0;QD=17.64;ReadPosRankSum=-1.102;SOR=1.022      GT:AD:DP:GQ:PL  0/0:31,0:31:81:0,81,1215        0/0:15,0:15:36:0,36,540      0/1:9,11:20:99:365,0,282        0/0:30,0:30:78:0,78,1170        0/0:18,0:18:48:0,48,720 0/0:21,0:21:54:0,54,810 0/0:17,0:17:36:0,36,540 0/0:18,0:18:51:0,51,765 0/0:23,0:23:60:0,60,900 0/0:22,0:22:60:0,60,900 0/0:17,0:17:39:0,39,585 0/0:21,0:21:54:0,54,810
SUPER_1 108383  .       A       T       4322.5  PASS    AC=12;AF=0.500;AN=24;BaseQRankSum=1.33;DP=259;ExcessHet=0.5921;FS=4.427;InbreedingCoeff=0.3333;MLEAC=12;MLEAF=0.5;MQ=60;MQRankSum=0;QD=27.19;ReadPosRankSum=0.232;SOR=0.391     GT:AD:DP:GQ:PL  0/0:29,0:29:81:0,81,1215        0/1:6,9:15:99:297,0,205      0/1:7,4:11:99:102,0,259 1/1:0,28:28:84:1089,84,0        1/1:0,24:24:72:934,72,0 0/0:25,0:25:69:0,69,1035        1/1:0,18:18:54:617,54,0 0/1:13,10:23:99:344,0,440       0/0:22,0:22:66:0,66,815 0/0:21,0:21:57:0,57,855 1/1:0,21:21:63:807,63,0 0/1:14,5:19:99:152,0,450

bcftools view -H SUPER_1_GATK_NOMISSING_NOINDELS.vcf.gz | head
SUPER_1 46182   .       C       .       .       PASS    DP=260  GT:AD:DP:RGQ    0/0:25:25:66    0/0:39:39:99    0/0:17:17:45    0/0:19:19:57    0/0:10:10:30    0/0:27:27:81    0/0:17:17:39    0/0:24:24:63    0/0:14:14:42    0/0:14:14:39    0/0:36:36:99    0/0:18:18:45
SUPER_1 46183   .       T       .       .       PASS    DP=262  GT:AD:DP:RGQ    0/0:24:24:66    0/0:41:41:99    0/0:17:17:45    0/0:19:19:57    0/0:10:10:30    0/0:29:29:87    0/0:18:18:42    0/0:22:22:63    0/0:14:14:42    0/0:14:14:39    0/0:36:36:99    0/0:18:18:45
SUPER_1 46184   .       A       .       .       PASS    DP=264  GT:AD:DP:RGQ    0/0:24:24:66    0/0:41:41:99    0/0:18:18:48    0/0:19:19:57    0/0:10:10:30    0/0:29:29:87    0/0:18:19:42    0/0:22:22:63    0/0:14:14:42    0/0:14:14:39    0/0:36:36:99    0/0:18:18:45
SUPER_1 46185   .       T       .       .       PASS    DP=268  GT:AD:DP:RGQ    0/0:24:24:66    0/0:42:42:99    0/0:18:18:45    0/0:19:19:57    0/0:10:10:30    0/0:31:31:90    0/0:19:19:42    0/0:23:23:57    0/0:14:14:42    0/0:14:14:39    0/0:36:36:99    0/0:18:18:45
SUPER_1 46186   .       C       .       .       PASS    DP=265  GT:AD:DP:RGQ    0/0:23:24:57    0/0:42:42:99    0/0:18:18:45    0/0:19:19:57    0/0:10:10:30    0/0:31:31:90    0/0:18:18:42    0/0:23:23:57    0/0:14:14:42    0/0:14:14:39    0/0:36:36:99    0/0:16:16:45
SUPER_1 46187   .       T       .       .       PASS    DP=265  GT:AD:DP:RGQ    0/0:24:24:66    0/0:41:41:99    0/0:17:17:45    0/0:19:19:57    0/0:10:10:30    0/0:31:31:90    0/0:18:18:39    0/0:23:23:57    0/0:14:14:42    0/0:14:14:39    0/0:37:37:99    0/0:17:17:48
SUPER_1 46188   .       C       .       .       PASS    DP=265  GT:AD:DP:RGQ    0/0:24:24:66    0/0:41:41:99    0/0:17:17:45    0/0:20:20:60    0/0:10:10:30    0/0:31:31:90    0/0:16:16:39    0/0:23:23:57    0/0:14:14:42    0/0:14:14:39    0/0:38:38:99    0/0:17:17:48
SUPER_1 46189   .       T       .       .       PASS    DP=269  GT:AD:DP:RGQ    0/0:24:24:66    0/0:41:41:99    0/0:17:17:45    0/0:21:21:63    0/0:10:10:30    0/0:33:33:96    0/0:16:16:39    0/0:23:23:57    0/0:15:15:45    0/0:14:14:39    0/0:38:38:99    0/0:17:17:48
SUPER_1 46190   .       G       .       .       PASS    DP=269  GT:AD:DP:RGQ    0/0:22:22:66    0/0:41:41:99    0/0:17:17:45    0/0:22:22:66    0/0:10:10:30    0/0:33:33:96    0/0:16:16:39    0/0:23:23:60    0/0:15:15:45    0/0:14:14:39    0/0:39:39:99    0/0:17:17:48
SUPER_1 46191   .       A       .       .       PASS    DP=269  GT:AD:DP:RGQ    0/0:22:22:66    0/0:42:42:99    0/0:17:17:45    0/0:22:22:66    0/0:10:10:30    0/0:32:33:57    0/0:16:16:39    0/0:23:23:60    0/0:15:15:39    0/0:14:14:39    0/0:39:39:99    0/0:16:16:48

bcftools view -H SUPER_1_GATK_NOREF_NOINDELS.vcf.gz | head
[W::vcf_parse_info] INFO 'AF' is not defined in the header, assuming Type=String
SUPER_1 46282   .       C       A       4367.99 PASS    AC=11;AF=0.500;AN=22;BaseQRankSum=-0.015;DP=238;ExcessHet=0.2176;FS=0;InbreedingCoeff=0.4955;MLEAC=13;MLEAF=0.542;MQ=60;MQRankSum=0;QD=31.2;ReadPosRankSum=-0.144;SOR=0.732     GT:AD:DP:FT:GQ:PL       1/1:0,17:17:PASS:51:676,51,0    0/1:4,18:22:PASS:97:617,0,97 0/1:8,5:13:PASS:99:136,0,237    0/0:22,0:22:PASS:57:0,57,855    0/0:22,0:22:PASS:63:0,63,945    1/1:0,26:26:PASS:78:1016,78,0   0/0:11,0:11:PASS:27:0,27,405    0/1:11,10:21:PASS:99:355,0,335  1/1:0,16:16:PASS:48:600,48,0    ./.:0,8:8:DPFILTER:24:322,24,0  0/0:35,0:35:PASS:90:0,90,13501/1:0,17:17:PASS:51:664,51,0
SUPER_1 46521   .       T       C       1408.96 PASS    AC=6;AF=0.273;AN=22;BaseQRankSum=1.17;DP=197;ExcessHet=1.7087;FS=2.991;InbreedingCoeff=0.1104;MLEAC=6;MLEAF=0.25;MQ=60;MQRankSum=0;QD=17.83;ReadPosRankSum=0.283;SOR=1.075      GT:AD:DP:FT:GQ:PL       0/0:15,0:15:PASS:36:0,36,540    0/1:10,4:14:PASS:99:125,0,298        ./.:9,0:9:DPFILTER:24:0,24,360  0/1:8,4:12:PASS:99:131,0,264    1/1:0,12:12:PASS:36:500,36,0    0/0:20,0:20:PASS:60:0,60,733    0/1:11,6:17:PASS:99:195,0,384   0/1:9,15:24:PASS:99:486,0,301   0/0:23,0:23:PASS:60:0,60,900    0/0:11,0:11:PASS:33:0,33,447    0/0:19,0:19:PASS:54:0,54,810 0/0:21,0:21:PASS:57:0,57,855
SUPER_1 55581   .       A       T       513.81  PASS    AC=1;AF=0.250;AN=4;BaseQRankSum=0;DP=87;ExcessHet=3.7566;FS=0;InbreedingCoeff=0.0493;MLEAC=3;MLEAF=0.15;MQ=51.7;MQRankSum=1.07;QD=16.06;ReadPosRankSum=-0.862;SOR=0.941 GT:AD:DP:FT:GQ:PL       ./.:9,0:9:DPFILTER:27:0,27,240  ./.:5,2:7:DPFILTER:63:63,0,122       ./.:8,0:8:DPFILTER:21:0,21,315  ./.:0,2:2:PASS:.:0,0,0  ./.:0,1:1:PASS:.:0,0,0  ./.:6,0:6:DPFILTER:18:0,18,239  0/1:6,10:16:PASS:99:219,0,165   ./.:2,7:9:DPFILTER:56:249,0,56  ./.:4,0:4:DPFILTER:12:0,12,133  ./.:9,0:9:DPFILTER:18:0,18,270  0/0:12,0:12:PASS:36:0,36,387    ./.:4,0:4:DPFILTER:9:0,9,135
SUPER_1 55916   .       G       A       311.54  PASS    AC=1;AF=0.071;AN=14;BaseQRankSum=-3.044;DP=138;ExcessHet=0.3218;FS=0;InbreedingCoeff=0.568;MLEAC=3;MLEAF=0.136;MQ=57.26;MQRankSum=0.674;QD=15.58;ReadPosRankSum=-0.746;SOR=0.892        GT:AD:DP:FT:GQ:PL       0/0:22,0:22:PASS:66:0,66,749    ./.:4,0:4:DPFILTER:12:0,12,155       0/0:15,0:15:PASS:45:0,45,558    0/0:13,0:13:PASS:36:0,36,540    ./.:0,6:6:DPFILTER:18:229,18,0  0/0:16,0:16:PASS:45:0,45,675    0/1:7,7:14:PASS:99:108,0,240    ./.:9,0:9:DPFILTER:27:0,27,318  ./.:4,0:4:DPFILTER:12:0,12,167  ./.:0,5:5:PASS:.:0,0,0  0/0:15,0:15:PASS:45:0,45,524 0/0:14,0:14:PASS:33:0,33,495
SUPER_1 56253   .       A       C       1503.09 PASS    AC=5;AF=0.250;AN=20;BaseQRankSum=1.66;DP=162;ExcessHet=0.0824;FS=1.413;InbreedingCoeff=0.6559;MLEAC=5;MLEAF=0.227;MQ=54.83;MQRankSum=0;QD=32.68;ReadPosRankSum=-0.354;SOR=0.557 GT:AD:DP:FT:GQ:PL       1/1:0,15:15:PASS:45:675,45,0    0/0:12,0:12:PASS:30:0,30,450 ./.:4,7:11:PASS:.:0,0,0 0/0:16,0:16:PASS:48:0,48,588    ./.:6,1:7:DPFILTER:5:0,5,167    1/1:0,17:17:PASS:51:721,51,0    0/0:12,0:12:PASS:33:0,33,495    0/0:17,0:17:PASS:51:0,51,596    0/0:15,1:16:PASS:33:0,33,553    0/0:10,0:10:PASS:30:0,30,357    0/0:10,1:11:PASS:18:0,18,356    0/1:10,4:14:PASS:99:138,0,408
SUPER_1 56344   .       T       A       3438.21 PASS    AC=12;AF=0.750;AN=16;BaseQRankSum=0.66;DP=150;ExcessHet=1.7087;FS=0;InbreedingCoeff=0.1081;MLEAC=18;MLEAF=0.75;MQ=54.34;MQRankSum=0;QD=24.91;ReadPosRankSum=0.244;SOR=0.659     GT:AD:DP:FT:GQ:PL       1/1:0,14:14:PASS:42:511,42,0    0/1:8,4:12:PASS:99:119,0,252 ./.:0,9:9:DPFILTER:27:296,27,0  1/1:0,12:12:PASS:36:401,36,0    ./.:0,6:6:DPFILTER:18:250,18,0  1/1:0,20:20:PASS:60:743,60,0    ./.:0,8:8:DPFILTER:24:270,24,0  0/1:9,7:16:PASS:99:232,0,192    0/1:11,2:13:PASS:25:25,0,453    ./.:8,0:8:DPFILTER:24:0,24,300  1/1:0,11:11:PASS:33:342,33,00/1:7,10:17:PASS:99:228,0,157
SUPER_1 65437   .       G       A       478.87  PASS    AC=1;AF=0.042;AN=24;BaseQRankSum=-2.027;DP=235;ExcessHet=3.0103;FS=1.603;InbreedingCoeff=-0.0437;MLEAC=1;MLEAF=0.042;MQ=60;MQRankSum=0;QD=19.15;ReadPosRankSum=-1.728;SOR=0.412 GT:AD:DP:GQ:PL  0/0:21,0:21:60:0,60,900 0/0:12,1:13:24:0,24,418 0/0:14,0:14:36:0,36,540      0/0:27,0:27:72:0,72,1080        0/0:21,0:21:60:0,60,900 0/0:19,0:19:54:0,54,810 0/0:20,0:20:54:0,54,810 0/0:23,0:23:66:0,66,990 0/0:19,0:19:57:0,57,760 0/0:15,0:15:45:0,45,592 0/0:18,0:18:51:0,51,765 0/1:9,16:25:99:491,0,301
SUPER_1 65899   .       A       C       36.57   PASS    AC=1;AF=0.056;AN=18;BaseQRankSum=0;DP=152;ExcessHet=3.0103;FS=3.256;InbreedingCoeff=-0.077;MLEAC=1;MLEAF=0.042;MQ=60;MQRankSum=0;QD=2.61;ReadPosRankSum=-1.184;SOR=2.303        GT:AD:DP:FT:GQ:PL       0/0:14,1:15:PASS:15:0,15,471    0/0:10,0:10:PASS:27:0,27,405 0/0:13,2:15:PASS:0:0,0,397      ./.:7,0:7:DPFILTER:21:0,21,280  0/0:14,0:14:PASS:39:0,39,585    0/0:19,1:20:PASS:44:0,44,612    0/0:11,2:13:PASS:8:0,8,361      0/1:12,2:14:PASS:48:48,0,463    0/0:16,0:16:PASS:42:0,42,630    ./.:9,0:9:DPFILTER:18:0,18,270  0/0:10,1:11:PASS:17:0,17,338./.:8,0:8:DPFILTER:24:0,24,296
SUPER_1 65904   .       A       C       37.02   PASS    AC=1;AF=0.056;AN=18;BaseQRankSum=-1.857;DP=154;ExcessHet=3.0103;FS=3.256;InbreedingCoeff=-0.0927;MLEAC=1;MLEAF=0.042;MQ=60;MQRankSum=0;QD=2.64;ReadPosRankSum=-0.494;SOR=2.303  GT:AD:DP:FT:GQ:PL       0/0:18,0:18:PASS:51:0,51,765    ./.:9,0:9:DPFILTER:24:0,24,360       0/0:14,0:14:PASS:36:0,36,540    ./.:7,0:7:DPFILTER:18:0,18,270  0/0:15,0:15:PASS:39:0,39,585    0/0:16,0:16:PASS:48:0,48,574    0/0:11,0:11:PASS:33:0,33,407    0/1:12,2:14:PASS:48:48,0,463    0/0:17,0:17:PASS:45:0,45,675    ./.:9,0:9:DPFILTER:18:0,18,270  0/0:12,2:14:PASS:0:0,0,379   0/0:9,1:10:PASS:0:0,0,346
SUPER_1 82574   .       T       G       388.87  PASS    AC=1;AF=0.045;AN=22;BaseQRankSum=1.24;DP=258;ExcessHet=3.0103;FS=1.603;InbreedingCoeff=-0.0437;MLEAC=1;MLEAF=0.042;MQ=60;MQRankSum=0;QD=17.68;ReadPosRankSum=1.42;SOR=0.412     GT:AD:DP:FT:GQ:PL       0/0:24,0:24:PASS:66:0,66,990    0/0:17,0:17:PASS:51:0,51,608 0/0:23,0:23:PASS:69:0,69,754    0/0:25,0:25:PASS:63:0,63,945    0/0:20,0:20:PASS:57:0,57,855    0/0:31,0:31:PASS:90:0,90,1350   ./.:9,0:9:DPFILTER:24:0,24,360  0/0:24,0:24:PASS:54:0,54,810    0/0:22,0:22:PASS:60:0,60,900    0/0:23,0:23:PASS:54:0,54,810    0/0:17,0:17:PASS:48:0,48,7200/1:10,12:22:PASS:99:401,0,278
'''