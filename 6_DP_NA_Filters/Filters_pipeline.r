#!/usr/bin/env Rscript

#==========================================================#
# Filters VCF files
#  Authors :  Elise Gay 
#  July 2022
#==========================================================#
#===============#
#===============#
# Load Libraries 
#===============#
#===============#
library(ggplot2)
library(reshape2)
library(gridExtra)
library(stringr)
library(vcfR)
# Load functions to apply filters on genotype and depth
#------------------------------------------------------#
# for filters on sequencing depth 
source("R_functions/Filter_DP.R")
# for filters on genotypes (heterozygous rate ; missing data)
# Each function can be used to filters individuals or positions based on genotype rates
source("R_functions/Filter_Het.r")
source("R_functions/Filter_Na.R")

#===============================#
#===============================#
# Input VCF Files + Description #
#===============================#
#===============================#

# Load VCF
#----------#

# Example :
#-----------------------------------------------------------------------------------------------#
# VCF on Y chromosome. Variant calling called with ploidy = 1 
# VCF already filters for lowqual, indel and repeat. Made from GATK filters on cluster account
# 25 male rats samples (data from Romuald Laso-Jadart data set)
#-----------------------------------------------------------------------------------------------#

VCF_P1 =read.vcfR("data/PPC_M_NC_051357.1_DP100_MQ40_GATK_TAG_Flowqual_Noindels_Norepeat_SNP.vcf.gz")

# Read metadata 
#----------------#
metadata=read.table("metadata/metadata.txt", header=TRUE)

#-----------------------------------------------------#
#-----------------------------------------------------#
# 1) Sequencing depth rate visualization and filter
#-----------------------------------------------------#
#-----------------------------------------------------#

# Get and format DP table for all position and all individuals
#--------------------------------------------------------------#
# get the position sequencing depth
DP<- extract.gt(VCF_P1, element='DP', as.numeric = TRUE) 

# get mean per column (individu)
moyenne_cover_par_indiv = apply(DP,2,mean,na.rm=T)
moyenne_cover_par_Ind_table=as.data.frame(moyenne_cover_par_indiv)

# get mean per line (position)
moyenne_cover_par_pos= apply(DP,1,mean,na.rm=T)
moyenne_cover_par_pos_table=as.data.frame(moyenne_cover_par_pos)

# Get the position vector in numeric
str_list=unlist(strsplit(rownames(DP)[1], '_'))
chrname=paste0(paste(str_list[-length(str_list)], collapse = "_"), "_")
position_vector=as.numeric(str_remove(rownames(DP), chrname))

# Add vector of position in the sequencing depth table 
moyenne_cover_par_pos_table_pos=cbind(moyenne_cover_par_pos_table, position_vector)

# PLOT detailed DP for each ind or position
#-------------------------------------------#
# plot DP by positions
ggplot()+
  geom_point(data=moyenne_cover_par_pos_table_pos, 
             aes(y=moyenne_cover_par_pos_table_pos$moyenne_cover_par_pos,
                 x=moyenne_cover_par_pos_table_pos$position_vector),
             color="black", pch=20, size=1) +
  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y="mean depth", x = "PB")+
  ggtitle("Average depth sequencing of SNPs")+
  scale_x_continuous(breaks = seq(0, 18e+06, by = 1000000)) +
  # vizualise only range of sequencing depth detph
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100))

# boxplot DP by individuals
# format table 
DP_melt=melt(DP)

ggplot()+
  geom_boxplot(data=DP_melt, aes(x=DP_melt$Var2, 
                                 y=DP_melt$value,
                                 colour="purple"),
               color="black") +
  theme(axis.title.x=element_blank())+
  labs(y="SNP sequencing depth", x = "samples")+
  ggtitle("depth") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # vizualise only range of sequencing depth detph
  scale_y_continuous(breaks = seq(0, 260, by = 20), limits = c(0, 260))

# FILTER FOR DP 
vcf_P1_DP = filter_vcf_by_dp(vcf = VCF_P1, min_dp_pos = 10, max_dp_pos = 50, min_dp_ind = 10, max_dp_ind = 50)
vcf_P2_DP = filter_vcf_by_dp(vcf = vcf_P2, min_dp_pos = 10, max_dp_pos = 50, min_dp_ind = 10, max_dp_ind = 50)

#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
# 2) Vizualisation of chosen genotype frequency NA or het
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#

#--------------------------------#
# Plot count and proportion
#--------------------------------#

# Extract genotype
#-------------------#
geno_table<-extract.gt(vcf_P2_na,element="GT",mask=F,as.numeric=F,return.alleles = F,
                 convertNA = F,extract = T)

# 1. Raw Count of genotype 
#---------------------------#
# In rowsum and colsum the number of '.' '1' and '0'and compute genotype count 
Ploidy=1  #or 1
GT = "Het" #or Na

if(Ploidy == 1 && GT =="Na") { gt="." 
}else if(Ploidy == 1 && GT =="Het") { gt="1" 
}else if(Ploidy == 2 && GT =="Na") {gt="./."
}else if(Ploidy == 2 && GT =="Het") {gt="0/1"}

# By position
sum_pos_geno_table = as.data.frame(rowSums(geno_table==gt))
# By individuals
sum_ind_geno_table = as.data.frame(colSums(geno_table==gt))

# 2. Genotype Frequencies
#---------------------------#

# By position
Prop_Pos_geno=((sum_pos_geno_table/dim(geno_table)[2])*100)
Prop_Pos_geno_table=as.data.frame(Prop_Pos_geno)

# Get the position vector in numeric
str_list=unlist(strsplit(rownames(geno_table)[1], '_'))
chrname=paste0(paste(str_list[-length(str_list)], collapse = "_"), "_")
position_vector=as.numeric(str_remove(rownames(geno_table), chrname))

# Add vector of position in the sequencing depth table 
Prop_Pos_geno_table_pos=cbind(Prop_Pos_geno_table, position_vector)

# By individuals
Prop_Ind_geno=((sum_ind_geno_table/dim(geno_table)[1])*100)
Prop_Ind_geno_table=as.data.frame(Prop_Ind_geno)

#----------------------------#
# PLOT genotype frequencies
#----------------------------#

# By position
ggplot()+
  geom_point(data=Prop_Pos_geno_table_pos, 
             aes(y=Prop_Pos_geno_table_pos[,c(1)] , 
                 x=Prop_Pos_geno_table_pos$position_vector), 
             color="purple",
             size=1)+
  labs(x="pos", y = paste("prop", gt)) +
  ggtitle("Prop gt site by position") +
  scale_x_continuous(breaks = seq(0, position_vector[length(position_vector)], 
                                  by = 1000000)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# By indivduals
ggplot()+
  geom_point(data=Prop_Ind_geno_table, 
             aes(y=Prop_Ind_geno_table[,c(1)] , 
                 x=row.names(Prop_Ind_geno_table)),
             color="purple",
             size=2) +
  labs(x="Ind", y =  paste("prop", gt)) +
  ggtitle("Prop gt site by individuals") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100))

# FILTER FOR NA RATE
vcf_P1_na = Filters_na_Ind_Pos(vcf = VCF_P1_DP, rate_na_max_POS = 50, rate_na_max_Ind = 10, ploidy = 1)

# FILTER for HET RATE 
vcf_P1_het = Filters_Het_Ind_Pos(vcf = vcf_P1_na, rate_het_max_POS = 80, rate_het_max_Ind = 10, ploidy = 1)
