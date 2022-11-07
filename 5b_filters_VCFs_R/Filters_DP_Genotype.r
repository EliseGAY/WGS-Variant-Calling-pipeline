#!/usr/bin/env Rscript

#==========================================================#
# Compute Vizualisation and Filters on VCF files
# 
#
#  Authors : Stefano Mona / Romuald Laso-Jadart / Elise Gay 
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


# Stefano Mona libraries housemade
#----------------------------------#
source("function/libreria_filtri_VCF_e_SFS_unfolded.r")


# function to apply filters on genotype and depth
#--------------------------------------------------#
source("function/fonction_filtre_cover.R")
source("function/Filter_Ref.R")
source("function/Filter_Het.R")
source("function/Filter_Hom_ref.R")
source("function/Filter_Na.R")

# get/set dir
getwd()
setwd("~/../Desktop/capture_final/")

#===============================#
#===============================#
# Input VCF Files + Description #
#===============================#
#===============================#

# Clean VCF
#----------#
# Data with base filters on lowqual, indel and repeat. Made from GATK filters on cluster account
data=read.vcfR("data/DP_10_300_norepeat/VCF_214_capture_DP300_GATK_TAG_Flowqual_Noindels_No_repeat_SNP.vcf.gz")

#==============================================================#
#==============================================================#
# Visualization and fileters DP / MISSING DATA / Genotype RATE
#==============================================================#
#==============================================================#

#--------------------------------------------#
#--------------------------------------------#
# 1) Sequencing depth rate visualization
#--------------------------------------------#
#--------------------------------------------#

# Get and format DP table for all position and all individuals
#--------------------------------------------------------------#

# get the position sequencing depth
DP<- extract.gt(data, element='DP', as.numeric = TRUE) 

# get mean per column (individu)
moyenne_cover_par_indiv = apply(DP,2,mean,na.rm=T)
moyenne_cover_par_Ind_table=as.data.frame(moyenne_cover_par_indiv)

# summarise the data per individu
head(moyenne_cover_par_Ind_table)
summary(moyenne_cover_par_indiv)
boxplot(moyenne_cover_par_Ind_table, ylim = c(0, 300), ylab="DP_individu")

# get mean per line (position)
moyenne_cover_par_pos= apply(DP,1,mean,na.rm=T)
moyenne_cover_par_pos_table=as.data.frame(moyenne_cover_par_pos)
# summarise the data per position
summary(moyenne_cover_par_pos)
boxplot(moyenne_cover_par_pos_table, ylim = c(0, 400), ylab="DP_position")

# Get the position vector in numeric 
head(rownames(SNP_R1_R2_M_DP_mean_table))
# "SUPER_Y_207313" "SUPER_Y_235344" "SUPER_Y_399284" "SUPER_Y_400478" "SUPER_Y_401569" "SUPER_Y_403880"
position_vector=as.numeric(str_remove(rownames(SNP_R1_R2_M_DP_mean_table), "SUPER_Y"))

# Add vector of position in the sequencing depth table 
moyenne_cover_par_pos_table_pos=cbind(moyenne_cover_par_pos_table, position_vector)
head(SNP_R1_R2_M_DP_mean_table_pos)
dim(SNP_R1_R2_M_DP_mean_table_pos)


# PLOT detailed DP for each ind or position
#-------------------------------------------#

# plot DP by positions
ggplot()+
  geom_point(data=SNP_R1_R2_M_DP_mean_table_pos, aes(y=SNP_R1_R2_M_DP_mean_table_pos$SNP_R1_R2_M_DP_mean, 
                                                     x=SNP_R1_R2_M_DP_mean_table_pos$POS),
             color="blue") +
  
  theme(axis.title.x=element_blank())+
  labs(y="mean depth per position", x = "Nb position")+
  ggtitle("distribution SNPs (filtered, Na 20%) mean depth per position : R1+R2 in 11 Males ")+
  scale_y_continuous(breaks = seq(0, 50, by = 10), limits = c(0, 50)) +
  scale_x_continuous(breaks = seq(0, 5e+06, by = 500000))

# Plot DP Individu
moyenne_cover_par_Ind_table
ggplot()+
  geom_point(data=moyenne_cover_par_Ind_table, aes(y=moyenne_cover_par_Ind_table$moyenne_cover_par_indiv,
                                                   x=seq(1,dim(DP)[2])),
             color="blue") +
  labs(y="mean depth per individues", x = "Nb samples")+
  ggtitle("distribution position mean depth per ind")+
  scale_x_continuous(breaks = seq(0, dim(DP)[2], by=1), limits = c(0,dim(DP)[2]))+
  scale_y_continuous(breaks = seq(0, 300, by = 10), limits = c(0,300)) +
  geom_hline(yintercept = c(50,150), 
             color = "red", size=1)

#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
# 2) Vizualisation of any genotype frequency
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#

#--------------------------------#
# Plot count and proportion
#--------------------------------#
# Extract genotype
geno<-extract.gt(clean_vcf_214_DP10_200_50_150,element="GT",mask=F,as.numeric=F,return.alleles = F,
                   convertNA = F,extract = T)

# Replcae genotype selected by '1' and others by '0' in a new dataframe
# !!!!!! Put "1" in the gneotype line you want to count !!!!!!

for (i in 1:length(geno)) {
  if (geno[i]=="./.") {geno[i]<-1}
  if (geno[i]=="0/0") {geno[i]<-0}
  if (geno[i]=="0/1") {geno[i]<-0}
  if (geno[i]=="1/1") {geno[i]<-0} # if you want to count "1/1" genotype
  if (geno[i]=="0/2") {geno[i]<-0}
  if (geno[i]=="1/2") {geno[i]<-0}
  if (geno[i]=="2/2") {geno[i]<-0}
}
geno_table=as.data.frame(geno)

# Compute genotype count 
# position
sum_pos_geno_table = rowSums(geno_table==1)
summary(sum_pos_geno_table)

# individu
sum_ind_geno_table = colSums(geno_table==1)
summary(sum_ind_geno_table)

# Compute genotype proportion (either by divided by nb of pos : line , or nb of individus : column) and get dataframe
Prop_Pos_geno=((sum_pos_geno_table/dim(geno_table)[2])*100)
head(Prop_Pos_geno)
summary(Prop_Pos_geno)
Prop_Pos_geno_table=as.data.frame(Prop_Pos_geno)

Prop_Ind_geno=((sum_ind_geno_table/dim(geno_table)[1])*100)
head(Prop_Ind_geno)
summary(Prop_Ind_geno)
Prop_Ind_geno_table=as.data.frame(Prop_Ind_geno)

# PLOT
#-------------#
# Plot genotype distribution in a boxplot
boxplot(Prop_Pos_geno_table, ylim = c(0, 100), ylab="Proportion Na by position")
boxplot(Prop_Ind_geno_table, ylim = c(0, 100), ylab="Proportion Na by individu")

# plot either raw count or proportion of heterozygous
ggplot()+
  geom_point(data=Prop_Pos_geno_table, 
             aes(y=Prop_Pos_geno_table$Prop_Pos_geno , 
                 x=seq(1,dim(geno_table)[1])), 
             color="black",
             size=1) +
  labs(x="Ind", y = "Prop 0/1 site") +
  ggtitle("Prop 0/1 site per pos") +
  scale_x_continuous(breaks = seq(0, dim(geno_table)[1], by = 1000), 
                     limits = c(0,dim(geno_table)[1])) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100)) +
  #geom_vline(xintercept = c(3250,3750,7000,7250,8250,9000,12500,13000,14000,15000,15500,16000,19250,20750,22000,23000,23500,24000), 
  #           color = "red", size=1)+
  theme(axis.text=element_text(size=10))

# plot either raw count or proportion of heterozygous
ggplot()+
  geom_point(data=Prop_Ind_geno_table, 
             aes(y=Prop_Ind_geno_table$Prop_Ind_geno , 
                 x=seq(1,dim(geno_table)[2])), 
             color="black",
             size=1) +
  labs(x="Ind", y = "Prop ./. site") +
  ggtitle("Prop ./. site per ind") +
  scale_x_continuous(breaks = seq(0, dim(geno_table)[2], by = 1), limits = c(0,dim(geno_table)[2])) +
  #scale_y_continuous(breaks = seq(0, 5, by = 1), limits = c(0,5))
  geom_hline(yintercept = c(0,1), 
             color = "red", size=1)

#--------------------------------------------------#
# Combine all stats by individuals : DP, Het, Na
#--------------------------------------------------#

# get Na
Prop_Ind_Na_table=Prop_Ind_geno_table
colnames(Prop_Ind_Na_table)=c("Prop_Na")

# get Het
Prop_Ind_het_table=Prop_Ind_geno_table
colnames(Prop_Ind_het_table)=c("Prop_Het")
head(Prop_Ind_het_table)

# get DP
colnames(moyenne_cover_par_Ind_table)=c("Average_DP")

# compbine Na, Het and DP in one table
All_stat_VCF_213_ind=cbind(moyenne_cover_par_Ind_table,Prop_Ind_Na_table,Prop_Ind_het_table)
write.table(All_stat_VCF_213_ind, "data/DP_10_300_norepeat_NOHetRegion/VCF_214_capture_DP300_GATK_Flowqual_Noindels_No_repeat_SNP_NoHET_Region.STAT",
            sep = "\t",
            row.names = TRUE,
            col.names = TRUE)
			
#                     Average_DP  Prop_Na  Prop_Het
# GN1060_S4_L001     69.86732 1.632653 0.3363983
# GN1061_S5_L001     56.65768 2.112581 0.5786051
# GN10687_S7_L001    22.75972 6.400538 1.7896389


# Plot all stat :

# read_table:
table_stat=read.table("data/DP_10_300_norepeat_NOHetRegion/VCF_214_capture_DP300_GATK_Flowqual_Noindels_No_repeat_SNP_NoHET_Region.STAT", 
               row.names = 1, header = TRUE)
# select (if needed) the chosen ind :
WGS_12=c("GN1430","GN1434","GN1436","GN17525",
         "GN18404","GN18411","GN17582","GN17597",
         "GN17608","GN17768","GN17774", 
         "GN18167_S15_L001", "GN17539_S12_L001")

tabl_stat_12WGS=table_stat[WGS_12,]

# Plot table
data_plot=All_stat_VCF_146_ind
 
# x=row.names(data_plot)
# scale_x_continuous(breaks = seq(0, dim(data_plot)[1], by = 10), limits = c(0,dim(data_plot)[1])) 

# prop Na and Het
ggplot()+
  geom_point(data=data_plot, 
             aes(y=data_plot$Prop_Na, 
                 x=seq(1,dim(data_plot)[1])), 
             color="blue",
             size=1)+
  geom_point(data=data_plot, 
             aes(y=data_plot$Prop_Het, 
                 x=seq(1,dim(data_plot)[1])), 
             color="red",
             size=1) +
  geom_hline(yintercept = c(10,1), 
             color = c("blue","red"), size=1) +
  labs(x="Samples", y = "Percent") +
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle("Rate (%) of Heterozygous sites (red) and missing data sites (blue) per individu") +
  scale_x_continuous(breaks = seq(0, dim(data_plot)[1], by = 10), limits = c(0,dim(data_plot)[1])) +
  scale_y_continuous(breaks = seq(0, 30, by = 2), limits = c(0,30))

# DP 
ggplot()+
  geom_point(data=data_plot, 
             aes(y=data_plot$Average_DP , 
                 x=seq(1,dim(data_plot)[1])), 
             color="black",
             size=1) +
  labs(x="Ind", y = "Depth") +
  ggtitle("Average sequencing depth per samples") +
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_continuous(breaks = seq(0, dim(data_plot)[1], by = 10), limits = c(0,dim(data_plot)[1])) +
  scale_y_continuous(breaks = seq(0, 400, by = 50), limits = c(0,400))+
  geom_hline(yintercept = c(50,150), 
           color = "red", size=1)

#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
# 3) Apply filters 
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#

# a. Filter manually regions with low heterozygous rate per site
#-------------------------------------------------------------#

# Regions were selected from the heterozygous rate plot per position (see slide n°6)
data_REGION=data[c(0:3250, 3750:7000,7250:8250,9000:12500,13000:14000,15000:15500,16000:19250,
                   20750:22000,23000:23500,24000:28786),]
data_REGION

# b. Apply filters on sequencing Depth on ind and by positions
#--------------------------------------------------------------#
# ARGUMENTS : vcf,taux_min_cover_pos ,taux_max_cover_pos, taux_min_cover_ind, taux_max_cover_ind
clean_vcf_214_DP10_200_50_150<-fun_filtre_cover(data_REGION, 10.0, 200.0, 50.0, 150.0)
clean_vcf_214_DP10_200_50_150

#---------------------------#
# c. FILTERS BY INDIVIDUALS 
#---------------------------#

# Apply filters on missing data
#-------------------------------------#
# ARGUMENTS : vcf, taux max de NA par position, taux max per ind
# Method : 
# sum of "genotype" in vcf@gt has to be inferior or equal to the rate chosen in arg (which(rate_na_ind<=rate_na_max_Ind))
# The filters is first applied on individuals and then on position
clean_vcf_214_DP10_200_50_150_Na10Ind=Filters_na_Ind_Pos(clean_vcf_214_DP10_200_50_150, 100.0, 10.0) 
clean_vcf_214_DP10_200_50_150_Na10Ind


# Apply filters on Heterozygous site
#-------------------------------------#
# ARGUMENTS : 
# rate_het_max_POS (double) : Max percent of heterozygous site in one position
# rate_het_max_Ind (double) : Max percent of heterozygous site in one individus
# rate_het_min_pos (int) : Min nb of heterozygous site in one position

# Method : 
# sum of "genotype" in vcf@gt has to be inferior or equal to the rate chosen in arg (which(rate_het_ind<=rate_het_max_Ind))
# The filters is first applied on individuals and then on position

clean_vcf_214_DP10_200_50_150_Na10Ind_Het1<-Filters_Het_Ind_Pos(clean_vcf_214_DP10_200_50_150_Na10Ind,100.0,1.0,1.0) # discard ind with max 1% het sites and all positions with no het sites at all (might be created by Ind filtering)
clean_vcf_214_DP10_200_50_150_Na10Ind_Het1

# Apply filters on hom ref sites
#-------------------------------------#
# FUNCTION = Filters_hom_Ind_Pos<-function(vfc,rate_homRef_max_POS,rate_homRef_max_Ind)
# final rate have to be sctrictly inferior (not equal) to the chosen rate

clean_vcf_214_DP10_200_50_150_Na10Ind_Het1_NoREF= Filters_hom_Ind_Pos(clean_vcf_214_DP10_200_50_150_Na10Ind_Het1, 100.0,100.0)
clean_vcf_214_DP10_200_50_150_Na10Ind_Het1_NoREF

#--------------------------------#
# d. Apply filter on individuals
#--------------------------------#

# Apply filters on missing data
#-------------------------------------#
clean_vcf_214_DP10_200_50_150_Na10Ind_Het1_NoREF=Filters_na_Ind_Pos(clean_vcf_214_DP10_200_50_150, 0.0 "or 0.20", 10.0) 
clean_vcf_214_DP10_200_50_150_Na10Ind_20_Het1_NoREF
clean_vcf_214_DP10_200_50_150_Na10Ind_0_Het1_NoREF

# Apply filters on Heterozygous site
#-------------------------------------#

clean_vcf_214_DP10_200_50_150_Na10Ind_0_Het1_080_NoREF<-Filters_Het_Ind_Pos(clean_vcf_214_DP10_200_50_150_Na10Ind,80.0,1.0,1.0) 
clean_vcf_214_DP10_200_50_150_Na10Ind_0_Het1_080_NoREF
clean_vcf_214_DP10_200_50_150_Na10Ind_20_Het1_080_NoREF

# Saving VCF files :
#--------------------#

# Write VCF filtered
detach("package:pegas", unload = TRUE)

# 150 samples
write.vcf(clean_vcf_214_DP10_200_50_150_Na10Ind_Het1, 
          "data/DP_10_300_norepeat_NOHetRegion/150_REGION_DP10_200_50_150_NaInd10_HetInd1.vcf.gz")

write.vcf(clean_vcf_214_DP10_200_50_150_Na10Ind_Het1_NoREF, 
          "data/DP_10_300_norepeat_NOHetRegion/150_REGION_DP10_200_50_150_NaInd10_HetInd1_noRef.vcf.gz")


write.vcf(clean_vcf_214_DP10_200_50_150_Na10Ind_0_Het1_080_NoREF, 
          "data/DP_10_300_norepeat_NOHetRegion/150_REGION_DP10_200_50_150_NaInd10_HetInd1__noRef_Na0_Het80.vcf.gz")

write.vcf(clean_vcf_214_DP10_200_50_150_Na10Ind_20_Het1_NoREF, 
          "data/DP_10_300_norepeat_NOHetRegion/150_REGION_DP10_200_50_150_NaInd10_HetInd1__noRef_Na20.vcf.gz")


#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# Vizualisation : 
# Sequencing depth and genotype frequencies along the Y chromosome ()
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

######### INPUT DATA #############################################################
# read SNP vcf in male in R1 and R2
# SNP were filtered 
# - to discard all missing data (Na =0)
# - For Sequencing depth to be between 5 and 30 by position and by indiviudals
#--------------------------------------------------------------------------------#
# R1
SNP_R1_M_Na20=read.vcfR("raw_data/WGS_YP1_DP5_25_GATK_TAG_Flowqual_Noindels_Norepeat_Male_order_DP5_30_Na20_20.vcf.gz")
SNP_R1_M_Na20

POS_R1_2_Na20=read.vcfR("raw_data/WGS_YP1_DP5_25_GATK_TAG_Flowqual_Noindels_Norepeat_Male_order_DP5_30_Na20_20.vcf.gz")
POS_R1_2_Na20

head(Pos_12_Male)
dim(Pos_12_Male)

#R2
SNP_R2_M_Na20=read.vcfR("raw_data/WGS_YP1_DP5_30_GATK_TAG_Flowqual_Noindels_Norepeat_SNP_12Morder_R2_Na20.vcf.gz")
SNP_R2_M_Na20
# Because one individuals were discarded in R1 we have to remove it in R2 VCF
SNP_R2_M_Na20=SNP_R2_M_Na20[,colnames(SNP_R1_M_Na20@gt)]
SNP_R2_M_Na20

######### Compute mean DEPTH ###################

# extract sequencing depth table from for each R1 and R2 VCF
SNP_R2_11M_DP<- extract.gt(SNP_R2_M_Na20, element='DP', as.numeric = TRUE) 
SNP_R1_M_DP<- extract.gt(SNP_R1_M_Na20, element='DP', as.numeric = TRUE) 
POS_R1_2_Na20_DP<- extract.gt(POS_R1_2_Na20, element='DP', as.numeric = TRUE) 

# Merge in one unique table the SNP info in R1 and R2 in 11 Males
SNP_R1_R2_M_DP=rbind(SNP_R1_M_DP, SNP_R2_11M_DP)
dim(SNP_R1_R2_M_DP)
head(SNP_R1_R2_M_DP)

# get mean depth per line (position)
POS_R1_2_Na20_DP_mean= apply(POS_R1_2_Na20_DP,1,mean,na.rm=T)
POS_R1_2_Na20_DP_mean_table=as.data.frame(POS_R1_2_Na20_DP_mean)
head(POS_R1_2_Na20_DP_mean_table)
dim(POS_R1_2_Na20_DP_mean_table)

######### Get all position in a vector ###################
library(stringr)


# Get the position vector in numeric 
head(rownames(SNP_R1_R2_M_DP_mean_table))
# [1] "SUPER_Y_207313" "SUPER_Y_235344" "SUPER_Y_399284" "SUPER_Y_400478" "SUPER_Y_401569" "SUPER_Y_403880"
position_vector=as.numeric(str_remove(rownames(SNP_R1_R2_M_DP_mean_table), "SUPER_Y"))


######### Add vector of position in the sequencing depth table ###################

SNP_R1_R2_M_DP_mean_table_pos=cbind(SNP_R1_R2_M_DP_mean_table, position_vector)
head(SNP_R1_R2_M_DP_mean_table_pos)
dim(SNP_R1_R2_M_DP_mean_table_pos)
# [1] 685   2 = Na 0 
# [1] 2661    2 = Na 20

# plot Seuquencing depth of each SNPs relatively to their position in the Y chromosome

ggplot()+
  geom_point(data=SNP_R1_R2_M_DP_mean_table_pos, aes(y=SNP_R1_R2_M_DP_mean_table_pos$SNP_R1_R2_M_DP_mean, 
                                                     x=SNP_R1_R2_M_DP_mean_table_pos$POS),
             color="blue") +
  
  theme(axis.title.x=element_blank())+
  labs(y="mean depth per position", x = "Nb position")+
  ggtitle("distribution SNPs (filtered, Na 20%) mean depth per position : R1+R2 in 11 Males ")+
  scale_y_continuous(breaks = seq(0, 50, by = 10), limits = c(0, 50)) +
  scale_x_continuous(breaks = seq(0, 5e+06, by = 500000))


