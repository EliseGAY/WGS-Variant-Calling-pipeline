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
source("R_functions_ploidy1/libreria_filtri_VCF_e_SFS_unfolded.r")


# Load functions to apply filters on genotype and depth
#------------------------------------------------------#
# for filters on sequencing depth 
source("R_functions_ploidy1/fonction_filtre_cover.R")
# for filters on genotypes (heterozygous rate ; homozygous ; missing data)
# Each function can be used to filters individuals or positions based on genotype rates
source("R_functions_ploidy1/Filter_Het.R")
source("R_functions_ploidy1/Filter_Hom_ref.R")
source("R_functions_ploidy1/Filter_Na.R")

# get/set dir
getwd()
setwd("")

#===============================#
#===============================#
# Input VCF Files + Description #
#===============================#
#===============================#

# Load VCF
#----------#
# Data with base filters on lowqual, indel and repeat. Made from GATK filters on cluster account
R1_M=read.vcfR("raw_data/WGS_YP1_DP5_25_GATK_TAG_Flowqual_Noindels_Norepeat_SNP_12Morder_R1.vcf.gz")
R1_F=read.vcfR("raw_data/WGS_YP1_DP5_25_GATK_TAG_Flowqual_Noindels_Norepeat_SNP_9Forder_R1.vcf.gz")

# create vector with ID, in the same order in VCF file
samples_All=c("GN17775", "GN17768", "GN17774", "GN1430", "GN1434", "GN18972_S1_L004", "GN18973_S2_L004", "GN1436", "GN18983_S3_L004", "GN18986_S4_L004", "GN18411", "GN18417_S6_L004", "GN18404", "GN17526_S5_L004", "GN18410_S9_L004", "GN18424_S7_L004", "GN18426_S8_L004", "GN17525", "GN17582", "GN17597", "GN17608")
Male=c("GN17775", "GN17768", "GN1430", "GN1434", "GN18972_S1_L004", "GN18973_S2_L004", "GN18411", "GN18417_S6_L004", "GN18404", "GN17526_S5_L004", "GN17582", "GN17597")
Female=c("GN17774", "GN1436", "GN18983_S3_L004", "GN18986_S4_L004", "GN18410_S9_L004", "GN18424_S7_L004", "GN18426_S8_L004", "GN17525", "GN17608")

sampling_sites_all=c("WIO", "WIO", "WIO", "SPO", "SPO", "NPO", "NPO", "NPO", "NAO", "NAO", "NAO", "SPO", "SPO", "SPO", "SPO", "SPO", "SPO", "WIO", "WIO", "WIO", "WIO")
sampling_sites_Male=c("NAO", "NAO", "WIO", "WIO", "WIO", "WIO", "SPO", "SPO", "SPO", "SPO", "NPO", "NPO")
sampling_sites_Female=c("NAO", "WIO", "WIO", "WIO", "SPO", "SPO", "SPO", "SPO", "NPO")

sex=c("M", "M", "F", "F", "M", "M", "M", "F", "M", "F", "M", "M", "F", "M", "M", "F", "F", "M", "M", "F", "F")

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
DP<- extract.gt(R1_M, element='DP', as.numeric = TRUE) 
head(DP)
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
                                                   x=row.names(moyenne_cover_par_Ind_table)),
             color="blue") +
  labs(y="mean depth per individues", x = "Nb samples")+
  ggtitle("distribution position mean depth per ind : R2_8Male") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #scale_x_continuous(breaks = seq(0, dim(DP)[2], by=1), limits = c(0,dim(DP)[2]))
  scale_y_continuous(breaks = seq(0, 30, by = 2), limits = c(0,30))
  #geom_hline(yintercept = c(50,150), 
  #           color = "red", size=1)

#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
# 2) Vizualisation of any genotype frequency
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#

#--------------------------------#
# Plot count and proportion
#--------------------------------#

# Data : 
# As the Female were discarded from the VCF on Y, a run of dsicarding newly 
# homozygous reference site have been done before other filtering steps

# Extract genotype
head(R2_All@gt)
geno<-extract.gt(R2_8Male_DP5_30_Na20,element="GT",mask=F,as.numeric=F,return.alleles = F,
                   convertNA = F,extract = T)
head(geno)
#             GN1430 GN1434 GN17526_S5_L004 GN17582 GN17597 GN17768 GN17775 GN18404 GN18411 GN18417_S6_L004 GN18972_S1_L004 GN18973_S2_L004
# SUPER_Y_333 "."    "."    "."             "."     "."     "."     "."     "."     "."     "."             "."             "."            
# SUPER_Y_393 "0"    "0"    "."             "0"     "0"     "."     "."     "0"     "."     "."             "0"             "."            
# SUPER_Y_589 "."    "."    "."             "."     "1"     "."     "0"     "."     "0"     "0"             "0"             "0"            
# SUPER_Y_596 "."    "."    "."             "."     "0"     "."     "0"     "."     "0"     "0"             "0"             "1"            
# SUPER_Y_597 "."    "."    "."             "."     "0"     "."     "0"     "."     "0"     "0"             "0"             "1"            
# SUPER_Y_730 "."    "1"    "."             "."     "."     "."     "."     "."     "."     "."             "."             "."        

geno_table=as.data.frame(geno)

# Count in rowsum and colsum the number of '.' '1' and '0'and compute genotype count 
table(geno)

# position
sum_pos_geno_table = rowSums(geno==".") # put ".", "0" or "1" depending on the genotype you want to count
head(sum_pos_geno_table)

# individu
sum_ind_geno_table = colSums(geno==".")
sum_ind_geno_table # 
summary(sum_ind_geno_table)

# Compute genotype proportion (either by divided by nb of pos : line , or nb of individus : column) and get dataframe
Prop_Pos_geno=((sum_pos_geno_table/dim(geno_table)[2])*100)
head(Prop_Pos_geno)
summary(Prop_Pos_geno)
Prop_Pos_geno_table=as.data.frame(Prop_Pos_geno)

# Get the position vector in numeric 
head(rownames(Prop_Pos_geno_table))
# "SUPER_Y_207313" "SUPER_Y_235344" "SUPER_Y_399284" "SUPER_Y_400478" "SUPER_Y_401569" "SUPER_Y_403880"
position_vector=as.numeric(str_remove(rownames(Prop_Pos_geno_table), "SUPER_Y"))

# Add vector of position in the sequencing depth table 
Prop_Pos_geno_table_pos=cbind(Prop_Pos_geno_table, position_vector)
head(Prop_Pos_geno_table_pos)
dim(Prop_Pos_geno_table_pos)

# frequencies by individuals
Prop_Ind_geno=((sum_ind_geno_table/dim(geno_table)[1])*100)
head(Prop_Ind_geno)
summary(Prop_Ind_geno)
Prop_Ind_geno_table=as.data.frame(Prop_Ind_geno)
Prop_Ind_geno_table

# PLOT genotype
#---------------#

# plot either raw count or proportion of heterozygous
ggplot()+
  geom_point(data=Prop_Pos_geno_table_pos, 
             aes(y=Prop_Pos_geno_table_pos$Prop_Pos_geno , 
                 x=Prop_Pos_geno_table_pos$position_vector), 
             color="black",
             size=1)+
  labs(x="pos", y = "Prop '.' site") +
  ggtitle("Prop '.' site per pos : R2_8male") +
  scale_x_continuous(breaks = seq(0, dim(geno)[1], by = 100), 
                     limits = c(0,dim(geno)[1])) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100)) +
  theme(axis.text=element_text(size=10))

# plot either raw count or proportion of heterozygous
colnames(Prop_Ind_geno_table)
samples_All
Prop_Ind_geno_table
ggplot()+
  geom_point(data=Prop_Ind_geno_table, 
             aes(y=Prop_Ind_geno_table$Prop_Ind_geno , 
                 x=row.names(Prop_Ind_geno_table)),
             color="black",
             size=1) +
  labs(x="Ind", y = "Prop '.' site") +
  ggtitle("Prop '.' site per ind : R2_8male") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  #scale_x_continuous(breaks = seq(0, dim(geno_table)[2], by = 1), limits = c(0,dim(geno_table)[2]))
  #scale_y_continuous(breaks = seq(0, 5, by = 1), limits = c(0,5))
  #geom_hline(yintercept = c(0,1), 
  #           color = "red", size=1)

#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
# 3) Apply filters 
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#

#--------------------------------------------------------------#
# a. Apply filters on sequencing Depth on ind and/or by positions
#--------------------------------------------------------------#

# ARGUMENTS : vcf,taux_min_cover_pos ,taux_max_cover_pos, taux_min_cover_ind, taux_max_cover_ind

Males_R1_R2_DP5_30<-fun_filtre_cover(Males_R1_R2, 5, 30, 5, 30)
Males_R1_R2_DP5_30

#---------------------------#
# b. FILTERS BY INDIVIDUALS 
#---------------------------#

# Apply filters on missing data
#-------------------------------------#
# ARGUMENTS :function(vfc,rate_na_max_POS,rate_na_max_Ind)
# Method : 
# sum of "genotype" in vcf@gt has to be inferior or equal to the rate chosen in arg (which(rate_na_ind<=rate_na_max_Ind))
# The filters is first applied on individuals and then on position

Males_R1_R2_DP5_30_NaInd20=Filters_na_Ind_Pos(Males_R1_R2_DP5_30, 100.0, 20.0) # here filters only on individuals (max Na by positions is set to 100% = no filters)
R1_8Male_DP5_30_Na20

#----------------------------------------------------------------------#
# Apply filters on hom ref sites if individuals have been discarded 
#----------------------------------------------------------------------#
# FUNCTION = Filters_hom_Ind_Pos<-function(vfc,rate_homRef_max_POS,rate_homRef_max_Ind)
# final rate will be sctrictly inferior (not equal) to the chosen rate
# to discard all position remaining with 100 % homozygous genotype, set option to 100%

Males_R1_R2_DP5_30_NaInd20_NoREF= Filters_hom_Ind_Pos(Males_R1_R2_DP5_30_NaInd20, 100.0,100.0)
Males_R1_R2_DP5_30_NaInd20_NoREF

#---------------------------#
# FILTERS BY Positions 
#---------------------------#

# a. Apply filters on missing data
#-------------------------------------#
# ARGUMENTS :function(vfc,rate_na_max_POS,rate_na_max_Ind)
# Method : 
# sum of "genotype" in vcf@gt has to be inferior or equal to the rate chosen in arg (which(rate_na_ind<=rate_na_max_Ind))
# The filters is first applied on individuals and then on position

# Apply 20% max rate of missing data by positions
Males_R1_R2_DP5_30_NaInd20_NoREF_Na20=Filters_na_Ind_Pos(Males_R1_R2_DP5_30_NaInd20_NoREF, 20.0, 100.0) 
Males_R1_R2_DP5_30_NaInd20_NoREF_Na20

# Apply 0% max rate of missing data by positions
Males_R1_R2_DP5_30_NaInd20_NoREF_Na0=Filters_na_Ind_Pos(Males_R1_R2_DP5_30_NaInd20_NoREF, 0.0, 100.0) 
Males_R1_R2_DP5_30_NaInd20_NoREF_Na0

# Apply filters on Heterozygous site
#-------------------------------------#
# ARGUMENTS : 
# rate_het_max_POS (double) : Max percent of heterozygous site in one position
# rate_het_max_Ind (double) : Max percent of heterozygous site in one individus
# rate_het_min_pos (int) : Min nb of heterozygous site in one position

# Method : 
# sum of "genotype" in vcf@gt has to be inferior or equal to the rate chosen in arg (which(rate_het_ind<=rate_het_max_Ind))
# The filters is first applied on individuals and then on position

Males_R1_R2_DP5_30_NaInd20_NoREF_Na0_Het80<-Filters_Het_Ind_Pos(Males_R1_R2_DP5_30_NaInd20_NoREF_Na0,80.0,100.0,1.0)

Males_R1_R2_DP5_30_NaInd20_NoREF_Na0_Het80

# Saving VCF files :
#--------------------#

# Write VCF filtered
detach("package:pegas", unload = TRUE)
colnames(R1_12all_DP5_30_Na0@gt)
Males_R1_R2_DP5_30_Na20_20

write.vcf(Males_R1_R2_DP5_30_NaInd20_NoREF, 
          "raw_data/WGS_YP1_DP5_25_GATK_TAG_Flowqual_Noindels_Norepeat_Male_order_DP5_30_Na20_20_SNP.vcf.gz")


# Filtering done
# Follow analysis with Structure.R script ..