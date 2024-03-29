# ##############################################################
# #### Filter Individues and Position by alternative rate #####
# ##############################################################

# Library
#----------#
library("vcfR")

# input
# ----------#
# A VCF file with Ploidy 1

# arg
# -------#
# rate_het_max_POS (double) : Max percent of alternative site in one position
# rate_het_max_Ind (double) : Max percent of alternative site in one individus
# rate_het_min_pos (int) : Min nb of alternative site in one position

# output
# -------#
# VCF file filtered for alternative rate by Individu and then by position

# Usage :
# --------#
# Filter_VCF=Filters_Het_Ind_Pos(vcf_file, rate_het_max_POS,rate_het_max_Ind)

# Function 
# -------------------------------------------------------------------------------------------------------------------------------------------#
# Get genotype from vfc with the vcfR package
# Count number of "1" per column and compute alternative rate (percent) by divided by the number of line (position)
# Create new vcf with only samples which pass the filter
# On the new VCF : Count number of "1" per line and compute alternative rate (percent) by divided by the number of column (individus)
# Create new vcf with only samples AND position which pass the filter
# -------------------------------------------------------------------------------------------------------------------------------------------#

Filters_Het_Ind_Pos<-function(vfc, rate_het_max_POS,rate_het_max_Ind, rate_het_min_pos){
	# Filters_Het_Ind_Pos(vcf_file, rate_het_max_POS,rate_het_max_Ind)
	# Get genotypes data frame
	gt<-extract.gt(vfc,element="GT",mask=F,as.numeric=F,return.alleles = F, convertNA = F,extract = T)
	# get ind rate het (column)
	somme_het_ind = colSums(gt=="1")
	rate_het_ind = (somme_het_ind/dim(gt)[1])*100
	samples_kept<-names(which(rate_het_ind<=rate_het_max_Ind))
	vcf_filtered_ind<-vfc[,c("FORMAT",samples_kept)]

	# get pos nb het min (lines)
	gt_1<-extract.gt(vcf_filtered_ind,element="GT",mask=F,as.numeric=F,return.alleles = F, convertNA = F,extract = T)
	somme_het_Pos = rowSums(gt_1=="1")
	Pos_kept<-which(somme_het_Pos>=rate_het_min_pos)
	vcf_filtered_Het_min<-vcf_filtered_ind[c(Pos_kept),]
	
	
	# get pos rate het (lines)
	gt_1<-extract.gt(vcf_filtered_Het_min,element="GT",mask=F,as.numeric=F,return.alleles = F, convertNA = F,extract = T)
	somme_het_Pos = rowSums(gt_1=="1")
	rate_het_Pos = (somme_het_Pos/dim(gt_1)[2])*100
	Pos_kept<-which(rate_het_Pos<=rate_het_max_POS)
	vcf_filtered_Het<-vcf_filtered_Het_min[c(Pos_kept),]
	
	return(vcf_filtered_Het)
}
