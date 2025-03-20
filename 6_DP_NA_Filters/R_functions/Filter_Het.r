#' Filter Individuals and Positions by Alternative Rate
#'
#' This function filters a VCF file based on the alternative allele rate both by individuals and by positions. 
#' The function calculates the percentage of heterozygous (alternative) genotypes for each individual and for each position,
#' and filters the VCF file according to the provided thresholds. The function also handles ploidy (1 or 2) by adjusting
#' the genotypes accordingly ("0/1" for ploidy 2, "1" for ploidy 1).
#'
#' @param vfc A VCF object (vcfR format) to be filtered.
#' @param ploidy (integer) The ploidy level of the dataset. Either 1 or 2.
#' @param rate_het_max_POS (double) The maximum allowed alternative allele rate (percent) for positions.
#' @param rate_het_max_Ind (double) The maximum allowed alternative allele rate (percent) for individuals.
#' @param rate_het_min_pos (integer) The minimum number of alternative alleles allowed per position.
#'
#' @return A filtered VCF object where positions and individuals are removed based on the alternative allele rate filters.
#' @import vcfR
#' @export
#'
#' @examples
#' # Filter a VCF file based on alternative rate thresholds for ploidy 2
#' filtered_vcf <- Filters_Het_Ind_Pos(vcf_file, 2, 30.0, 50.0, 0.0)
#'
#' # Filter a VCF file based on alternative rate thresholds for ploidy 1
#' filtered_vcf <- Filters_Het_Ind_Pos(vcf_file, 1,  30.0, 50.0, 0.0)

Filters_Het_Ind_Pos <- function(vfc, ploidy, rate_het_max_POS, rate_het_max_Ind, rate_het_min_pos) {
  
  # Extract genotypes from the VCF based on the ploidy level
  if (ploidy == 1) {
    gt <- extract.gt(vfc, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, convertNA = FALSE, extract = TRUE)
  } else if (ploidy == 2) {
    gt <- extract.gt(vfc, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, convertNA = FALSE, extract = TRUE)
  } else {
    stop("Ploidy must be either 1 or 2")
  }

  # Calculate alternative allele rate for individuals (columns)
  if (ploidy == 1) {
    somme_het_ind <- colSums(gt == "1")
  } else {
    somme_het_ind <- colSums(gt == "0/1")
  }
  
  rate_het_ind <- (somme_het_ind / dim(gt)[1]) * 100
  samples_kept <- names(rate_het_ind[rate_het_ind <= rate_het_max_Ind])
  
  # Filter VCF by individuals
  vcf_filtered_ind <- vfc[, c("FORMAT", samples_kept)]

  # Extract genotypes for filtered VCF
  gt_1 <- extract.gt(vcf_filtered_ind, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, convertNA = FALSE, extract = TRUE)

  # Calculate alternative allele count for positions (rows)
  if (ploidy == 1) {
    somme_het_pos <- rowSums(gt_1 == "1")
  } else {
    somme_het_pos <- rowSums(gt_1 == "0/1")
  }
  
  pos_kept <- which(somme_het_pos >= rate_het_min_pos)
  
  # Filter VCF by positions with enough alternative alleles
  vcf_filtered_het_min <- vcf_filtered_ind[pos_kept, ]
  
  # Calculate alternative allele rate for positions (rows)
  gt_1 <- extract.gt(vcf_filtered_het_min, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, convertNA = FALSE, extract = TRUE)
  
  if (ploidy == 1) {
    somme_het_pos <- rowSums(gt_1 == "1")
  } else {
    somme_het_pos <- rowSums(gt_1 == "0/1")
  }
  
  rate_het_pos <- (somme_het_pos / dim(gt_1)[2]) * 100
  pos_kept <- which(rate_het_pos <= rate_het_max_POS)
  
  # Final filtered VCF by both position and individual heterozygosity rates
  vcf_filtered_het <- vcf_filtered_het_min[pos_kept, ]
  
  # Return filtered VCF
  return(vcf_filtered_het)
}
