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

Filters_Het_Ind_Pos <- function(vcf, ploidy, rate_het_max_POS, rate_het_max_Ind) {
  
  # Validate ploidy value
  if (!ploidy %in% c(1, 2)) {
    stop("Ploidy must be either 1 or 2.")
  }
  
  # Define genotype pattern based on ploidy
  het_pattern <- if (ploidy == 1) "1" else "0/1"
  
	# Extract genotypes from VCF
	gt <- extract.gt(VCF_P1, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, convertNA = FALSE)

	# Calculate alternative allele rate for individuals (columns)
	somme_het_ind <- colSums(gt == het_pattern)
	rate_het_ind <- (somme_het_ind / nrow(gt)) * 100
	samples_kept <- names(rate_het_ind[rate_het_ind <= rate_het_max_Ind])

	# Filter VCF by individuals
	vcf_filtered_ind <- VCF_P1[, c("FORMAT", samples_kept)]

	# Re-extract genotypes for the filtered VCF
	gt_1 <- extract.gt(vcf_filtered_ind, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, convertNA = FALSE)

	# Calculate alternative allele count for positions (rows)
	somme_het_pos <- rowSums(gt_1 == het_pattern)
	rate_het_pos <- (somme_het_pos / ncol(gt)) * 100
	# Apply minimum threshold for alternative alleles at positions
	pos_kept_min <- which(rate_het_pos >= rate_het_max_POS)
	vcf_filtered_het <- vcf_filtered_ind[pos_kept_min, ]

  # Return the filtered VCF
  return(vcf_filtered_het)
}
