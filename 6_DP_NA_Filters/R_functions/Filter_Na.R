#' Filters VCF based on missing data rate for individuals and positions
#'
#' This function filters a VCF object by the missing data rate for both individuals (samples)
#' and positions (variants). The user can specify the ploidy level (1 or 2), and the function
#' will compute the missing data rate based on either the `.` (for ploidy 1) or `./.` (for ploidy 2)
#' genotype representation.
#'
#' @param vcf A VCF object (vcfR object) to be filtered.
#' @param rate_na_max_POS Maximum allowed percentage of missing data per position.
#' @param rate_na_max_Ind Maximum allowed percentage of missing data per individual.
#' @param ploidy Ploidy level (1 or 2) to determine the missing data representation (`.` or `./.`).
#' 
#' @return A filtered VCF object.
#' @export

Filters_na_Ind_Pos <- function(vcf, rate_na_max_POS, rate_na_max_Ind, ploidy) {
  
  # Check for valid ploidy value
  if (!ploidy %in% c(1, 2)) {
    stop("Ploidy must be either 1 or 2.")
  }
  
  # Get genotype data from the VCF
  gt <- extract.gt(vcf, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, convertNA = FALSE)
  
  # Set the missing genotype string based on ploidy
  missing_str <- if (ploidy == 1) "." else "./."
  
  # Calculate missing data rate for individuals (columns)
  missing_individuals <- colSums(gt == missing_str)
  rate_na_ind <- (missing_individuals / nrow(gt)) * 100
  samples_kept <- names(rate_na_ind[rate_na_ind <= rate_na_max_Ind])
  
  # Filter the VCF by individuals
  vcf_filtered_ind <- vcf[, c("FORMAT", samples_kept)]
  
  # Re-extract genotype data for the filtered individuals
  gt_filtered <- extract.gt(vcf_filtered_ind, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, convertNA = FALSE)
  
  # Calculate missing data rate for positions (rows)
  missing_positions <- rowSums(gt_filtered == missing_str)
  rate_na_pos <- (missing_positions / ncol(gt_filtered)) * 100
  positions_kept <- which(rate_na_pos <= rate_na_max_POS)
  
  # Filter the VCF by positions
  vcf_filtered_na <- vcf_filtered_ind[positions_kept, ]
  
  # Return the filtered VCF
  return(vcf_filtered_na)
}
