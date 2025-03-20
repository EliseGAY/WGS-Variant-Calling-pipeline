#' Filter Inds and Pos by DP Thresholds
#'
#' This function filters a VCF file based on the depth of coverage (DP) for both positions (pos) and individuals (ind). 
#' It removes positions and individuals whose DP is outside the specified minimum and maximum thresholds.
#'
#' @param vcf A VCF object (vcfR format) to be filtered.
#' @param min_dp_pos (double) The minimum allowed DP for positions.
#' @param max_dp_pos (double) The maximum allowed DP for positions.
#' @param min_dp_ind (double) The minimum allowed DP for individuals.
#' @param max_dp_ind (double) The maximum allowed DP for individuals.
#'
#' @return A filtered VCF object where positions and individuals are removed based on the DP thresholds.
#' @import vcfR
#' @export
#'
#' @examples
#' # Filter a VCF file based on DP thresholds for both pos and ind
#' filtered_vcf <- filter_vcf_by_dp(vcf_file, 10, 100, 20, 200)

filter_vcf_by_dp <- function(vcf, min_dp_pos, max_dp_pos, min_dp_ind, max_dp_ind) {
  
  # Extract the depth of coverage (DP) for all positions
  dp_data <- extract.gt(vcf, element = 'DP', as.numeric = TRUE)
  
  # Calculate the mean DP for each position (rows)
  mean_dp_pos <- rowMeans(dp_data, na.rm = TRUE)
  
  # Filter positions based on the DP thresholds
  filtered_vcf_pos <- vcf[mean_dp_pos >= min_dp_pos & 
                           mean_dp_pos <= max_dp_pos, ]
  
  # Extract the depth of coverage (DP) for ind after pos filtering
  dp_data_ind <- extract.gt(filtered_vcf_pos, element = 'DP', as.numeric = TRUE)
  
  # Calculate the mean DP for each ind (columns)
  mean_dp_ind <- colMeans(dp_data_ind, na.rm = TRUE)
  
  # Filter ind based on the DP thresholds
  valid_inds <- names(mean_dp_ind[mean_dp_ind >= min_dp_ind & 
                                  mean_dp_ind <= max_dp_ind])
  
  # Filter the VCF by ind
  filtered_vcf_ind_pos <- filtered_vcf_pos[, c("FORMAT", valid_inds)]
  
  return(filtered_vcf_ind_pos)
}
