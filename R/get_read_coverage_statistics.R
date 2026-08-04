#' Extract gene IDs from coverage matrix row names
#'
#' @param x A vector of the row names of the coverage matrix 
#'
#' @returns A character vector with the genes IDs only
#' @keywords internal
get_gene_ids <- function(x) {
  as.character(unlist(strsplit(x, ":", fixed = TRUE)))[[3]]
}

#' Extract chromosome names from coverage matrix row names
#'
#' @param x A vector of the row names of the coverage matrix 
#'
#' @returns A character vector with the chromosome names only
#' @keywords internal
get_chrom <- function(x) {
  as.character(unlist(strsplit(x, ":", fixed = TRUE)))[[1]]
}

#' Extract region's genomic coordinates from coverage matrix row names
#'
#' @param x A vector of the row names of the coverage matrix 
#'
#' @returns A character vector with the region's genomic coordinates only
#' @keywords internal
get_pos <- function(x) {
  as.character(unlist(strsplit(x, ":", fixed = TRUE)))[[2]]
}

#' Extract region's start position from a vector of genomic coordinates
#'
#' @param x A vector of the genomic coordinates 
#'
#' @returns A numeric vector with the region's start positions only
#' @keywords internal
get_start_pos = function(x) {
  as.character(unlist(strsplit(x, "-", fixed = TRUE)))[[1]]
}

#' Extract region's end position from a vector of genomic coordinates
#'
#' @param x A vector of the genomic coordinates 
#'
#' @returns A numeric vector with the region's end positions only
#' @keywords internal
get_end_pos = function(x) {
  as.character(unlist(strsplit(x, "-", fixed = TRUE)))[[2]]
}

#' Extract genes genomic coordinates from the row names of the read coverage 
#' matrix
#'
#' @inheritParams perform_coverage_qc coverage_matrix
#' @returns A list of two elements: the coverage data matrix where the row names
#'    are only composed of the gene IDs, and a data frame with the gene's
#'    genomic coordinates.
#' @keywords internal
get_genomic_coordinates <- function(coverage_matrix) {
  gene_ids <- as.character(lapply(rownames(coverage_matrix), get_gene_ids))
  chrom <- as.character(lapply(rownames(coverage_matrix), get_chrom))
  pos <- as.character(lapply(rownames(coverage_matrix),get_pos))
  start <- as.numeric(as.character(lapply(pos, get_start_pos)))
  end <- as.numeric(as.character(lapply(pos, get_end_pos)))
  genomic_coordinates <- data.frame(cbind(chrom, start, end))
  names(genomic_coordinates) <- c("chr", "start", "end")
  genomic_coordinates[["start"]] <- as.numeric(genomic_coordinates[["start"]])
  genomic_coordinates[["end"]] <- as.numeric(genomic_coordinates[["end"]])

  # update matrix row names
  rownames(coverage_matrix) <- gene_ids
  return(list(
    coverage_matrix = coverage_matrix,
    genomic_coordinates = genomic_coordinates
  ))
}

#' Get read coverage statistics per sample and per region
#' 
#' The function returns: the total number of good quality reads (accounted for
#' during variant calling process) per sample across all genes, the mean number
#' of good reads per sample across all genes, the gene length, the median number
#' of reads covered by each gene, and mean read length per sample.
#'
#' @param @param coverage_matrix A numeric matrix with the coverage data. This
#'    is the first element of the output from the `calculate_read_coverage()`
#'    function.
#' @param mean_read_length A data frame with the mean read length per sample.
#'    This is the first element of the output from the
#'    `calculate_read_coverage()`
#'    function.
#'
#' @returns A list of three data frames: one with the total, mean coverage, and
#'    the mean read length per sample. Another one with the region's genomic
#'    coordinates, lengths, and median number of read across all samples. The
#'    third one corresponds to the read coverage matrix.
#' @export
#'
#' @examples
#' \dontrun{
#'   cov_stats <- get_coverage_stats(
#'     coverage_matrix = per_gene_coverage,
#'     mean_read_length = sample_mean_read_length
#'   )
#'   coverage_matrix <- cov_stats[["coverage_matrix"]]
#'   samples_cov_stats <- cov_stats[["meta"]]
#'   regions_cov_stats <- cov_stats[["details"]]
#' }
get_coverage_stats <- function(coverage_matrix, mean_read_length) {
  checkmate::assert_matrix(coverage_matrix, mode = "numeric",
                           null.ok = FALSE)

  tmp <- get_genomic_coordinates(coverage_matrix)
  coverage_matrix <- tmp[["coverage_matrix"]]
  details <- tmp[["genomic_coordinates"]]

  # add region coverage stats
  details[["length"]] <- details[["end"]] - details[["start"]]
  details[["meadian_cov"]] <- as.numeric(apply(coverage_matrix, 1, median))
  details[["total_cov"]] <- as.numeric(apply(coverage_matrix, 1, sum))
  details[["mean_cov"]] <- as.numeric(apply(coverage_matrix, 1, mean))

  # add sample coverage stats
  meta <- mean_read_length
  meta[["total_cov"]] <- as.numeric(apply(coverage_matrix, 2, sum))
  meta[["mean_cov"]] <- as.numeric(apply(coverage_matrix, 2, mean))

  return(list(
    coverage_matrix <- coverage_matrix,
    meta <- meta,
    details <- details
  ))
}
