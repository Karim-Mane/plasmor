#' Calculate the mean read length for all samples and get the read coverage
#' across each gene in the input BED file
#'
#' @param bam_files_dir A character with the full path to the input BAM files
#' @param bam_files_suffix A character with the BAM file suffix pattern
#' @param bed_file A character with the full path to the input BED file
#' @param gene_level_bed A logical used to determine whether the BED file only
#'    contains gene-level features or not (has several lines per gene ID,
#'    accounting for exon, CDS, etc genomic coordinates). If your BED file has
#'    annotation other than gene-level (contains exon, CDS coordinates), set the
#'    value of this parameter to `FALSE`.
#' @param count_dups A logical used to determine whether to account for
#'    duplicated reads or no when performing read counting. Set this to `FALSE`
#'    when determining coverage for variant calling from whole genome sequencing
#'    data.
#'
#' @returns A list of two elements: a matrix with the number of reads that
#'    overlap any given region in the BED file for all samples, and a data frame
#'    with the mean read length per sample.
#' @export 
#'
#' @examples
#' \dontrun{
#'   res <- calculate_read_coverage(
#'     bam_files_dir = "data/alignment_outputs",
#'     bam_files_suffix = "*_recal.bam",
#'     bed_file = "ref_genomes/PlasmoDB-68_Pfalciparum3D7_genes_only.bed",
#'     gene_level_bed = TRUE,
#'     count_dups = FALSE
#'   )
#'   per_gene_coverage <- res[["coverage_matrix"]]
#'   sample_mean_read_length <- res[["per_sample_mean_read_length"]]
#' }
calculate_read_coverage <- function(bam_files_dir,
                                    bam_files_suffix,
                                    bed_file,
                                    gene_level_bed = TRUE,
                                    count_dups = FALSE) {
  checkmate::assert_directory_exists(bam_files_dir, access = "r")
  checkmate::assert_character(bam_files_suffix, any.missing = FALSE,
                              null.ok = FALSE, len = 1)
  checkmate::assert_file_exists(bed_file)
  checkmate::assert_logical(gene_level_bed, len = 1, null.ok = FALSE)
  checkmate::assert_logical(count_dups, len = 1, null.ok = FALSE)

  # If input BED file contains all levels feature,
  # create gene-level only features BED file from input file
  if (!gene_level_bed) {
    cli::cli_inform(
      "Creating BED file with gene-level features only"
    )
    system(sprintf("awk '$8 ~ /_gene$/' %s > %s",
                   bed_file, file.path(dirname(bed_file), "genes_only.bed")))
    bed_file <- file.path(dirname(bed_file), "genes_only.bed")
  }
    
  # check whether the input directory contains BAM files with the specified
  # suffix
  bam_files <- list.files(
    path = bam_files_dir,
    pattern = bam_files_suffix,
    full.names = TRUE
  )
  if (!length(bam_files)) {
    cli::cli_abort(c(
      i = "No BAM file with suffix {bam_files_suffix} found in \\\
      {.dir {bam_files_dir}}",
      x = "Incorrect BAM file suffix provided"
    ))
  }
  sample_names <- gsub(bam_files_suffix, " ", basename(bam_files))
  
  # read in the BED file and define an IRanges object with the gene genomic
  # coordinates
  bed <- data.table::fread(bed_file, sep = "\t", nThread = 5, header = FALSE)
  genes_genomic_coordinates <- GenomicRanges::GRanges(
    seqnames = bed[["V1"]],
    ranges = IRanges::IRanges(start = bed[["V2"]], end = bed[["V3"]])
  )
  genes_genomic_coordinates <- sort(genes_genomic_coordinates)
  
  # define the output matrix
  coverage_matrix <- matrix(nrow = nrow(bed), ncol = length(bam_files))
  gene_id <- as.character(
    lapply(bed$V10, function(x) {unlist(strsplit(x, ";", fixed = TRUE))[[1]]})
  )
  gene_id <- as.character(
    lapply(gene_id, function(x) {unlist(strsplit(x, "=", fixed = TRUE))[[2]]})
  )
  rownames(coverage_matrix) <- paste0(
    GenomicRanges::seqnames(genes_genomic_coordinates), ":",
    GenomicRanges::start(genes_genomic_coordinates), "-",
    GenomicRanges::end(genes_genomic_coordinates), ":",
    gene_id
  )
  colnames(coverage_matrix) <- sample_names
  
  # define the data frame of mean read length per sample
  per_sample_mean_read_length <- data.frame(
    sample = sample_names,
    mean_read_length = NA,
    stringsAsFactors = FALSE
  )
    
  # loop trough the BAM files and calculate the mean read coverage
  cli::cli_progress_bar("Getting coverage...", total = length(bam_files))
  for (i in seq_len(length(bam_files))) {
    Sys.sleep(1 / length(bam_files))
    target_bam <- bam_files[i]
    sample <- sample_names[i]
    
    # get the mean read length for the sample
    what <- c("rname", "pos", "mapq", "qwidth")
    flag <- Rsamtools::scanBamFlag(
      isDuplicate = count_dups,
      isUnmappedQuery = FALSE,
      isNotPassingQualityControls = FALSE,
      isFirstMateRead = TRUE
    )
    param <- Rsamtools::ScanBamParam(what = what, flag = flag)
    bam <- Rsamtools::scanBam(target_bam, param = param)[[1]]
    mean_read_length <- round(mean(bam[["qwidth"]]))
    if (is.na(mean_read_length)) {
      flag <- Rsamtools::scanBamFlag(isDuplicate = count_dups,
                                     isUnmappedQuery = FALSE,
                                     isNotPassingQualityControls = FALSE)
      param <- Rsamtools::ScanBamParam(what = what, flag = flag)
      bam <- Rsamtools::scanBam(target_bam, param = param)[[1]]
      mean_read_length <- round(mean(bam[["qwidth"]], na.rm = TRUE))
    }
    per_sample_mean_read_length[["mean_read_length"]][i] <- mean_read_length

    # check if chromosome names in BED and BAM files are the same
    # if (!all(unique(bed[["V1"]]) %in% as.character(unique(bam[["rname"]])))) {
    #   cli::cli_abort(c(
    #     i = "Chromosome names in BED and BAM files are different",
    #     x = "Chromosome names in BED file must be the same as in BAM file."
    #   ))
    # }

    # get read coverage per gene region
    bam_ref <- GenomicRanges::GRanges(
      seqnames = bam[["rname"]],
      ranges = IRanges::IRanges(start = bam[["pos"]], width = bam[["qwidth"]])
    )
    # only consider reads with MAPQ >= 20
    bam_ref <- bam_ref[bam[["mapq"]] >= 20]
    # determine how many reads in BAM file overlap each region in BED file
    read_count <- GenomicRanges::countOverlaps(
      genes_genomic_coordinates,
      bam_ref
    )
    read_count <- as.matrix(read_count, ncol = 1)
    coverage_matrix[, i] <- read_count
    cli::cli_progress_update()
  }
  cli::cli_progress_done()
  
  return(list(
    coverage_matrix = coverage_matrix,
    per_sample_mean_read_length = per_sample_mean_read_length
  ))
}


#' Normalise coverage matrix
#' 
#' The function performs relative read count abundance normalization i.e. for
#' each sample, divides a region read count by the sum read count across all
#' regions (also known as total sum scaling). This gives the relative coverage
#' of each gene as a proportion of total coverage per sample.
#' Then computes the z-scores on the normalised coverage data.
#'
#' @param coverage_matrix A numeric matrix with the read count per region across
#'    all samples. This is the output of the `calculate_read_coverage()`
#'    function.
#' @returns A numeric matrix obtained from the normalization process explained
#'    in the function description.
#' @export
#' @examples
#' \dontrun{
#'   res <- calculate_read_coverage(
#'     bam_files_dir = "data/alignment_outputs",
#'     bam_files_suffix = "*_recal.bam",
#'     bed_file = "ref_genomes/PlasmoDB-68_Pfalciparum3D7_genes_only.bed",
#'     gene_level_bed = TRUE,
#'     count_dups = FALSE
#'   )
#'   per_gene_coverage <- res[["coverage_matrix"]]
#'   res <- normalize_coverage_matrix(per_gene_coverage)
#' }
normalize_coverage_matrix <- function(coverage_matrix){
  # calculate the relative read count abundance normalization i.e.
  # for each sample, divide a region read count by the sum read count across
  # all regions (also known as total sum scaling). This gives the relative
  # coverage of each gene as a proportion of total coverage per sample.
  normalised_coverage <- sweep(
    coverage_matrix, 2, colSums(coverage_matrix), FUN = "/"
  )
  
  # calculate the z-scores on the normalised coverage data
  # scale() applies Z-score standardization to each column of the matrix i.e.
  # z-score = (Xij - mean(X)) / sd(X)
  # Positive z-score: gene has higher coverage than average across samples
  # Negative z-score: gene has lower coverage than average across samples
  # Zero: gene coverage is exactly average
  scaled_coverage <- scale(normalised_coverage)
  scaled_coverage[is.nan(scaled_coverage)] <- 0
  
  return(scaled_coverage)
}


