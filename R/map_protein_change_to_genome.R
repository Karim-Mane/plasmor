#' Map Plasmodium falciparum protein change to a Single Nucleotide Polymorphism
#' (SNP) 
#'
#' @param ref_genome_file A character with the path to the reference genome FASTA
#'    file
#' @param gff_file A character with the path the to reference genome GFF file
#' @param gene_name A character with the target gene name or gene ID
#' @param protein_position An integer with the protein change position
#' @param ref_aa A character with the one-letter amino acid code for the
#'    reference amino acid
#' @param alt_aa A character with the one-letter amino acid code for the
#'    alternative amino acid
#'
#' @returns
#' @export
#'
#' @examples
#' \dontrun{
#'   snp <- map_protein_change_to_genome(
#'     ref_genome_file = "PlasmoDB-68_Pfalciparum3D7_Genome.fasta",
#'     gff_file = "PlasmoDB-68_Pfalciparum3D7.gff",
#'     gene_name = "PF3D7_0523000",
#'     protein_position = 86,
#'     ref_aa = "N",
#'     alt_aa = "Y"
#'   )
#' }
map_protein_change_to_genome <- function(ref_genome_file,
                                         gff_file,
                                         gene_name,
                                         protein_position,
                                         ref_aa,
                                         alt_aa) {
  checkmate::assert_file_exists(ref_genome_file, access = "r")
  checkmate::assert_file_exists(gff_file, access = "r")
  checkmate::assert_character(gene_name, any.missing = FALSE, len = 1,
                              null.ok = FALSE)
  checkmate::assert_integer(protein_position, any.missing = FALSE, len = 1,
                            null.ok = FALSE)
  checkmate::assert_character(ref_aa, any.missing = FALSE, len = 1, n.chars = 1,
                              null.ok = FALSE)
  checkmate::assert_character(alt_aa, any.missing = FALSE, len = 1, n.chars = 1,
                              null.ok = FALSE)

  # reading in the genome FASTA file
  cli::cli_progress_step("Loading genome FASTA (this may take a moment)...")
  genome <- Biostrings::readDNAStringSet(ref_genome_file)
  names(genome) <- sub(" .*", "", names(genome))

  # Parse GFF file
  cli::cli_progress_step("Parsing GFF file...")
  gff_df <- parse_gff(gff_file)

  # define some variables
  cds_cache <- list() # cache CDS by gene ID to avoid re-parsing for same gene
  results_rows <- list() # results collector

  cli::cli_progress_step("Performing variant mapping...")
  # resolve gene ID
  gene_id <- resolve_gene_id(gene_name, gene_lookup)

  # get CDS for the target gene (use cache)
  if (!gene_id %in% names(cds_cache)) {
    cds_cache[[gene_id]] <- get_cds(gff_df, gene_id)
  }
  cds_df <- cds_cache[[gene_id]]
  if (is.null(cds_df)) {
    cli::cli_abort(c(
      x = "No CDS found in GFF for {gene_id}",
      "!" = "Did you specify the correct reference GFF file ?"
    ))
  }

  # extract the strand information
  strand <- unique(cds_df[["strand"]])

  # get codon genomic positions
  codon_bases <- get_codon_positions(cds_df, protein_position)
  if (is.null(codon_bases)) {
    cli::cli_abort(c(
      x = "Could not map protein position to the gene coding sequences",
      "!" = "Did you specify the correct protein change ?"
    ))
  }

  # get the codon sequence based on its coordinates on the genome
  codon_seq <- get_ref_codon(genome, codon_bases, strand)
  ref_codon_aa <- as.character(Biostrings::translate(codon_seq))
  chrom <- codon_bases[[1]]["chrom"]
  cli::cli_alert_info(
    "This is a protein change on chromosome {.val {chrom}}, strand \\\
    {.val {strand}}"
  )

  # send a warning if identified amino acid (AA) is not the same as reference AA
  if (ref_codon_aa != ref_aa) {
    cli::cli_inform(c(
      "!" = "Reference amino acid (AA) mismatch",
      i = "Got {.val {ref_codon_aa}}, but expected AA is {.val {ref_aa}}"
    ))
  }

  # get the corresponding nucleotide of the detected codon
  nts <- character(3)
  for (k in seq_along(codon_bases)) {
    chrom <- codon_bases[[k]]["chrom"]
    pos <- as.integer(codon_bases[[k]]["pos"])
    nts[k] <- as.character(
      Biostrings::subseq(genome[[chrom]], start = pos, end = pos)
    )
  }
  codon_seq <- Biostrings::DNAString(paste(nts, collapse = ""))
  ref_codon_str <- as.character(codon_seq)
  ref_codon_aa  <- as.character(Biostrings::translate(codon_seq))
  codon_pos <- c(
    sprintf("%s:%s", codon_bases[[1]]["chrom"], codon_bases[[1]]["pos"]),
    sprintf("%s:%s", codon_bases[[2]]["chrom"], codon_bases[[2]]["pos"]),
    sprintf("%s:%s", codon_bases[[3]]["chrom"], codon_bases[[3]]["pos"])
  )
  cli::cli_inform(c(
    "*" = "The dectected reference codon {.val {ref_codon_str}} codes for \\\
    the amino acid {.val {ref_codon_aa}}",
    "*" = "Positions of the codon nucleotides on the reference genome are: \\\
    {.val {toString(codon_pos)}}"
    
  ))

  # scan every SNP at each codon position
  bases <- c("A", "T", "C", "G")
  cli::cli_alert_success(
    sprintf("All possible SNPs for %s%d%s", ref_aa, protein_position, alt_aa)
  )
  cat(sprintf("%-10s %-20s %-8s %-8s %-12s %-6s %s\n",
              "CodonPos", "GenomicPos", "RefNt", "AltNt", "MutCodon", "AltAA",
              ""))
  cat(strrep("-", 75), "\n")

  for (k in 1:3) {
    ref_nt <- substr(ref_codon_str, k, k)
    genomic_pos <- codon_bases[[k]]["pos"]
    chrom <- codon_bases[[k]]["chrom"]

    for (alt_nt in bases) {
      if (alt_nt == ref_nt) next
      
      mut_codon <- strsplit(ref_codon_str, "")[[1]]
      mut_codon[k] <- alt_nt
      mut_codon_str <- paste(mut_codon, collapse = "")
      mut_aa <- as.character(
        Biostrings::translate(Biostrings::DNAString(mut_codon_str))
      )
      match_flag <- if (mut_aa == alt_aa) "  *** MATCH ***" else ""

      cat(sprintf("%-10d %-20s %-8s %-8s %-12s %-6s %s\n",
                  k, genomic_pos, ref_nt, alt_nt,
                  mut_codon_str, mut_aa, match_flag))
    }
    cat("\n")
  }
}
