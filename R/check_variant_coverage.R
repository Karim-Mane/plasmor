#' Check whether a VCF file contains SNPs overlapping a set of genes of interest
#'
#' @param vcf_file A character with the path to the input `VCF` file
#' @param genes_file A character with the path to the `BED` file with the
#'    genomic coordinates of the genes of interest. It must contain at least the
#'    CHROM, START, and END columns
#' @param annotation_file A character with the path to file with the annotation
#'    columns for every variant
#'
#' @returns Writes out files with the variants overlapping each gene and the
#'    gene's variant counts.
#' @export
#'
#' @examples
#' \dontrun{
#'   get_variant_coverage(
#'     vcf_file = file.path("data", "variant_calling_outputs",
#'                          "biallelic_snps.vcf.gz"),
#'     genes_file = file.path("ref_genomes", "genes_cordinates_on_v68.txt"),
#'     annotation_file = file.path("data", "variant_calling_outputs",
#'                                 "annotation_fields.txt")
#'   )
#' }
get_variant_coverage <- function(vcf_file, genes_file, annotation_file) {
  checkmate::assert_file_exists(vcf_file, access = "r")
  checkmate::assert_file_exists(genes_file, access = "r")
  checkmate::assert_file_exists(annotation_file, access = "r")
  
  # extract the SNPs genomic coordinates (chrom and pos columns)
  if (Sys.info()["sysname"] == "Darwin") {
    genomic_coordinates <- data.table::fread(
      cmd = sprintf("bcftools query -f'%CHROM\t%POS\n' %s", vcf_file),
      nThread = 8,
      header = FALSE
    )
  } else if (Sys.info()["sysname"] == "Linux") {
    genomic_coordinates <- data.table::fread(
      cmd = sprintf("bcftools query -f'%%CHROM\t%%POS\n' %s", vcf_file),
      nThread = 8,
      header = FALSE
    )
  }
  names(genomic_coordinates) <- c("chrom", "pos")
  
  # read in the genes file
  genes <- data.table::fread(genes_file)
  
  # prepare a data frame to find overlap between the genomic coordinates and the
  # gene annotation
  pos_dt <- data.table::data.table(
    start = genomic_coordinates[["pos"]],
    end = genomic_coordinates[["pos"]]
  )
  
  # find the overlaps between the genes and the SNPs genomic coordinates
  data.table::setkey(genes, start, end)
  result <- data.table::foverlaps(
    pos_dt, genes, type = "within", nomatch = NULL
  )
  
  # remove the duplicated column introduced at line 50
  # and rename the column with the snp positions
  result[[ncol(result)]] <- NULL
  names(result)[ncol(result)] <- "pos"
  
  # join the genes to their annotations
  data.table::setkeyv(result, c("chrom", "pos"))
  annotatio_df <- data.table::fread(annotation_file)
  names(annotatio_df) <- c("annotation", "impact", "hgvs.c", "hgvs.p")
  annotatio_df <- cbind(genomic_coordinates, annotatio_df)
  result <- result |>
    dplyr::left_join(annotatio_df, by = c("chrom", "pos"))
  
  # get the count
  counts <- data.table::setDT(result)[, .N, by = .(chrom, start, end)]
  genes <- genes |>
    dplyr::left_join(counts, by = c("chrom", "start", "end")) |>
    cleanepi::standardize_column_names(
      rename = c("num_of_variants" = "N")
    )
  
  # reorder the genes by chromosome
  genes[["chrom"]] <- factor(
    genes[["chrom"]],
    levels = c("Pf3D7_01_v3", "Pf3D7_02_v3", "Pf3D7_03_v3", "Pf3D7_04_v3",
               "Pf3D7_05_v3", "Pf3D7_06_v3", "Pf3D7_07_v3", "Pf3D7_08_v3",
               "Pf3D7_09_v3", "Pf3D7_10_v3", "Pf3D7_11_v3", "Pf3D7_12_v3",
               "Pf3D7_13_v3", "Pf3D7_14_v3", "Pf3D7_MIT_v3")
  )
  genes <- genes |> dplyr::arrange(chrom)
  
  # save the per gene variant count
  data.table::fwrite(
    genes,
    file.path(dirname(vcf_file), "per_gene_variant_count.txt"),
    sep = "\t",
    nThread = 8
  )
  
  # save the all variant annotations
  data.table::fwrite(
    result,
    file.path(dirname(vcf_file), "variants_in_target_genes.txt"),
    sep = "\t",
    nThread = 8
  )
  cli::cli_alert_success("Output files {.strong per_gene_variant_count.txt} \\\
                         and {.strong variants_in_target_genes.txt} have been \\\
                         successfully saved in {.path {dirname(vcf_file)}}")
}

# tmp <- result |>
#   dplyr::select(chrom, pos, gene_description, annotation, impact, hgvs.c,
#                 hgvs.p)
# saveRDS(tmp, "data/variant_calling_outputs/snps_genomic_cordinates.RDS")
