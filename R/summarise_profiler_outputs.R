#' Summarise `malaria-profiler` outputs
#'
#' @param target_dir A character with the full path to the folder where the
#'    output files (.json files) from malaria-profiler are stored
#' @param variants_file A character with the path to the file with the confirmed
#'    mutations associated to the resistance to the drugs of interest
#' @param metadata_file A character with the full path to file with the sample
#'    metadata
#' @param mutation_type A character that determines whether to summarise the
#'    outputs for the confirmed or potential drug resistance mutations. Possible
#'    values are `confirmed` or `potential`. Default is `confirmed`
#' @param geo_localisation A logical used to determine whether to summarise the
#'    outputs about the potential region of origin of the samples. If `TRUE`, a
#'    column with inferred sample potential regions of origin will be added to
#'    the sample metadata table. Default is `TRUE`
#'
#' @returns A list with the following three elements:
#'    \enumerate{
#'      \item variants: a data frame with the genomic coordinates of the
#'          detected drug resistance mutations and `x` columns of `0` and `1`
#'          used to determine whether the sample carries a particular mutation
#'          (`1`) or no (`0`), where `x` is the number of processed samples.
#'      \item metadata: a data frame with the metadata table augmented with two
#'          columns: the sample origins, and the fraction of detected
#'          SNP markers used for geographic origin prediction
#'      \item region_probabilities: a data frame with the probabilities obtained
#'          during the geographic origin prediction. This is empty when
#'          `geo_localisation = FALSE`
#'    }
#' @export
#'
#' @examples
#' \dontrun{
#'   # establish the connection to the database
#'   results <- summarise_mp_outputs(
#'     target_dir = tempdir(),
#'     variants_file = file.path(tempdir(), "variant.csv"),
#'     metadata_file = file.path(tempdir(), "metadata.txt"),
#'     mutation_type = "confirmed",
#'     geo_localisation = "TRUE"
#'   )
#' }
summarise_mp_outputs <- function(target_dir, variants_file, metadata_file,
                                 mutation_type = c("confirmed", "potential"),
                                 geo_localisation = TRUE) {
  checkmate::assert_directory(target_dir, access = "r")
  checkmate::assert_file_exists(variants_file, access = "r")
  checkmate::assert_file_exists(metadata_file, access = "r")
  checkmate::assert_choice(mutation_type,
                           choices = c("confirmed", "potential"),
                           null.ok = FALSE)
  mutation_type <- match.arg(mutation_type)
  checkmate::assert_logical(geo_localisation, any.missing = FALSE, len = 1,
                            null.ok = FALSE)
  
  # read in the metadata file
  metadata <- rio::import(metadata_file)
  if (!"samples" %in% names(metadata)) {
    cli::cli_abort(c(
      "i" = "Column {.field samples} not found in metadata file.",
      "x" = "The column with the sample IDs must be named as {.field samples}."
    ))
  }
  metadata[["region"]] <- NA # for the inferred region
  metadata[["fraction_genotyped_geo"]] <- NA # for the inferred fraction
  
  # create a data frame to store the probabilities of inferred regions
  region_proba <- data.frame(
    samples = character(),
    ca = numeric(),
    ea = numeric(),
    sa = numeric(),
    sas = numeric(),
    seas = numeric(),
    wa = numeric()
  )
  names(region_proba) <- c("samples", "Central Africa", "Eastern Africa",
                           "South America", "South Asia", "Southeast Asia",
                           "Western Africa")
  
  # get the list of output files
  target_files <- list.files(
    path = target_dir,
    full.names = TRUE,
    pattern = "*.json"
  )
  
  # get the sample names
  samples <- gsub(".results.json", "", basename(target_files))
  
  # read in the variants file
  variants <- read_confirmed_variants(variants_file = variants_file)
  
  # loop through the files and summarise the malaria-profiler output
  other_dr <- NULL
  for (file in target_files) {
    # read in the JSON file
    json <- jsonlite::read_json(file)
    
    # add a column for the target sample
    # all values are set to '0' by default
    sample_name <- gsub(".results.json", "", basename(file))

    # add resistance mutations found in the target sample to the
    # resistance_mutations data frame
    resistance_mutations <- switch(mutation_type,
      confirmed = get_confirmed_dr_variants(
        json = json,
        variants = variants,
        sample_name = sample_name
      ),
      potential = get_other_dr_variants(
        nested_list = json[["other_variants"]],
        sample_name = sample_name,
        tmp_others_dr = other_dr
      )
    )
    other_dr <- variants <- resistance_mutations
    
    # extract the geo-localisation information for the target sample
    if (geo_localisation) {
      geo_localisation_data <- get_inferred_region(
        json = json,
        metadata = metadata,
        region_proba = region_proba,
        sample_name = sample_name
      )
      metadata <- geo_localisation_data[["metadata"]]
      region_proba <- geo_localisation_data[["region_proba"]]
    }
    
  }
  
  # return a list with the 3 elements
  return(list(
    variants = variants,
    metadata = metadata,
    region_probabilities = region_proba
  ))
}




