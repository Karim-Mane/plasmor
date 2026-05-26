# required packages
# jsonlite, cli, dplyr, purrr, rio, checkmate, cleanepi

#' Extract drug name from variants file
#'
#' @param x A character that corresponds to an entry of the Info column in the
#'    variant file
#'
#' @returns A character with the drug name
#' @keywords internal
get_drug <- function(x) {
  drug <- unlist(strsplit(x, ";", fixed = TRUE))[[2]]
  return(unlist(strsplit(drug, "=", fixed = TRUE))[[2]])
}

#' Add a drug resistance profile for the specified sample to the data frame of
#' drug resistance profiles based on potential drug resistance mutations
#'
#' @param dat A data frame with the resistance profile of the target sample
#'    based on the potential drug resistance mutations
#' @inheritParams get_other_dr_variants
#'
#' @returns A data frame potential drug resistance profile augmented with one
#'    column corresponding the profile of the target sample
#' @export
#'
#' @examples
construct_dr_matrix <- function(dat, sample_name, tmp_others_dr) {
  # select the columns of interest
  dat <- dat |>
    dplyr::select(chrom, pos, ref, alt, gene_name, gene_id, protein_change) |>
    dplyr::rename(mutation = "protein_change") |>
    dplyr::distinct(.keep_all = TRUE)
  
  # if first iteration, simply add the sample column
  if (is.null(tmp_others_dr)) {
    tmp_others_dr <- dat
    tmp_others_dr[[sample_name]] <- 1
  } else {
    # if not first iteration:
    # check for common mutations between the new sample and the ones that are
    # already stored
    # if any new mutation, add it
    current_samples <- names(tmp_others_dr)[-c(1:7)]
    tmp_others_dr[[sample_name]] <- 0
    target <- as.character(
      apply(
        cbind(tmp_others_dr[["chrom"]], tmp_others_dr[["pos"]], tmp_others_dr[["mutation"]]),
        1,
        paste,
        collapse = "_"
      )
    )
    test <- as.character(
      apply(
        cbind(dat[["chrom"]], dat[["pos"]], dat[["mutation"]]),
        1,
        paste,
        collapse = "_"
      )
    )
    idx <- match(dat$mutation, tmp_others_dr$mutation)
    tmp_others_dr[[sample_name]][idx[!is.na(idx)]] <- 1
    if (any(is.na(idx))) {
      dat <- dat |>
        dplyr::mutate(sample = 1) |>
        dplyr::mutate(idx = idx) |>
        dplyr::filter(is.na(idx)) |>
        dplyr::mutate(idx = NULL)
      names(dat)[ncol(dat)] <- sample_name
      for (sample in current_samples) {
        dat[[sample]] <- 0
      }
      dat <- dat[, names(tmp_others_dr)]
      tmp_others_dr <- rbind(tmp_others_dr, dat)
    }
  }
  
  return(tmp_others_dr)
}

#' Get the sample profiles for the potential drug resistance mutations
#'
#' @inheritParams get_confirmed_dr_variants
#' @param nested_list A list with the information on the potential drug
#'    resistance mutations found on the given sample
#' @param tmp_others_dr A data frame used to store the sample profiles
#'
#' @returns A data frame with the sample profiles augmented with a new column
#'    corresponding to the profile of the specified sample
#' @export
#'
#' @examples
get_other_dr_variants <- function(nested_list, sample_name, tmp_others_dr) {
  # safe extraction — when no other drug resistance variants are found for a
  # given sample, fill it with NA
  df_scalars <- nested_list |>
    purrr::map_dfr(function(x) {
      scalars <- x[!sapply(x, is.list)]
      
      # replace NULL or empty values with NA
      scalars <- lapply(scalars, function(v) {
        if (is.null(v) || length(v) == 0) NA else v
      })
      
      as.data.frame(scalars)
    })
  
  # for the nested list of 8 elements
  df_nested <- nested_list |>
    purrr::map_dfr(function(x) {
      nested <- x[sapply(x, is.list)]
      as.data.frame(t(unlist(nested)))
    })
  
  # combine the two side by side
  df_final <- dplyr::bind_cols(df_scalars, df_nested)
  num_mutations <- nrow(df_final)
  num_genes <- length(unique(df_final[["gene_name"]]))
  cli::cli_alert_info(
    "Found {cli::no(num_mutations)} mutation{?s} across {cli::no(num_genes)} \\\
    protential drug resistance gene{?s} in sample {.val {sample_name}}."
  )
  
  # if no other mutations is found for the target sample, add a new column for
  # the sample and set all values at 0
  if (num_mutations == 0) {
    tmp_others_dr[[sample_name]] <- 0
    return(tmp_others_dr)
  }
  
  # when mutations are found, account for them by considering the ones that are
  # already detected from the processed samples
  tmp_others_dr <- construct_dr_matrix(df_final, sample_name, tmp_others_dr)
  
  return(tmp_others_dr)
}


#' Extract the prediction region of origin for the specified sample
#'
#' @inheritParams get_confirmed_dr_variants
#' @inheritParams summarise_mp_outputs 
#' @param region_proba A data frame used to store the prediction probabilities
#'    for each region
#'
#' @returns A list with the updated sample metadata and region probabilities
#'    data frame
#' @export
#'
#' @examples
get_inferred_region <- function(json, metadata, region_proba, sample_name) {
  # populate the inferred region column in the sample metadata
  regions <- json[["geo_classification"]][["probabilities"]] |>
    dplyr::bind_rows()
  # when the region was not inferred, set it to NA in the metadata
  idx <- match(sample_name, metadata[["samples"]])
  metadata[["fraction_genotyped_geo"]][idx] <- json[["geo_classification"]][["fraction_genotyped"]]
  if (nrow(regions) == 0) {
    metadata[["region"]][idx] <- NA
  } else {
    regions <- regions |>
      dplyr::mutate(probability = format(probability, scientific = FALSE))
    metadata[["region"]][idx] <- regions[["region"]][which.max(regions[["probability"]])]
    
    # populate the data with the region probabilities
    regions <- as.data.frame(t(regions))
    regions <- regions[-1, ]
    rownames(regions) <- NULL
    regions <- regions |>
      dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))
    names(regions) <- names(region_proba)[-1]
    regions <- regions |>
      dplyr::mutate(samples = sample_name) |>
      dplyr::relocate(samples, .before = `Central Africa`)
    region_proba <- rbind(region_proba, regions)
  }
  
  return(list(
    metadata = metadata,
    region_proba = region_proba
  ))
}

#' Summarise the malaria-profiler outputs for confirmed drug resistance
#' mutations
#'
#' @param json A list object obtained from reading an output file in JSon format
#' @param variants A data frame that will be used to store the variant profile
#'    for each sample
#' @param sample_name A character with a sample name
#'
#' @returns An updated variant table where the profile of the specified sample
#'    is added
#' @export
#'
#' @examples
get_confirmed_dr_variants <- function(json, variants, sample_name) {
  # add a column for the target sample
  # all values are set to '0' by default
  variants[[sample_name]] <- 0
  
  # identify the drug resistance mutations found on the target sample
  # and set the corresponding row to 1 (to mark its presence on that sample)
  resistance_mutations <- json[["dr_variants"]] |>
    dplyr::bind_rows()
  num_mutations <- nrow(resistance_mutations)
  cli::cli_alert_info(
    "Found {cli::no(num_mutations)} confirmed drug resistance mutation{?s} \\\
    in sample {.val {sample_name}}."
  )
  if (num_mutations == 0) {
    return(variants) # next
  }
  idx <- match(resistance_mutations[["change"]], variants[["mutation"]])
  variants[[sample_name]][idx] <- 1
  selection_table <- data.frame(
    idx = idx,
    chrom = resistance_mutations[["chrom"]],
    pos = resistance_mutations[["pos"]],
    ref = resistance_mutations[["ref"]],
    alt = resistance_mutations[["alt"]],
    gene_name = resistance_mutations[["gene_name"]]
  )
  are_missing <- which(is.na(variants[["chrom"]][idx]))
  if (length(are_missing) > 0) {
    selection_table <- selection_table[are_missing, ]
    variants[["chrom"]][selection_table[["idx"]]] <- selection_table[["chrom"]]
    variants[["pos"]][selection_table[["idx"]]] <- selection_table[["pos"]]
    variants[["ref"]][selection_table[["idx"]]] <- selection_table[["ref"]]
    variants[["alt"]][selection_table[["idx"]]] <- selection_table[["alt"]]
    variants[["gene_name"]][selection_table[["idx"]]] <- selection_table[["gene_name"]]
  }
  
  return(variants)
}