#' Parse Fastq-Screen output file and check coverage of the target genome
#'
#' @param output_file A character with the path to the Fastq-Screen output file
#'
#' @returns A list of two elements: A numeric with number of reads that do not
#'    map to any of the specified genome, and a data frame with screening result
#' @keywords internal
parse_output_file <- function(output_file) {
  # read all lines
  lines <- readLines(output_file)
  
  # get the proportion of reads that do not map to any of the genomes
  no_hits <- as.numeric(unlist(
    strsplit(lines[startsWith(lines, "%")], ":", fixed = TRUE)
  )[[2]])
  
  # get the mapping stats 
  lines <- lines[!startsWith(lines, "#") | !startsWith(lines, "%")]
  fields <- strsplit(lines, "\t")
  fields <- fields[lengths(fields) >= 12]
  res <- data.frame(
    sapply(fields, `[[`, 1),
    sapply(fields, `[[`, 2),
    sapply(fields, `[[`, 3),
    sapply(fields, `[[`, 4),
    sapply(fields, `[[`, 5),
    sapply(fields, `[[`, 6),
    sapply(fields, `[[`, 7),
    sapply(fields, `[[`, 8),
    sapply(fields, `[[`, 9),
    sapply(fields, `[[`, 10),
    sapply(fields, `[[`, 11),
    sapply(fields, `[[`, 12),
    stringsAsFactors = FALSE
  )
  headers <- res[1, ]
  res <- res[-1, ]
  names(res) <- headers

  return(list(
    no_hits = no_hits,
    hits = res
  ))
}


#' Parse Fastq-Screen output file to extract genome mapping statistics.
#' 
#' If the coverage of the target genome is less than 50% or the percent of
#' reads mapping to none of the genomes is greater than 5%, a comment is
#' associated to that read.
#' Also, the %reads mapping to the Human genome is reported. Files where there
#' are reads mapping uniquely to other Plasmodium genome will have the name of
#' these organisms in the 'others' column.
#'
#' @param dir A character with the path to the folder where the Fastq-Screen
#'    files are stored
#' @param target_genome A character with the name of the target genome as
#'    specified in the Fastq-Screen config file
#'
#' @returns A data frame with seven columns when there are samples with a
#'    coverage is less than 50% or the percent of no hits is greater than 5%
#' @export
#'
#' @examples
#' \dontrun{
#'   screening_summary <- summarise_screening_outputs(
#'     dir = "/Volumes/NO NAME/fastq_screen_outputs",
#'     target_genome = "falciparum"
#'   )
#' }
summarise_screening_outputs <- function(dir, target_genome) {
  checkmate::assert_directory_exists(dir, access = "r")
  checkmate::assert_character(target_genome, any.missing = FALSE, len = 1,
                              null.ok = FALSE)
  files <- list.files(dir, pattern = "*_screen.txt", full.names = TRUE)
  if (!length(files)) {
    cli::cli_abort(c(
      x = "No file with a {.val _screen.txt} suffix find",
      i = "Please make sure to specify the correct directory"
    ))
  }

  # loop over files to collect the screening stats
  alert <- 0
  output <- NULL
  cli::cli_progress_bar("Checking files", total = length(files))
  for (i in seq_len(length(files))) {
    Sys.sleep(4 / length(files))
    output_file = files[i]
    tmp_res <- parse_output_file(output_file = output_file)
    no_hits <- tmp_res[["no_hits"]]
    hits <- tmp_res[["hits"]]
    hits[["%One_hit_one_genome"]] <- as.numeric(hits[["%One_hit_one_genome"]])
    hits[["%Multiple_hits_one_genome"]] <- as.numeric(hits[["%Multiple_hits_one_genome"]])
    hits[["#Reads_processed"]] <- as.numeric(hits[["#Reads_processed"]])
    # check if target genome name is correct
    if (!target_genome %in% hits[["Genome"]]) {
      cli::cli_abort(c(
        x = "Specified target genome {.val {target_genome}} was not used \\\
      during the mapping process.",
        i = "Please make sure to use the same name in the Fastq-screen config \\\
      file."
      ))
    }

    # extract the stats for the target genome
    # report the sample if the %mapped reads to the target genome is <= 50% or
    # %not-mapped reads is >= 5%
    sample_name <- gsub("_screen.txt", "", basename(output_file))
    comment <- "PASS"
    idx <- which(hits[["Genome"]] == target_genome)
    target <- as.numeric(hits[idx, ]["%One_hit_one_genome"])
    target_multi <- as.numeric(hits[idx, ]["%Multiple_hits_one_genome"])
    all_target_hits <- target + target_multi
    num_processed_reads <- as.numeric(hits[idx, ]["#Reads_processed"])
    if ((target + target_multi) <= 50 || no_hits >= 5) {
      comment <- dplyr::case_when(
        all_target_hits <= 50 & no_hits >= 5 ~ "Less than 50% of the reads map to the target genome and there is more than 5% of reads that do not map to any
      genome.",
        all_target_hits <= 50 ~ "Less than 50% of the reads map to the target genome.",
        no_hits >= 5 ~ "There is more than 5% of reads that do not map to any genome."
      )
      alert <- alert + 1
    }

    # extract the stats for the Human genome
    idx_human <- which(hits[["Genome"]] == "Human")
    human_hit <- as.numeric(
      hits[idx_human, ]["%One_hit_one_genome"] +
        hits[idx_human, ]["%Multiple_hits_one_genome"]
    )

    # extract the stats for the remaining genome
    # these will be reported if the %hits to a specific genome > 0
    others <- NA
    tmp <- hits[hits[["Genome"]] != "Human" & hits[["Genome"]] != target_genome, ]
    if (any(tmp[["%One_hit_one_genome"]] > 0)) {
      others <- toString(tmp[["Genome"]][tmp[["%One_hit_one_genome"]] > 0])
    }
    
    # build the output data frame
    output <- rbind(
      output,
      data.frame(
        sample = sample_name,
        processed_reads = num_processed_reads,
        no_hits = no_hits,
        hits = all_target_hits,
        comment = comment,
        human = human_hit,
        others = others,
        stringsAsFactors = FALSE
      )
    )
    cli::cli_progress_update()
  }
  cli::cli_progress_done()
  
  if (length(alert)) {
    cli::cli_inform(c(
      i = "Found {.val {alert}} files to be looked at carefully."
    ))
  }
  
  return(output)
}


#' Plot the Fastq-Screen output for a specific sample
#'
#' @inheritParams summarise_screening_outputs dir 
#' @param sample A character with the name of a target sample
#'
#' @returns Invisibly returns the screening output statistics for the specified
#'    sample
#' @export
#'
#' @examples
#' \dontrun{
#'   plot_sample_screening(
#'     dir = "/Volumes/NO NAME/fastq_screen_outputs",
#'     sample = "29504"
#'   )
#' }
plot_sample_screening <- function(dir, sample) {
  checkmate::assert_directory_exists(dir, access = "r")
  checkmate::assert_character(sample, any.missing = FALSE, len = 1,
                              null.ok = FALSE)
  files <- list.files(dir, pattern = "*_screen.txt", full.names = TRUE)
  idx <- which(grepl(sample, basename(files), fixed = TRUE))
  if (!length(idx)) {
    cli::cli_abort(c(
      x = "Found no file associated to the provided sample ID {.val {sample}}",
      i = "Please make sure to specify the correct directory and sample ID"
    ))
  }
  
  files <- files[idx]
  output <- NULL
  for (file in files) {
    # parse a target file
    tmp_res <- parse_output_file(output_file = file)
    no_hits <- tmp_res[["no_hits"]]
    hits <- data.frame(t(tmp_res[["hits"]]))
    names(hits) <- as.character(hits[1, ])
    hits <- hits[-1, ]
    
    # add a column and transform the data
    hits <- cbind(
      stats = rownames(hits),
      hits
    )
    rownames(hits) <- NULL
    idx <- which(
      hits[["stats"]] %in% c("%One_hit_one_genome", "%Multiple_hits_one_genome",
                             "%One_hit_multiple_genomes",
                             "%Multiple_hits_multiple_genomes", "No hits")
    )
    hits <- hits[idx, ]
    hits <- cbind(
      hits,
      "No hits" = as.character(c(no_hits, rep(NA, (nrow(hits) - 1))))
    )
    hits <- tidyr::pivot_longer(
      hits,
      cols = names(hits)[!names(hits) == "stats"]
    )
    hits[["value"]] <- as.numeric(hits[["value"]])
    hits[["file"]] <- gsub("_screen.txt", "", basename(file))
    output <- rbind(output, hits)
  }
  
  # define your manual colors
  idx <- which(output[["name"]] == "No hits")
  output[["stats"]][idx] <- "No_hit"
  manual_colors <- c(
    "%One_hit_one_genome" = "lightblue",
    "%Multiple_hits_one_genome" = "steelblue",
    "%One_hit_multiple_genomes" = "salmon",
    "%Multiple_hits_multiple_genomes" = "red",
    "No_hit" = "darkgray"
  )
  
  output[["stats"]] <- factor(
    output[["stats"]],
    levels = c("No_hit", "%Multiple_hits_multiple_genomes",
               "%One_hit_multiple_genomes", "%Multiple_hits_one_genome",
               "%One_hit_one_genome")
  )
  output[["name"]] <- factor(
    output[["name"]],
    levels = c("Human", "falciparum", "vivax", "malariae",
               "ovale_curtisi", "ovale_wallikeri", "No hits")
  )
  
  # build ggplot with tooltip aesthetic
  p <- ggplot(output, aes(
    fill = stats,
    y = value,
    x = name,
    text = paste0(
      "<b>", stats, "</b><br>",
      "Proportion: ", round(value, 2), "%"
    )
  )) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = manual_colors) +
    facet_wrap(~ file) +
    theme_minimal() +
    labs(x = "target genomes", y = "% mapped reads") +
    theme(
      axis.title  = element_text(face = "bold"),
      axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5),
      axis.text.y = element_text(size = 8)
    )
  
  # convert to plotly
  p <- plotly::ggplotly(p, tooltip = "text")
  # print(p)
  
  # invisibly return the output data frame
  return(invisible(p))
}
