#' Choose the target images based on the values of the function's arguments
#'
#' @inheritParams convert_to_png
#'
#' @returns A character vector of the path to the images to be converted
#' @keywords internal
choose_files <- function(images_path, images) {
  # abort if both arguments are NULL
  if (all(is.null(c(images_path, images)))) {
    cli::cli_abort(c(
      i = "No argument was provided.",
      x = "You must provide either a {.cls vector} of the target images with \\\
      the {.emph images} argument or the path to the folder where the images \\\
      are stored using the {.emph images_path} argument."
    ))
  }
  # when both image path and images are provided, only consider the target
  # images
  if (!is.null(images) && !is.null(images_path)) {
    cli::cli_inform(c(
      "!" = "You provided both {.emph images_path} and {.emph images}.",
      i = "Only images in {.emph images} will be considered."
    ))
  }
  
  # get image list from directory
  target_files <- images
  if (is.null(target_files)) {
    # list all .heic, .jpg, .jpeg, .tif, .tiff, .webp files
    target_files <- list.files(
      path = images_path,
      pattern = "\\.(heic|jpg|jpeg|tif|tiff|webp)$",
      ignore.case = TRUE,
      full.names = TRUE
    )

    # send an error message if no file is found
    if (length(target_files) == 0) {
      cli::cli_abort(c(
        i = "Could not find any file to convert in {.url {images_path}}",
        x = "Target images must have one of the following extensions: \\\
        {.strong heic}, {.strong jpg}, {.strong jpeg}, {.strong tif}, \\\
        {.strong tiff}, {.strong webp},"
      ))
    }
  }
  return(target_files)
}

#' Convert images to PNG format
#'
#' @param images_path A character with the full path to the directory where the
#'    images are stored
#' @param images A character vector with the path to the images to be converted
#'
#' @returns Invisibly return a vector of the converted images
#' @export
#'
#' @examples
#' \dontrun{
#'   # convert a specific image
#'   convert_to_png(
#'     images <- system.file("images", "Image003.heic", package = "plasmor")
#'   )
#' }
convert_to_png <- function(images = NULL, images_path = NULL) {
  if (!is.null(images_path)) {
    checkmate::assert_directory(images_path, access = "r")
  }
  checkmate::assert_vector(images, unique = TRUE, null.ok = TRUE,
                           any.missing = FALSE)
  if (!is.null(images)) {
    if (any(!file.exists(images))) {
      cli::cli_abort(c(
        i = "Unable to find some files specified in argument {.emph images}.",
        "!" = "Did you specify the correct image file paths ?",
        x = "{.emph images} must be a {.cls vector} of {.cls character} with \\\
        the full path the target images."
      ))
    }
  }

  # get the vector of target images
  target_files <- choose_files(images_path, images)

  # loop through each detected file and convert it into PNG
  png_images <- NULL
  for (file in target_files) {
    cli::cli_progress_step("converting file {.file {file}}")
    img <- magick::image_read(file)
    output_file <- sub("\\.(heic|jpg|jpeg|tif|tiff|webp)$", ".png", file,
                       ignore.case = TRUE)
    magick::image_write(img, output_file, format = "png")
    unlink(file)
    png_images <- c(png_images, output_file)
  }
  
  return(invisible(png_images))
}
  