#' Read delimited files into tibble
#'
#' @param path path to file
#' @param delim file delimiter
#' @param header whether file contains header
#'
#' @return raw file contents as tibble
#' @export
#'
#' @examples
#' file <- system.file("extdata", "ex_anno.csv", package = "glacier")
#' common_raw <- read_common(file, ",", FALSE)
read_common <- function(path, delim, header) {
  readr::read_delim(
    path, delim, col_types = readr::cols(.default = readr::col_character()),
    col_names = header, trim_ws = T
  )
}

#' Pre-process common file type tibble
#'
#' @param common_raw raw file tibble
#' @param cols data column boundaries
#' @param info info column
#'
#' @return pre-processed common file type tibble
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @examples \dontrun{
#' file <- system.file("extdata", "ex_anno.csv", package = "glacier")
#' common_raw <- read_common(file, ",", FALSE)
#' common_pre <- import_common(common_raw, c(2, 9), 10)
#' }
import_common <- function(common_raw, cols, info) {
  if (!dim(common_raw)[1] * dim(common_raw)[2]) stop("File is empty")

  if (is.null(cols)) {
    start <- stop <- NULL
  } else {
    start <- cols[1]
    stop <- cols[2]
  }

  # validate
  max <- ncol(common_raw)
  if (is.null(start) || !dplyr::between(start, 1, max)) {
    start <- min(2, max)
  }
  if (is.null(stop) || !dplyr::between(stop, 1, max)) {
    stop <- max
  }
  if (is.null(info) || !dplyr::between(info, 1, max)) {
    info <- character(1)
  } else {
    . <- NULL
    info <- common_raw %>% dplyr::pull(info) %>% replace(is.na(.), "")
  }

  # process
  common_raw %>%
    dplyr::select(
      name = if (ncol(common_raw)) 1,
      data_ = dplyr::all_of(start:stop)
    ) %>%
    tibble::add_column(info = info)
}

#' Import annotation file
#'
#' @param anno_raw raw file tibble or output of \code{\link{import_common}}
#' @param cols data column boundaries
#' @param info info column
#'
#' @return gene set annotations and information
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' file <- system.file("extdata", "ex_anno.csv", package = "glacier")
#' anno_raw <- read_common(file, ",", FALSE)
#' anno_pre <- import_annotations(anno_raw, c(2, 9), 10)
import_annotations <- function(anno_raw, cols = NULL, info = NULL) {
  anno_raw %>%
    import_common(cols, info) %>%
    dplyr::rename(anno_ = dplyr::starts_with("data_"))
}

#' Import MSigDB XML file
#'
#' @param path path to file
#'
#' @return gene set contents and information
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' file <- system.file("extdata", "ex_msig.xml", package = "glacier")
#' data_pre <- import_msigdb_xml(file)
import_msigdb_xml <- function(path) {
  data <- path %>% xml2::read_xml() %>% xml2::xml_children()
  catc <- data %>% xml2::xml_attr("CATEGORY_CODE")
  cats <- data %>% xml2::xml_attr("SUB_CATEGORY_CODE")

  # extract information
  info <- tibble::tibble(
    name = data %>% xml2::xml_attr("STANDARD_NAME"),
    info = data %>% xml2::xml_attr("DESCRIPTION_BRIEF"),
    desc = data %>% xml2::xml_attr("DESCRIPTION_FULL") %>% as.factor(),
    category = stringr::str_c(catc, cats, sep = " ") %>%
      stringr::str_squish() %>%
      as.factor(),
    organism = data %>% xml2::xml_attr("ORGANISM") %>% as.factor()
  )

  # extract contents
  genes <- data %>%
    xml2::xml_attr("MEMBERS_SYMBOLIZED") %>%
    stringr::str_split(",") %>%
    purrr::map(~.[. != ""]) %>%
    purrr::set_names(info$name)

  list(gs_genes = genes, gs_info = info)
}

#' Import database file
#'
#' @param data_raw raw file tibble or output of \code{\link{import_common}}
#' @param data data column boundaries
#' @param info info column
#'
#' @return gene set contents and information
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' file <- system.file("extdata", "ex_data.csv", package = "glacier")
#' data_raw <- read_common(file, ",", FALSE)
#' data_pre <- import_database(data_raw, c(2, 4), 0)
import_database <- function(data_raw, data = NULL, info = NULL) {
  none <- as.factor("Not assigned")
  proc <- data_raw %>%
    import_common(data, info) %>%
    tibble::add_column(category = none, organism = none)

  # extract
  info <- proc %>% dplyr::select(-dplyr::starts_with("data_"))
  genes <- proc %>%
    dplyr::select(dplyr::starts_with("data_")) %>%
    purrr::pmap(c, use.names = F) %>%
    purrr::map(~.[!is.na(.)]) %>%
    purrr::set_names(info$name)

  list(gs_genes = genes, gs_info = info)
}
