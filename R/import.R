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
#' \dontrun{
#' read.common('path/to/common.csv', ',', T)
#' }
read.common <- function(path, delim, header) {
  readr::read_delim(path, delim, col_types = readr::cols(.default = readr::col_character()), 
    col_names = header, trim_ws = T)
}

#' Pre-process common file type tibble
#'
#' @param common.raw raw file tibble
#' @param cols data column boundaries
#' @param info info column
#'
#' @return pre-processed common file type tibble
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' common.raw <- read.common('path/to/common.csv', ',', T)
#' import.common(common.raw, c(2, 9), 10)
#' }
import.common <- function(common.raw, cols, info) {
  if (is.null(cols)) {
    start <- stop <- NULL
  } else {
    start <- cols[1]
    stop <- cols[2]
  }
  
  # validate
  max <- ncol(common.raw)
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
    info <- common.raw %>% dplyr::pull(info) %>% replace(is.na(.), "")
  }
  
  # process
  common.raw %>% dplyr::select(name = 1, data. = dplyr::all_of(start:stop)) %>% 
    tibble::add_column(info = info)
}

#' Import annotation file
#'
#' @param anno.raw raw file tibble or output of \code{\link{import.common}}
#' @param cols data column boundaries
#' @param info info column
#'
#' @return gene set annotations and information
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' anno.raw <- read.common('path/to/anno.csv', ',', T)
#' import.annotations(anno.raw, c(2, 9), 10)
#' }
import.annotations <- function(anno.raw, cols = NULL, info = NULL) {
  if (dim(anno.raw)[1] == 0 || dim(anno.raw)[2] == 0) {
    stop("anno.raw is empty")
  }
  
  anno.raw %>% import.common(cols, info) %>% dplyr::rename(anno. = dplyr::starts_with("data."))
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
#' \dontrun{
#' import.msigdb_xml('path/to/msigdb.xml')
#' }
import.msigdb_xml <- function(path) {
  data <- path %>% xml2::read_xml() %>% xml2::xml_children()
  catc <- data %>% xml2::xml_attr("CATEGORY_CODE")
  cats <- data %>% xml2::xml_attr("SUB_CATEGORY_CODE")
  
  # extract information
  info <- tibble::tibble(name = data %>% xml2::xml_attr("STANDARD_NAME"), info = data %>% 
    xml2::xml_attr("DESCRIPTION_BRIEF"), desc = data %>% xml2::xml_attr("DESCRIPTION_FULL") %>% 
    as.factor(), category = stringr::str_c(catc, cats, sep = " ") %>% stringr::str_squish() %>%
    as.factor(), organism = data %>% xml2::xml_attr("ORGANISM") %>% as.factor())
  
  # extract contents
  genes <- data %>% xml2::xml_attr("MEMBERS_SYMBOLIZED") %>% stringr::str_split(",") %>% 
    purrr::map(~.[. != ""]) %>% purrr::set_names(info$name)
  
  list(gs_genes = genes, gs_info = info)
}

#' Import database file
#'
#' @param data.raw raw file tibble or output of \code{\link{import.common}}
#' @param data data column boundaries
#' @param info info column
#'
#' @return gene set contents and information
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' data.raw <- read.common('path/to/data.csv', ',', T)
#' import.database(data.raw, c(2, 9), 10)
#' }
import.database <- function(data.raw, data = NULL, info = NULL) {
  if (dim(data.raw)[1] == 0 || dim(data.raw)[2] == 0) {
    stop("data.raw is empty")
  }
  
  none <- as.factor("Not assigned")
  proc <- data.raw %>% import.common(data, info) %>% tibble::add_column(category = none, 
    organism = none)
  
  # extract
  info <- proc %>% dplyr::select(-dplyr::starts_with("data."))
  genes <- proc %>% dplyr::select(dplyr::starts_with("data.")) %>% purrr::pmap(c, 
    use.names = F) %>% purrr::map(~.[!is.na(.)]) %>% purrr::set_names(info$name)
  
  list(gs_genes = genes, gs_info = info)
}
