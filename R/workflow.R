#' Perform computations on imported and processed data
#'
#' @param input output of \code{\link{process_input}}
#' @param anno output of \code{\link{import_annotations}}
#' @param data output of \code{\link{import_database}} or
#'   \code{\link{import_msigdb}}
#' @param universe number of genes in universe
#' @param info_from optional: \code{"annotation"} or (default) \code{"database"}
#' @param anno_opts optional: \code{"name"} for gene set names, \code{"syms"}
#'   for gene set symbols, \code{"info"} for gene set descriptions,
#'   \code{"auto"} for automatically generated annotations and/or (default)
#'   \code{"file"} for manual annotations
#' @param categories optional: categories to include; default all
#' @param organisms optional: organisms to include; default all
#'
#' @return
#' \code{stats} tibble: overlap statistics
#'
#' \code{matches} list: names: annotations vector: matched genes
#' @export
#'
#' @importFrom magrittr %>%
#' @examples
#' anno_path <- system.file("extdata", "ex_anno.csv", package = "glacier")
#' anno <- import_annotations(anno_path, ",", FALSE, c(2, 10), 11)
#' data_path <- system.file("extdata", "ex_data.csv", package = "glacier")
#' data <- import_database(data_path, ",", FALSE, c(2, 4), 0)
#'
#' input <- process_input('CYP1A1 0.2 CYP1B1 NQO1 0.3 SODD 9.0')
#' results <- compute(input, anno, data, 10000)
compute <- function(input, anno, data, universe, info_from = "database",
                    anno_opts = "file", categories = NULL, organisms = NULL) {
  source <- if (info_from == "annotations") anno else data$gs_info
  info <- source %>% dplyr::select("name", "info")

  anno_proc <- process_annotations(anno, info, anno_opts)
  data_proc <- process_database(data, categories, organisms)

  # calculate statistics
  calculate(input, anno_proc$annos, anno_proc$gs_annos, data_proc$gs_genes,
            10000)
}

#' Start web application
#'
#' @export
webstart <- function() {
  if (!requireNamespace("shiny", quietly = F) ||
      !requireNamespace("shinyjs", quietly = F) ||
      !requireNamespace("shinythemes", quietly = F) ||
      !requireNamespace("shinyWidgets", quietly = F)) {
    stop("Shiny-related packages are not available, install and try again")
  }

  shiny::runApp(system.file("shiny", package = "glacier"))
}
