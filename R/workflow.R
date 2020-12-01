#' Perform computations on imported and processed data
#'
#' @param input output of \code{\link{process_input_text}} or
#'   \code{\link{process_input_seurat}}
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
#' @param save optional: path to save overlap statistics as \code{.csv}
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
#' data_path <- system.file("extdata", "ex_data.csv", package = "glacier")
#' anno <- import_annotations(anno_path, ",", TRUE, c(2, 4), 5)
#' data <- import_database(data_path, ",", FALSE, c(2, 4), 0)
#'
#' input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
#' results <- compute(input, anno, data)
compute <- function(input, anno, data, universe = NULL, info_from = "database",
                    anno_opts = "file", categories = FALSE, organisms = FALSE,
                    save = NULL) {
  source <- if (info_from == "annotations") anno else data$gs_info
  info <- source %>% dplyr::select("name", "info")

  anno_proc <- process_annotations(anno, info, anno_opts)
  data_proc <- process_database(data, categories, organisms)

  # calculate statistics
  calc <- calculate(input, anno_proc$annos, anno_proc$gs_annos,
                    data_proc$gs_genes, universe)

  if (!is.null(save)) readr::write_csv(calc$stats, save)
  return(calc)
}

#' Launch web interface
#'
#' @export
#'
#' @examples
#' vignette("web-app", package = "glacier")
webstart <- function() {
  uses(c("shiny", "shinyjs", "shinythemes", "shinyWidgets"), stop,
       "'shiny', 'shinyjs', 'shinythemes' and 'shinyWidgets' are required")
  uses("Seurat", warning, "'Seurat' is required to handle Seurat data")
  uses("limma", message, "'limma' is recommended to accelerate analysis")

  shiny::runApp(system.file("shiny", package = "glacier"),
                host = "0.0.0.0",
                launch.browser = T)
}
