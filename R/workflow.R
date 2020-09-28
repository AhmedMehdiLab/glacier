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
#' input <- process_input_text('CYP1A1 0.2 CYP1B1 NQO1 0.3 SODD 9.0')
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
  if (!requireNamespace("shiny", quietly = F) ||
      !requireNamespace("shinyjs", quietly = F) ||
      !requireNamespace("shinythemes", quietly = F) ||
      !requireNamespace("shinyWidgets", quietly = F))
    stop("Packages 'shiny', 'shinyjs', 'shinythemes' and 'shinyWidgets' are \
         required for this feature")

  if (!requireNamespace("Seurat", quietly = F))
    warning("Package 'Seurat' is required for single-cell analysis")

  if (!requireNamespace("limma", quietly = F))
    message("Package 'limma' from Bioconductor is recommended to accelerate \
            single-cell analysis")

  shiny::runApp(system.file("shiny", package = "glacier"), launch.browser = T)
}
