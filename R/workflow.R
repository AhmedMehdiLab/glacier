#' Calculate overlap statistics for select genes
#'
#' @param input_raw string input with genes and values
#' @param anno_pre output of \code{\link{import_annotations}}
#' @param anno_opts 'name' for gene set names, 'syms' for gene set symbols,
#'   'info' for gene set information, 'auto' for automatically generated
#'   annotations, 'file' for manual annotations
#' @param data_pre output of \code{\link{import_database}}
#' @param info_from "annotation" or "database"
#' @param categories optionally filter categories
#' @param organisms optionally filter organisms
#' @param universe number of genes in universe
#'
#' @return overlap statistics
#' @export
#'
#' @examples
#' anno_file <- system.file("extdata", "ex_anno.csv", package = "glacier")
#' anno_raw <- read_common(anno_file, ",", FALSE)
#' anno_pre <- import_annotations(anno_raw, c(2, 9), 10)
#'
#' data_file <- system.file("extdata", "ex_data.csv", package = "glacier")
#' data_raw <- read_common(data_file, ",", FALSE)
#' data_pre <- import_database(data_raw, c(2, 4), 0)
#'
#' calc <- workflow('CYP1A1 0.2 CYP1B1 NQO1 0.3 SODD 9.0', anno_pre,
#'                  data_pre = data_pre, universe = 100)
workflow <- function(input_raw, anno_pre, anno_opts = "file", data_pre,
                     info_from = "database", categories = NULL,
                     organisms = NULL, universe) {
  input <- process_input(input_raw)
  info <- if (info_from == "database")
    data_pre$gs_info %>% dplyr::select("name", "info")
  else anno_pre %>% dplyr::select("name", "info")
  anno <- process_annotations(anno_pre, info, anno_opts)
  data <- process_database(data_pre, categories, organisms)

  calculate(input$input, anno$annos, anno$gs_annos, data$gs_genes, universe)
}
