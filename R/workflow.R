#' Calculate overlap statistics for select genes
#'
#' @param input.raw string input with genes and values
#' @param anno.pre output of \code{\link{import.annotations}}
#' @param anno.opts 'name' for gene set names, 'syms' for gene set symbols,
#'   'info' for gene set information, 'auto' for automatically generated
#'   annotations, 'file' for manual annotations
#' @param data.pre output of \code{\link{import.database}}
#' @param info.from "annotation" or "database"
#' @param categories optionally filter categories
#' @param organisms optionally filter organisms
#' @param universe number of genes in universe
#'
#' @return overlap statistics
#' @export
#'
#' @examples \dontrun{
#' anno.raw <- read.common('path/to/anno.csv', ',', T)
#' anno.pre <- import.annotations(anno.raw, c(2, 9), 10)
#' 
#' data.raw <- read.common('path/to/data.csv', ',', T)
#' data.pre <- import.database(data.raw, c(2, 9), 10)
#' 
#' workflow('GENE1 0.1 GENE2 0.2 GENE3 0.3', anno.pre, data.pre = data.pre,
#'          universe = 10000)
#' }
workflow <- function(input.raw, anno.pre, anno.opts = "file", data.pre, info.from = "database", 
  categories = NULL, organisms = NULL, universe) {
  input <- process.input(input.raw)
  info <- if (info.from == "database") 
    data.pre$gs_info %>% dplyr::select("name", "info") else anno.pre %>% dplyr::select("name", "info")
  anno <- process.annotations(anno.pre, info, anno.opts)
  data <- process.database(data.pre, categories, organisms)
  
  calculate(input$input, anno$annos, anno$gs_annos, data$gs_genes, universe)
}
