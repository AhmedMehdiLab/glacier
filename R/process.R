#' Process annotations
#'
#' @param anno.pre output of \code{\link{import.annotations}}
#' @param info tibble of gene set name 'name' and information 'info'
#' @param options 'name' for gene set names, 'syms' for gene set symbols,
#'   'info' for gene set information, 'auto' for automatically generated
#'   annotations, 'file' for manual annotations
#'
#' @return tibble mapping gene sets to annotations and vector of annotations
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' anno.raw <- read.common('path/to/anno.csv', ',', T)
#' anno.pre <- import.annotations(anno.raw, c(2, 9), 10)
#' info.pre <- anno.pre[c('name', 'info')]
#' process.annotations(anno.pre, info.pre, c('file', 'auto'))
#' }
process.annotations <- function(anno.pre, info, options) {
  anno <- anno.pre %>% dplyr::select(.data$name)
  info <- anno %>% dplyr::left_join(info, by = "name") %>% dplyr::pull(info)
  
  # constants
  SPLIT_1 <- stringr::str_c(sep = "|", " in comparison between ", " in comparison of ", 
    " in ", " during ", " after ", " upon ")
  SPLIT_2 <- stringr::str_c(sep = "|", " versus ", " vs ", " before ", " after ", 
    " compared to ")
  FORM_LHS <- stringr::str_glue("(?:{SPLIT_1})(.*)(?:{SPLIT_2})")
  FORM_RHS <- stringr::str_glue("(?:{SPLIT_2})(.*?)\\.?$")
  
  # extract annotations
  if ("name" %in% options) {
    anno$anno.name <- anno$name
  }
  if ("syms" %in% options) {
    anno$anno.syms <- stringr::str_match(anno$name, "_(.*)")[, 2]
  }
  if ("info" %in% options) {
    anno$anno.info <- info
  }
  if ("auto" %in% options) {
    anno$anno.auto <- stringr::str_detect(anno$name, "_UP") %>% ifelse(stringr::str_match(info, 
      FORM_LHS)[, 2], stringr::str_match(info, FORM_RHS)[, 2])
  }
  if ("file" %in% options) {
    anno %<>% tibble::add_column(dplyr::select(anno.pre, dplyr::starts_with("anno.")))
  }
  
  # generate annotation list
  annos <- anno %>% dplyr::select(dplyr::starts_with("anno.")) %>% unlist(use.names = F) %>% unique()
  list(gs_annos = anno, annos = annos[!is.na(annos) & annos != ""])
}

#' Process database
#'
#' @param data.pre output of \code{\link{import.database}} or
#'   \code{\link{import.msigdb_xml}}
#' @param categories optionally filter categories
#' @param organisms optionally filter organisms
#'
#' @return list of gene sets to genes, gene set information and vector of genes
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' data.raw <- read.common('path/to/data.csv', ',', T)
#' data.pre <- import.database(data.raw, c(2, 9), 10)
#' process.database(data.pre, 'Not assigned', 'Not assigned')
#' }
process.database <- function(data.pre, categories = NULL, organisms = NULL) {
  info <- data.pre$gs_info %>% dplyr::filter(.data$category %in% categories, .data$organism %in% 
    organisms)
  sets <- data.pre$gs_genes[info$name]
  gene <- sets %>% unlist(use.names = F) %>% unique()
  
  list(gs_genes = sets, gs_info = info, genes = gene)
}

#' Pre-process gene string input
#'
#' Removes duplicate genes. If duplicate values are found, the first value will
#' be used.
#'
#' @param input.raw string input with genes and values
#'
#' @return pre-processed input
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' process.input.pre('GENE1 GENE2 GENE3')
#' process.input.pre('GENE1 0.1 GENE2 0.2 GENE3 0.3')
#' }
process.input.pre <- function(input.raw) {
  tokens <- input.raw %>% stringr::str_split("[ \t\r\n,;]") %>% unlist() %>% purrr::discard(~. == 
    "")
  values <- suppressWarnings(as.numeric(tokens))
  
  # process
  tokens[!is.na(values)] <- NA
  values <- values[-1] %>% c(NA)
  
  # store results
  genes <- tibble::tibble(gene = tokens, value = values) %>% tidyr::drop_na(.data$gene) %>% 
    dplyr::distinct(.data$gene, .keep_all = T)
}

#' Post-process gene string input
#'
#' @param input.pre output of \code{\link{process.input.pre}}
#'
#' @return gene names, values and bitmap index of missing values
#'
#' @examples
#' \dontrun{
#' input.pre <- process.input.pre('GENE1 GENE2 GENE3')
#' process.input.post(input.pre)
#' }
process.input.post <- function(input.pre) {
  missing <- is.na(input.pre$value)
  input.pre$value[missing] <- 0
  
  list(input = input.pre, missing = missing)
}

#' Process gene string input
#'
#' @param input.raw string input with genes and values
#'
#' @return gene names, values and bitmap index of missing values
#' @export
#'
#' @examples
#' process.input('GENE1 GENE2 GENE3')
#' process.input('GENE1 0.1 GENE2 0.2 GENE3 0.3')
process.input <- function(input.raw) {
  input.raw %>% process.input.pre() %>% process.input.post()
}

#' Map annotation to gene sets and genes
#'
#' @param annotation annotation to map
#' @param gs_annos value from \code{\link{process.annotations}}
#' @param gs_genes value from \code{\link{process.database}}
#'
#' @return vector of associated gene sets and vector of associated genes
#' @export
#'
#' @examples
#' \dontrun{
#' anno.raw <- read.common('path/to/anno.csv', ',', T)
#' anno.pre <- import.annotations(anno.raw, c(2, 9), 10)
#' info.pre <- anno.pre[c('name', 'info')]
#' anno <- process.annotations(anno.pre, info.pre, c('file', 'auto'))
#' 
#' data.raw <- read.common('path/to/data.csv', ',', T)
#' data.pre <- import.database(data.raw, c(2, 9), 10)
#' data <- process.database(data.pre, 'Not assigned', 'Not assigned')
#' 
#' annotation.map('annotation', anno$gs_annos, data$gs_genes)
#' }
annotation.map <- function(annotation, gs_annos, gs_genes) {
  index <- (gs_annos == annotation) %>% rowSums(na.rm = T) %>% as.logical()
  names <- gs_annos$name[index]
  genes <- gs_genes[names] %>% unlist(use.names = F) %>% unique()
  
  list(gs_names = names, genes = genes)
}

#' Pre-calculate overlap statistics
#'
#' @param input output of \code{\link{process.input}}
#' @param annos value from \code{\link{process.annotations}}
#' @param gs_annos value from \code{\link{process.annotations}}
#' @param gs_genes value from \code{\link{process.database}}
#'
#' @return partially-computed statistics and list of annotation matches
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' anno.raw <- read.common('path/to/anno.csv', ',', T)
#' anno.pre <- import.annotations(anno.raw, c(2, 9), 10)
#' info.pre <- anno.pre[c('name', 'info')]
#' anno <- process.annotations(anno.pre, info.pre, c('file', 'auto'))
#' 
#' data.raw <- read.common('path/to/data.csv', ',', T)
#' data.pre <- import.database(data.raw, c(2, 9), 10)
#' data <- process.database(data.pre, 'Not assigned', 'Not assigned')
#' 
#' input <- process.input('GENE1 0.1 GENE2 0.2 GENE3 0.3')
#' calculate.pre(input, anno$annos, anno$gs_annos, data$gs_genes)
#' }
calculate.pre <- function(input, annos, gs_annos, gs_genes) {
  stat <- tibble::tibble(name = annos, n_sets = 0L, n_gene = 0L, n_hits = 0L, g_hits = "")
  hits <- list()
  
  # iterate over annotations
  for (i in seq_along(annos)) {
    # get related genes and find overlap
    overlap <- annotation.map(annos[i], gs_annos, gs_genes)
    matches <- intersect(overlap$gene, input$gene)
    
    # store information
    stat$n_sets[i] <- length(overlap$sets)
    stat$n_gene[i] <- length(overlap$gene)
    stat$n_hits[i] <- length(matches)
    stat$g_hits[i] <- matches %>% stringr::str_c(collapse = ", ")
    hits[[annos[i]]] <- matches
  }
  
  list(stats.pre = stat, matches = hits)
}

#' Finish statistical calculations
#'
#' @param stats.pre value from \code{\link{calculate.pre}}
#' @param input.count number of genes in input
#' @param universe number of genes in universe
#'
#' @return overlap statistics
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' anno.raw <- read.common('path/to/anno.csv', ',', T)
#' anno.pre <- import.annotations(anno.raw, c(2, 9), 10)
#' info.pre <- anno.pre[c('name', 'info')]
#' anno <- process.annotations(anno.pre, info.pre, c('file', 'auto'))
#' 
#' data.raw <- read.common('path/to/data.csv', ',', T)
#' data.pre <- import.database(data.raw, c(2, 9), 10)
#' data <- process.database(data.pre, 'Not assigned', 'Not assigned')
#' 
#' input <- process.input('GENE1 0.1 GENE2 0.2 GENE3 0.3')
#' calc.pre <- calculate.pre(input, anno$annos, anno$gs_annos, data$gs_genes)
#' calculate.post(calc.pre$stats.pre, nrow(input), 10000)
#' }
calculate.post <- function(stats.pre, input.count, universe) {
  stats.pre %<>% tibble::add_column(pvalue = 0, odds_r = 0)
  
  # calculate statistics
  for (i in seq_len(nrow(stats.pre))) {
    # fisher's exact test: [1, ] belongs to annotation [, 1] entered in list
    data <- matrix(nrow = 2, ncol = 2)
    data[1, 1] <- stats.pre$n_hits[i]
    data[1, 2] <- stats.pre$n_gene[i] - data[1, 1]
    data[2, 1] <- input.count - data[1, 1]
    data[2, 2] <- universe - data[1, 1] - data[1, 2] - data[2, 1]
    
    # assign statistics
    test <- data %>% stats::fisher.test(alternative = "greater")
    stats.pre$pvalue[i] <- test$p.value
    stats.pre$odds_r[i] <- test$estimate[[1]]
  }
  
  # post-processing
  stats.pre$enrich <- (stats.pre$n_hits/input.count)/(stats.pre$n_gene/universe)
  stats.pre$adj_pv <- stats.pre$pvalue %>% stats::p.adjust(method = "fdr")
  stats.pre$adj_fe <- stats.pre$enrich/-log(stats.pre$adj_pv)
  
  stats.pre %>% dplyr::select(Annotation = .data$name, `# gene sets` = .data$n_sets, 
    `# genes` = .data$n_gene, `# matches` = .data$n_hits, `P-value` = .data$pvalue, 
    `Adjusted P-value` = .data$adj_pv, `Odds Ratio` = .data$odds_r, `Fold Enrichment` = .data$enrich, 
    `Adjusted Fold Enrichment` = .data$adj_fe, Matches = .data$g_hits)
}

#' Calculate overlap statistics
#'
#' @param input value from \code{\link{process.input}}
#' @param annos value from \code{\link{process.annotations}}
#' @param gs_annos value from \code{\link{process.annotations}}
#' @param gs_genes value from \code{\link{process.database}}
#' @param universe number of genes in universe
#'
#' @return tibble of overlap statistics and list of annotation matches
#' @export
#'
#' @examples
#' \dontrun{
#' anno.raw <- read.common('path/to/anno.csv', ',', T)
#' anno.pre <- import.annotations(anno.raw, c(2, 9), 10)
#' info.pre <- anno.pre[c('name', 'info')]
#' anno <- process.annotations(anno.pre, info.pre, c('file', 'auto'))
#' 
#' data.raw <- read.common('path/to/data.csv', ',', T)
#' data.pre <- import.database(data.raw, c(2, 9), 10)
#' data <- process.database(data.pre, 'Not assigned', 'Not assigned')
#' 
#' input <- process.input('GENE1 0.1 GENE2 0.2 GENE3 0.3')
#' calculate(input$input, anno$annos, anno$gs_annos, data$gs_genes, 10000)
#' }
calculate <- function(input, annos, gs_annos, gs_genes, universe) {
  calc.pre <- calculate.pre(input, annos, gs_annos, gs_genes)
  calc.post <- calculate.post(calc.pre$stats.pre, nrow(input), universe)
  
  list(stats = calc.post, matches = calc.pre$matches)
}
