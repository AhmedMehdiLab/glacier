#' Process annotations
#'
#' @param anno_pre output of \code{\link{import_annotations}}
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
#' anno_raw <- read_common('path/to/anno.csv', ',', T)
#' anno_pre <- import_annotations(anno_raw, c(2, 9), 10)
#' info <- anno_pre[c('name', 'info')]
#' process_annotations(anno_pre, info, c('file', 'auto'))
#' }
process_annotations <- function(anno_pre, info, options) {
  anno <- anno_pre %>% dplyr::select(.data$name)
  info <- anno %>% dplyr::left_join(info, by = "name") %>% dplyr::pull(info)

  # constants
  split_1 <- stringr::str_c(sep = "|", " in comparison between ",
                            " in comparison of ", " in ", " during ", " after ",
                            " upon ")
  split_2 <- stringr::str_c(sep = "|", " versus ", " vs ", " before ",
                            " after ", " compared to ")
  form_lhs <- stringr::str_glue("(?:{split_1})(.*)(?:{split_2})")
  form_rhs <- stringr::str_glue("(?:{split_2})(.*?)\\.?$")

  # extract annotations
  if ("name" %in% options) {
    anno$anno_name <- anno$name
  }
  if ("syms" %in% options) {
    anno$anno_syms <- stringr::str_match(anno$name, "_(.*)")[, 2]
  }
  if ("info" %in% options) {
    anno$anno_info <- info
  }
  if ("auto" %in% options) {
    anno$anno_auto <- stringr::str_detect(anno$name, "_UP") %>%
      ifelse(
        stringr::str_match(info, form_lhs)[, 2],
        stringr::str_match(info, form_rhs)[, 2]
      )
  }
  if ("file" %in% options) {
    anno %<>%
      tibble::add_column(dplyr::select(anno_pre, dplyr::starts_with("anno_")))
  }

  # generate annotation list
  annos <- anno %>%
    dplyr::select(dplyr::starts_with("anno_")) %>%
    unlist(use.names = F) %>%
    unique()
  list(gs_annos = anno, annos = annos[!is.na(annos) & annos != ""])
}

#' Process database
#'
#' @param data_pre output of \code{\link{import_database}} or
#'   \code{\link{import_msigdb_xml}}
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
#' data_raw <- read_common('path/to/data.csv', ',', T)
#' data_pre <- import_database(data_raw, c(2, 9), 10)
#' process_database(data_pre, 'Not assigned', 'Not assigned')
#' }
process_database <- function(data_pre, categories = NULL, organisms = NULL) {
  info <- data_pre$gs_info %>%
    dplyr::filter(.data$category %in% categories, .data$organism %in% organisms)
  sets <- data_pre$gs_genes[info$name]
  gene <- sets %>% unlist(use.names = F) %>% unique()

  list(gs_genes = sets, gs_info = info, genes = gene)
}

#' Pre-process gene string input
#'
#' Removes duplicate genes. If duplicate values are found, the first value will
#' be used.
#'
#' @param input_raw string input with genes and values
#'
#' @return pre-processed input
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' process_input_pre('GENE1 GENE2 GENE3')
#' process_input_pre('GENE1 0.1 GENE2 0.2 GENE3 0.3')
#' }
process_input_pre <- function(input_raw) {
  tokens <- input_raw %>%
    stringr::str_split("[ \t\r\n,;]") %>%
    unlist() %>%
    purrr::discard(~. == "")
  values <- suppressWarnings(as.numeric(tokens))

  # process
  tokens[!is.na(values)] <- NA
  values <- values[-1] %>% c(NA)

  # store results
  genes <- tibble::tibble(gene = tokens, value = values) %>%
    tidyr::drop_na(.data$gene) %>%
    dplyr::distinct(.data$gene, .keep_all = T)
}

#' Post-process gene string input
#'
#' @param input_pre output of \code{\link{process_input_pre}}
#'
#' @return gene names, values and bitmap index of missing values
#'
#' @examples
#' \dontrun{
#' input_pre <- process_input_pre('GENE1 GENE2 GENE3')
#' process_input_post(input_pre)
#' }
process_input_post <- function(input_pre) {
  missing <- is.na(input_pre$value)
  input_pre$value[missing] <- 0

  list(input = input_pre, missing = missing)
}

#' Process gene string input
#'
#' @param input_raw string input with genes and values
#'
#' @return gene names, values and bitmap index of missing values
#' @export
#'
#' @examples
#' process_input('GENE1 GENE2 GENE3')
#' process_input('GENE1 0.1 GENE2 0.2 GENE3 0.3')
process_input <- function(input_raw) {
  input_raw %>% process_input_pre() %>% process_input_post()
}

#' Map annotation to gene sets and genes
#'
#' @param annotation annotation to map
#' @param gs_annos value from \code{\link{process_annotations}}
#' @param gs_genes value from \code{\link{process_database}}
#'
#' @return vector of associated gene sets and vector of associated genes
#' @export
#'
#' @examples
#' \dontrun{
#' anno_raw <- read_common('path/to/anno.csv', ',', T)
#' anno_pre <- import_annotations(anno_raw, c(2, 9), 10)
#' info_pre <- anno_pre[c('name', 'info')]
#' anno <- process_annotations(anno_pre, info_pre, c('file', 'auto'))
#'
#' data_raw <- read_common('path/to/data.csv', ',', T)
#' data_pre <- import_database(data_raw, c(2, 9), 10)
#' data <- process_database(data_pre, 'Not assigned', 'Not assigned')
#'
#' annotation_map('annotation', anno$gs_annos, data$gs_genes)
#' }
annotation_map <- function(annotation, gs_annos, gs_genes) {
  index <- (gs_annos == annotation) %>% rowSums(na.rm = T) %>% as.logical()
  names <- gs_annos$name[index]
  genes <- gs_genes[names] %>% unlist(use.names = F) %>% unique()

  list(gs_names = names, genes = genes)
}

#' Pre-calculate overlap statistics
#'
#' @param input output of \code{\link{process_input}}
#' @param annos value from \code{\link{process_annotations}}
#' @param gs_annos value from \code{\link{process_annotations}}
#' @param gs_genes value from \code{\link{process_database}}
#'
#' @return partially-computed statistics and list of annotation matches
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' anno_raw <- read_common('path/to/anno.csv', ',', T)
#' anno_pre <- import_annotations(anno_raw, c(2, 9), 10)
#' info_pre <- anno_pre[c('name', 'info')]
#' anno <- process_annotations(anno_pre, info_pre, c('file', 'auto'))
#'
#' data_raw <- read_common('path/to/data.csv', ',', T)
#' data_pre <- import_database(data_raw, c(2, 9), 10)
#' data <- process_database(data_pre, 'Not assigned', 'Not assigned')
#'
#' input <- process_input('GENE1 0.1 GENE2 0.2 GENE3 0.3')
#' calculate_pre(input, anno$annos, anno$gs_annos, data$gs_genes)
#' }
calculate_pre <- function(input, annos, gs_annos, gs_genes) {
  stat <- tibble::tibble(
    name = annos,
    n_sets = 0L,
    n_gene = 0L,
    n_hits = 0L,
    g_hits = ""
  )
  hits <- list()

  # iterate over annotations
  for (i in seq_along(annos)) {
    # get related genes and find overlap
    overlap <- annotation_map(annos[i], gs_annos, gs_genes)
    matches <- intersect(overlap$gene, input$gene)

    # store information
    stat$n_sets[i] <- length(overlap$sets)
    stat$n_gene[i] <- length(overlap$gene)
    stat$n_hits[i] <- length(matches)
    stat$g_hits[i] <- matches %>% stringr::str_c(collapse = ", ")
    hits[[annos[i]]] <- matches
  }

  list(stats_pre = stat, matches = hits)
}

#' Finish statistical calculations
#'
#' @param stats_pre value from \code{\link{calculate_pre}}
#' @param input_count number of genes in input
#' @param universe number of genes in universe
#'
#' @return overlap statistics
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' anno_raw <- read_common('path/to/anno.csv', ',', T)
#' anno_pre <- import_annotations(anno_raw, c(2, 9), 10)
#' info_pre <- anno_pre[c('name', 'info')]
#' anno <- process_annotations(anno_pre, info_pre, c('file', 'auto'))
#'
#' data_raw <- read_common('path/to/data.csv', ',', T)
#' data_pre <- import_database(data_raw, c(2, 9), 10)
#' data <- process_database(data_pre, 'Not assigned', 'Not assigned')
#'
#' input <- process_input('GENE1 0.1 GENE2 0.2 GENE3 0.3')
#' calc_pre <- calculate_pre(input, anno$annos, anno$gs_annos, data$gs_genes)
#' calculate_post(calc_pre$stats_pre, nrow(input), 10000)
#' }
calculate_post <- function(stats_pre, input_count, universe) {
  stats_pre %<>% tibble::add_column(pvalue = 0, odds_r = 0)

  # calculate statistics
  for (i in seq_len(nrow(stats_pre))) {
    # fisher's exact test: [1, ] belongs to annotation [, 1] entered in list
    data <- matrix(nrow = 2, ncol = 2)
    data[1, 1] <- stats_pre$n_hits[i]
    data[1, 2] <- stats_pre$n_gene[i] - data[1, 1]
    data[2, 1] <- input_count - data[1, 1]
    data[2, 2] <- universe - data[1, 1] - data[1, 2] - data[2, 1]

    # assign statistics
    test <- data %>% stats::fisher.test(alternative = "greater")
    stats_pre$pvalue[i] <- test$p.value
    stats_pre$odds_r[i] <- test$estimate[[1]]
  }

  # post-processing
  stats_pre$enrich <- (stats_pre$n_hits / input_count) /
    (stats_pre$n_gene / universe)
  stats_pre$adj_pv <- stats_pre$pvalue %>% stats::p.adjust(method = "fdr")
  stats_pre$adj_fe <- stats_pre$enrich / -log(stats_pre$adj_pv)

  stats_pre %>%
    dplyr::select(
      Annotation = .data$name,
      `# gene sets` = .data$n_sets,
      `# genes` = .data$n_gene,
      `# matches` = .data$n_hits,
      `P-value` = .data$pvalue,
      `Adjusted P-value` = .data$adj_pv,
      `Odds Ratio` = .data$odds_r,
      `Fold Enrichment` = .data$enrich,
      `Adjusted Fold Enrichment` = .data$adj_fe,
      Matches = .data$g_hits
    )
}

#' Calculate overlap statistics
#'
#' @param input value from \code{\link{process_input}}
#' @param annos value from \code{\link{process_annotations}}
#' @param gs_annos value from \code{\link{process_annotations}}
#' @param gs_genes value from \code{\link{process_database}}
#' @param universe number of genes in universe
#'
#' @return tibble of overlap statistics and list of annotation matches
#' @export
#'
#' @examples
#' \dontrun{
#' # import annotations
#' anno_raw <- read_common('path/to/anno.csv', ',', T)
#' anno_pre <- import_annotations(anno_raw, c(2, 9), 10)
#' info_pre <- anno_pre[c('name', 'info')]
#' anno <- process_annotations(anno_pre, info_pre, c('file', 'auto'))
#'
#' data_raw <- read_common('path/to/data.csv', ',', T)
#' data_pre <- import_database(data_raw, c(2, 9), 10)
#' data <- process_database(data_pre, 'Not assigned', 'Not assigned')
#'
#' input <- process_input('GENE1 0.1 GENE2 0.2 GENE3 0.3')
#' calculate(input$input, anno$annos, anno$gs_annos, data$gs_genes, 10000)
#' }
calculate <- function(input, annos, gs_annos, gs_genes, universe) {
  calc_pre <- calculate_pre(input, annos, gs_annos, gs_genes)
  calc_post <- calculate_post(calc_pre$stats_pre, nrow(input), universe)

  list(stats = calc_post, matches = calc_pre$matches)
}
