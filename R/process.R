# constants
split_1 <- stringr::str_c(sep = "|", " in comparison between ", " upon ",
                          " in comparison of ", " in ", " during ", " after ")
split_2 <- stringr::str_c(sep = "|", " versus ", " vs ", " before ", " after ",
                          " compared to ")
form_lhs <- stringr::str_glue("(?:{split_1})(.*)(?:{split_2})")
form_rhs <- stringr::str_glue("(?:{split_2})(.*?)\\.?$")

#' Process annotations
#'
#' Choose annotations to include in annotation list, from gene set names, gene
#' set symbols, gene set descriptions, annotations automatically extracted from
#' gene set descriptions or manual annotations.
#'
#' @param anno output of \code{\link{import_annotations}}
#' @param info tibble: "name" gene set name "info" gene set descriptions
#' @param options \code{"name"} for gene set names, \code{"syms"} for gene set
#'   symbols, \code{"info"} for gene set descriptions, \code{"auto"} for
#'   automatically generated annotations and/or \code{"file"} for manual
#'   annotations
#'
#' @return
#' \code{gs_annos} tibble: gene sets and annotations
#'
#' \code{annos} vector: annotations
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @examples \dontrun{
#' anno_path <- system.file("extdata", "ex_anno.csv", package = "glacier")
#' anno <- import_annotations(anno_path, ",", FALSE, c(2, 10), 11)
#' info <- anno[c("name", "info")]
#'
#' anno_proc <- process_annotations(anno, info, "file")
#' }
process_annotations <- function(anno, info, options) {
  anno_proc <- anno %>% dplyr::select("name")
  info <- anno_proc %>%
    dplyr::left_join(info, by = "name") %>%
    dplyr::pull("info")

  # extract annotations
  if ("name" %in% options) anno_proc$anno_name <- anno_proc$name
  if ("syms" %in% options) anno_proc$anno_syms <-
    stringr::str_match(anno_proc$name, "_(.*)")[, 2]
  if ("info" %in% options) anno_proc$anno_info <- info
  if ("auto" %in% options) anno_proc$anno_auto <-
    ifelse(
      stringr::str_detect(anno_proc$name, "_UP"),
      stringr::str_match(info, form_lhs)[, 2],
      stringr::str_match(info, form_rhs)[, 2]
    )
  if ("file" %in% options) anno_proc <- anno_proc %>%
    tibble::add_column(dplyr::select(anno, dplyr::starts_with("anno_")))

  # generate annotation list
  annos <- anno_proc %>%
    dplyr::select(dplyr::starts_with("anno_")) %>%
    unlist(use.names = F) %>%
    unique()

  list(gs_annos = anno_proc, annos = annos[!is.na(annos) & annos != ""])
}

#' Process database
#'
#' Choose database category and organism for use in analysis.
#'
#' @param data output of \code{\link{import_database}} or
#'   \code{\link{import_msigdb}}
#' @param categories optional: categories to include; default all
#' @param organisms optional: organisms to include; default all
#'
#' @return
#' \code{gs_genes} list: names: gene set names vector: genes
#'
#' \code{gs_info} tibble: gene set information
#'
#' \code{genes} vector: list of genes
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples \dontrun{
#' data_path <- system.file("extdata", "ex_data.csv", package = "glacier")
#' data <- import_database(data_path, ",", FALSE, c(2, 4), 0)
#'
#' data_proc <- process_database(data, "Not assigned", "Not assigned")
#' }
process_database <- function(data, categories = FALSE, organisms = FALSE) {
  # filter categories and organisms
  gs_info <- data$gs_info
  if (is.null(categories) || categories != FALSE)
    gs_info <- gs_info %>% dplyr::filter(.data$category %in% categories)
  if (is.null(organisms) || organisms != FALSE)
    gs_info <- gs_info %>% dplyr::filter(.data$organism %in% organisms)

  # extract gene sets and genes
  gs_genes <- data$gs_genes[gs_info$name]
  genes <- gs_genes %>% unlist(use.names = F) %>% unique()

  list(gs_genes = gs_genes, gs_info = gs_info, genes = genes)
}

#' Process text input
#'
#' Removes duplicate genes. If multiple values for the same gene are found, only
#' the first value will be kept.
#'
#' @param text character: input with genes and optionally values
#'
#' @return tibble: "gene" gene names "value" gene values
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples
#' input <- process_input_text('CYP1A1 CYP1B1 NQO1 SODD')
#' input <- process_input_text('CYP1A1 0.2 CYP1B1 NQO1 0.3 SODD 9.0')
process_input_text <- function(text) {
  tokens <- text %>%
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

#' Extract differentially expressed genes from Seurat object
#'
#' Finds differentially expressed genes, records adjusted P-value and filters
#' for values less than \code{max_p}.
#'
#' @param seurat Seurat object
#' @param clst_1 first cluster
#' @param clst_2 optional: second cluster; default all others
#' @param max_p P-value cutoff, only genes with P-values less than this will be
#'   returned
#'
#' @return tibble: "gene" gene names "value" gene values
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples \dontrun{
#' seu <- readRDS("path/to/seurat.rds")
#' input <- process_input_seurat(seu, 1)
#' }
process_input_seurat <- function(seurat, clst_1, clst_2 = NULL, max_p = 0.05) {
  if (!requireNamespace("Seurat", quietly = T))
    stop("Library 'Seurat' is required for this feature")
  if (!is.null(clst_2) && clst_1 == clst_2)
    return(tibble::tibble(gene = character(), value = numeric()))

  seurat %>%
    Seurat::FindMarkers(ident.1 = clst_1, ident.2 = clst_2) %>%
    tibble::rownames_to_column("gene") %>%
    tibble::tibble() %>%
    dplyr::select("gene", value = "p_val_adj") %>%
    dplyr::filter(.data$value <= max_p)
}

#' Find an annotation's associated gene sets and genes
#'
#' @param annotation annotation to explore
#' @param gs_annos value from \code{\link{process_annotations}}
#' @param gs_genes value from \code{\link{process_database}}
#' @param genes optional: genes to match, or (default) all
#'
#' @return
#' \code{"names"} vector: gene set names
#'
#' \code{"genes"} vector: genes
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @examples \dontrun{
#' anno_path <- system.file("extdata", "ex_anno.csv", package = "glacier")
#' anno <- import_annotations(anno_path, ",", FALSE, c(2, 10), 11)
#' data_path <- system.file("extdata", "ex_data.csv", package = "glacier")
#' data <- import_database(data_path, ",", FALSE, c(2, 4), 0)
#'
#' info <- anno[c("name", "info")]
#' anno_proc <- process_annotations(anno, info, "file")
#' data_proc <- process_database(data, "Not assigned", "Not assigned")
#' anno_assoc <- explore_annotation("Carcinogen", anno_proc$gs_annos,
#'                                  data_proc$gs_genes)
#' }
explore_annotation <- function(annotation, gs_annos, gs_genes, genes = NULL) {
  . <- NULL
  index <- (gs_annos == annotation) %>% rowSums(na.rm = T) %>% as.logical()
  match <- gs_genes[gs_annos$name[index]]
  if (!is.null(genes))
    match <- match %>% purrr::map(~intersect(., genes)) %>% purrr::compact()

  names <- names(match)
  genes <- match %>% unlist(use.names = F) %>% unique()
  list(names = names, genes = genes)
}

#' Begin calculating overlap statistics
#'
#' @param input output of \code{\link{process_input_text}} or
#'   \code{\link{process_input_seurat}}
#' @param annos value from \code{\link{process_annotations}}
#' @param gs_annos value from \code{\link{process_annotations}}
#' @param gs_genes value from \code{\link{process_database}}
#'
#' @return \code{stats_pre} tibble: overlap statistics (incomplete)
#'
#' \code{matches} list: names: annotations vector: matched genes
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @examples \dontrun{
#' anno_path <- system.file("extdata", "ex_anno.csv", package = "glacier")
#' anno <- import_annotations(anno_path, ",", FALSE, c(2, 10), 11)
#' data_path <- system.file("extdata", "ex_data.csv", package = "glacier")
#' data <- import_database(data_path, ",", FALSE, c(2, 4), 0)
#'
#' info <- anno[c("name", "info")]
#' anno_proc <- process_annotations(anno, info, "file")
#' data_proc <- process_database(data, "Not assigned", "Not assigned")
#'
#' input <- process_input('CYP1A1 0.2 CYP1B1 NQO1 0.3 SODD 9.0')
#' calc_pre <- calculate_pre(input, anno_proc$annos, anno_proc$gs_annos,
#'                           data_proc$gs_genes)
#' }
calculate_pre <- function(input, annos, gs_annos, gs_genes) {
  stat <- tibble::tibble(name = annos, n_sets = 0L, n_gene = 0L, n_hits = 0L,
                         g_hits = "")
  hits <- list()

  # iterate over annotations
  for (i in seq_along(annos)) {
    # get related genes and find overlap
    overlap <- explore_annotation(annos[i], gs_annos, gs_genes)
    matches <- intersect(overlap$genes, input$gene)

    # store information
    stat$n_sets[i] <- length(overlap$names)
    stat$n_gene[i] <- length(overlap$genes)
    stat$n_hits[i] <- length(matches)
    stat$g_hits[i] <- matches %>% stringr::str_c(collapse = ", ")
    hits[[annos[i]]] <- matches
  }

  list(stats_pre = stat, matches = hits)
}

#' Finish calculating overlap statistics
#'
#' @param stats_pre value from \code{\link{calculate_pre}}
#' @param input_size number of genes in input
#' @param universe number of genes in universe
#'
#' @return tibble: overlap statistics
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples \dontrun{
#' anno_path <- system.file("extdata", "ex_anno.csv", package = "glacier")
#' anno <- import_annotations(anno_path, ",", FALSE, c(2, 10), 11)
#' data_path <- system.file("extdata", "ex_data.csv", package = "glacier")
#' data <- import_database(data_path, ",", FALSE, c(2, 4), 0)
#'
#' info <- anno[c("name", "info")]
#' anno_proc <- process_annotations(anno, info, "file")
#' data_proc <- process_database(data, "Not assigned", "Not assigned")
#'
#' input <- process_input('CYP1A1 0.2 CYP1B1 NQO1 0.3 SODD 9.0')
#' calc_pre <- calculate_pre(input, anno_proc$annos, anno_proc$gs_annos,
#'                           data_proc$gs_genes)
#' calc <- calculate_post(calc_pre$stats_pre, nrow(input), 10000)
#' }
calculate_post <- function(stats_pre, input_size, universe) {
  stat <- stats_pre %>% tibble::add_column(pvalue = 0, odds_r = 0)

  # calculate statistics
  for (i in seq_len(nrow(stat))) {
    # fisher's exact test: [1, ] belongs to annotation [, 1] entered in list
    data <- matrix(nrow = 2, ncol = 2)
    data[1, 1] <- stat$n_hits[i]
    data[1, 2] <- stat$n_gene[i] - data[1, 1]
    data[2, 1] <- input_size - data[1, 1]
    data[2, 2] <- universe - data[1, 1] - data[1, 2] - data[2, 1]

    # assign statistics
    test <- data %>% stats::fisher.test(alternative = "greater")
    stat$pvalue[i] <- test$p.value
    stat$odds_r[i] <- test$estimate[[1]]
  }

  # post-processing
  stat$enrich <- (stat$n_hits / input_size) / (stat$n_gene / universe)
  stat$adj_pv <- stat$pvalue %>% stats::p.adjust(method = "fdr")
  stat$adj_fe <- stat$enrich / -log(stat$adj_pv)

  stat %>% dplyr::select(
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
#' @param input output of \code{\link{process_input_text}} or
#'   \code{\link{process_input_seurat}}
#' @param annos value from \code{\link{process_annotations}}
#' @param gs_annos value from \code{\link{process_annotations}}
#' @param gs_genes value from \code{\link{process_database}}
#' @param universe optional: number of genes in universe; default calculate from
#'   \code{gs_genes}
#'
#' @return \code{stats} tibble: overlap statistics
#'
#' \code{matches} list: names: annotations vector: matched genes
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @examples \dontrun{
#' anno_path <- system.file("extdata", "ex_anno.csv", package = "glacier")
#' anno <- import_annotations(anno_path, ",", FALSE, c(2, 10), 11)
#' data_path <- system.file("extdata", "ex_data.csv", package = "glacier")
#' data <- import_database(data_path, ",", FALSE, c(2, 4), 0)
#'
#' info <- anno[c("name", "info")]
#' anno_proc <- process_annotations(anno, info, "file")
#' data_proc <- process_database(data, "Not assigned", "Not assigned")
#'
#' input <- process_input('CYP1A1 0.2 CYP1B1 NQO1 0.3 SODD 9.0')
#' calc <- calculate(input, anno_proc$annos, anno_proc$gs_annos,
#'                   data_proc$gs_genes, 10000)
#' }
calculate <- function(input, annos, gs_annos, gs_genes, universe) {
  if (is.null(universe)) {
    genes <- gs_genes %>% unlist(use.names = F) %>% unique()
    universe <- max(nrow(input), length(genes[!is.na(genes) & genes != ""]))
  }

  calc_pre <- calculate_pre(input, annos, gs_annos, gs_genes)
  calc <- calculate_post(calc_pre$stats_pre, nrow(input), universe)

  list(stats = calc, matches = calc_pre$matches)
}
