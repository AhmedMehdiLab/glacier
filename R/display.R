#' Plot gene-annotation matches as heatmap
#'
#' @param input value from \code{\link{process_input}}
#' @param stats value from \code{\link{calculate}}
#' @param matches value from \code{\link{calculate}}
#' @param value option for tile color as 'Gene Value' or one of '# gene sets',
#'   '#genes', '# matches', 'P-value', 'Adjusted P-value', 'Odds Ratio',
#'   'Fold Enrichment' or 'Adjusted Fold Enrichment'
#' @param value_trans optional value transformation
#'
#' @return ggplot heatmap
#' @export
#'
#' @examples
#' \dontrun{
#' anno_raw <- read.common('path/to/anno.csv', ',', T)
#' anno_pre <- import.annotations(anno_raw, c(2, 9), 10)
#' info <- anno_pre[c('name', 'info')]
#' anno <- process_annotations(anno_pre, info, c('file', 'auto'))
#'
#' data_raw <- read.common('path/to/data.csv', ',', T)
#' data_pre <- import.database(data_raw, c(2, 9), 10)
#' data <- process_database(data_pre, 'Not assigned', 'Not assigned')
#'
#' input <- process_input('GENE1 0.1 GENE2 0.2 GENE3 0.3')
#' stats <- calculate(input, anno$annos, anno$gs_annos, data$gs_genes, 10000)
#' display_overlap(input$gene, stats$stats, stats$matches, 'Gene Value')
#' }
display_overlap <- function(input, stats, matches, value,
                            value_trans = "identity") {
  . <- NULL
  genes <- input$gene %>% factor(., levels = .)
  annos <- stats$Annotation %>% factor(., levels = .) %>% forcats::fct_rev()

  # helper function to find gene values
  get_value <- function(gene, anno) {
    if (value == "Gene Value") {
      input %>% dplyr::filter(.data$gene == gene) %>% dplyr::pull(.data$gene)
    } else {
      stats %>%
        dplyr::filter(.data$Annotation == anno) %>%
        dplyr::select(!!value) %>%
        dplyr::pull
    }
  }

  # construct data grid
  data <- expand.grid(Gene = genes, Anno = annos, Value = NA_real_)
  for (i in seq_len(nrow(data))) {
    anno <- data$Anno[i] %>% as.character()
    gene <- data$Gene[i] %>% as.character()

    if (gene %in% matches[[anno]]) {
      data$Value[i] <- get_value(gene, anno)
    }
  }

  # plot data
  ggplot2::ggplot(data,
                  ggplot2::aes(.data$Gene, .data$Anno, fill = .data$Value)) +
    ggplot2::geom_raster() +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::scale_fill_gradientn(na.value = "transparent", trans = value_trans,
      colours = grDevices::hcl.colors(3, palette = "Blue-Red 2")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0)) +
    if (nrow(data) != 0) ggplot2::scale_x_discrete(position = "top")
}

#' Plot overlap statistics on a bar graph
#'
#' @param stats value from \code{\link{calculate}}
#' @param value option for bar length as one of '# gene sets', '#genes',
#'   '# matches', 'P-value', 'Adjusted P-value', 'Odds Ratio', 'Fold Enrichment'
#'   or 'Adjusted Fold Enrichment'
#' @param color option for bar color as one of '# gene sets', '#genes',
#'   '# matches', 'P-value', 'Adjusted P-value', 'Odds Ratio', 'Fold Enrichment'
#'   or 'Adjusted Fold Enrichment'
#' @param value_trans optional value transformation
#' @param color_trans optional color transformation
#' @param sort whether to sort annotations alphabetically
#'
#' @return ggplot bar graph
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' anno_raw <- read.common('path/to/anno.csv', ',', T)
#' anno_pre <- import.annotations(anno_raw, c(2, 9), 10)
#' info <- anno_pre[c('name', 'info')]
#' anno <- process_annotations(anno_pre, info, c('file', 'auto'))
#'
#' data_raw <- read.common('path/to/data.csv', ',', T)
#' data_pre <- import.database(data_raw, c(2, 9), 10)
#' data <- process_database(data_pre, 'Not assigned', 'Not assigned')
#'
#' input <- process_input('GENE1 0.1 GENE2 0.2 GENE3 0.3')
#' stats <- calculate(input, anno$annos, anno$gs_annos, data$gs_genes, 10000)
#' display_stats(stats$stats, 'Fold Enrichment', 'Adjusted P-value')
#' }
display_stats <- function(stats, value, color, value_trans = "identity",
                          color_trans = "identity", sort = F) {
  # prepare axes
  . <- NULL
  value <- rlang::sym(value)
  color <- rlang::sym(color)

  # order Y axis
  if (sort) stats %<>% dplyr::arrange(dplyr::desc(!!value))
  stats$Annotation %<>% factor(., levels = .) %>% forcats::fct_rev()

  # plot data
  ggplot2::ggplot(stats,
                  ggplot2::aes(!!value, .data$Annotation, fill = !!color)) +
    ggplot2::geom_col(colour = "black") +
    ggplot2::labs(y = NULL) + ggplot2::scale_x_continuous(trans = value_trans) +
    ggplot2::scale_fill_gradient(low = "black", high = "white",
                                 trans = color_trans) +
    ggplot2::theme_minimal()
}
