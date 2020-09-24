#' Plot gene-annotation matches as heatmap
#'
#' Heatmap will plot genes on the X-axis against annotations on the Y-axis.
#' Matches will be shown with colored squares, and non-matches will be
#' transparent. If the \code{value} is "Gene Value", "Odds Ratio", "Fold
#' Enrichment" or "Adjusted Fold Enrichment", if the value is NaN or an
#' infinity, these will also be transparent.
#'
#' @param matches value from \code{\link{compute}}
#' @param value \code{"Gene Value"} for values from \code{input} or one of
#'   \code{"#gene sets"}, \code{"# genes"}, \code{"# matches"},
#'   \code{"P-value"}, \code{"Adjusted P-value"}, \code{"Odds Ratio"},
#'   \code{"Fold Enrichment"} or \code{"Adjusted Fold Enrichment"} for values
#'   from \code{stats}
#' @param input output of \code{\link{process_input}}
#' @param stats value from \code{\link{compute}}
#' @param val_trans optional: value transformation, see \code{trans} argument in
#'   \code{\link[ggplot2]{scale_continuous}}; default \code{"identity"}
#'
#' @return ggplot2: heatmap of overlap
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples
#' anno_path <- system.file("extdata", "ex_anno.csv", package = "glacier")
#' anno <- import_annotations(anno_path, ",", FALSE, c(2, 10), 11)
#' data_path <- system.file("extdata", "ex_data.csv", package = "glacier")
#' data <- import_database(data_path, ",", FALSE, c(2, 4), 0)
#'
#' input <- process_input('CYP1A1 0.2 CYP1B1 NQO1 0.3 SODD 9.0')
#' results <- compute(input, anno, data, 10000)
#' over <- plot_overlap(results$matches, "Gene Value", input, results$stats)
plot_overlap <- function(matches, value, input, stats, val_trans = "identity") {
  . <- NULL
  genes <- input$gene %>% factor(., levels = .)
  annos <- stats$Annotation %>% factor(., levels = .) %>% forcats::fct_rev()

  # helper function to find gene values
  get_value <- function(g, a) {
    if (value == "Gene Value") input %>%
      dplyr::filter(.data$gene == g) %>%
      dplyr::pull(.data$value)
    else stats %>%
      dplyr::filter(.data$Annotation == a) %>%
      dplyr::pull(!!value)
  }

  # construct data grid
  data <- expand.grid(Gene = genes, Anno = annos, Value = NA_real_)
  for (i in seq_len(nrow(data))) {
    anno <- data$Anno[i] %>% as.character()
    gene <- data$Gene[i] %>% as.character()

    if (gene %in% matches[[anno]]) data$Value[i] <- get_value(gene, anno)
  }

  # plot data
  ggplot2::ggplot(
    data, ggplot2::aes(.data$Gene, .data$Anno, fill = .data$Value)
  ) +
    ggplot2::geom_raster() +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::scale_fill_gradientn(
      na.value = "transparent", trans = val_trans,
      colours = grDevices::hcl.colors(3, palette = "Blue-Red 2")
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0)) +
    if (nrow(data) != 0) ggplot2::scale_x_discrete(position = "top")
}

#' Plot overlap statistics as a bar graph
#'
#' @param stats value from \code{\link{compute}}
#' @param value \code{"#gene sets"}, \code{"# genes"}, \code{"# matches"},
#'   \code{"P-value"}, \code{"Adjusted P-value"}, \code{"Odds Ratio"},
#'   \code{"Fold Enrichment"} or \code{"Adjusted Fold Enrichment"}
#' @param color \code{"#gene sets"}, \code{"# genes"}, \code{"# matches"},
#'   \code{"P-value"}, \code{"Adjusted P-value"}, \code{"Odds Ratio"},
#'   \code{"Fold Enrichment"} or \code{"Adjusted Fold Enrichment"}
#' @param val_trans optional: value transformation, see \code{trans} argument in
#'   \code{\link[ggplot2]{scale_continuous}}; default \code{"identity"}
#' @param col_trans optional: color transformation, see \code{trans} argument in
#'   \code{\link[ggplot2]{scale_continuous}}; default \code{"identity"}
#' @param sort_y whether to sort annotations by \code{value}
#'
#' @return ggplot2: bar chart of statistics
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @examples
#' anno_path <- system.file("extdata", "ex_anno.csv", package = "glacier")
#' anno <- import_annotations(anno_path, ",", FALSE, c(2, 10), 11)
#' data_path <- system.file("extdata", "ex_data.csv", package = "glacier")
#' data <- import_database(data_path, ",", FALSE, c(2, 4), 0)
#'
#' input <- process_input('CYP1A1 0.2 CYP1B1 NQO1 0.3 SODD 9.0')
#' results <- compute(input, anno, data, 10000)
#' stat <- plot_stats(results$stats, 'Fold Enrichment', 'Adjusted P-value')
plot_stats <- function(stats, value, color, val_trans = "identity",
                       col_trans = "identity", sort_y = FALSE) {
  # prepare axes
  . <- NULL
  value <- rlang::sym(value)
  color <- rlang::sym(color)

  # order Y axis
  if (sort_y) stats %<>% dplyr::arrange(dplyr::desc(!!value))
  stats$Annotation %<>% factor(., levels = .) %>% forcats::fct_rev()

  # plot data
  ggplot2::ggplot(
    stats, ggplot2::aes(!!value, .data$Annotation, fill = !!color)
  ) +
    ggplot2::geom_col(colour = "black") +
    ggplot2::labs(y = NULL) +
    ggplot2::scale_x_continuous(trans = val_trans) +
    ggplot2::scale_fill_gradient(low = "black", high = "white",
                                 trans = col_trans) +
    ggplot2::theme_minimal()
}
