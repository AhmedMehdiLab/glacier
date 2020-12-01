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
#' @param input output of \code{\link{process_input_text}} or
#'   \code{\link{process_input_seurat}}
#' @param stats value from \code{\link{compute}}
#' @param xpos x axis position, "top" or "bottom"
#'
#' @return ggplot2: heatmap of overlap
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_raster
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @examples
#' anno_path <- system.file("extdata", "ex_anno.csv", package = "glacier")
#' data_path <- system.file("extdata", "ex_data.csv", package = "glacier")
#' anno <- import_annotations(anno_path, ",", TRUE, c(2, 4), 5)
#' data <- import_database(data_path, ",", FALSE, c(2, 4), 0)
#'
#' input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
#' results <- compute(input, anno, data)
#' over <- plot_overlap(results$matches, "Gene Value", input, results$stats)
plot_overlap <- function(matches, value, input, stats, xpos = "bottom") {
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
  data <- expand.grid(Gene = genes, Annotation = annos, Value = NA_real_)
  for (i in seq_len(nrow(data))) {
    anno <- data$Annotation[i] %>% as.character()
    gene <- data$Gene[i] %>% as.character()

    if (gene %in% matches[[anno]]) data$Value[i] <- get_value(gene, anno)
  }
  colnames(data)[3] <- value

  # plot data
  ggplot(data, aes(.data$Gene, .data$Annotation, fill = .data[[value]])) +
    geom_raster() + if (nrow(data)) ggplot2::scale_x_discrete(position = xpos)
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
#' @param sort_y whether to sort annotations by \code{value}
#'
#' @return ggplot2: bar chart of statistics
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_col
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @examples
#' anno_path <- system.file("extdata", "ex_anno.csv", package = "glacier")
#' data_path <- system.file("extdata", "ex_data.csv", package = "glacier")
#' anno <- import_annotations(anno_path, ",", TRUE, c(2, 4), 5)
#' data <- import_database(data_path, ",", FALSE, c(2, 4), 0)
#'
#' input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
#' results <- compute(input, anno, data)
#' stat <- plot_stats(results$stats, 'Fold Enrichment', 'Adjusted P-value')
plot_stats <- function(stats, value, color, sort_y = FALSE) {
  # prepare axes
  . <- NULL
  value <- rlang::sym(value)
  color <- rlang::sym(color)

  # order Y axis
  if (nrow(stats) != 0 && sort_y) stats %<>% dplyr::arrange(dplyr::desc(!!value))
  stats$Annotation %<>% factor(., levels = .) %>% forcats::fct_rev()

  # plot data
  ggplot(stats, aes(!!value, .data$Annotation, fill = !!color)) + geom_col()
}

#' Plot calculated scores
#'
#' @param scores value from \code{\link{score_seurat}} or
#'   \code{\link{score_expr}}
#' @param x grouping variable e.g. "bin"
#' @param y plotting variable e.g. "pca"
#' @param color coloring variable e.g. "grp"
#' @param mode 'whiskers', 'box' or 'violin'
#'
#' @return ggplot2: scores
#' @export
#'
#' @examples
#' seu_path <- system.file("extdata", "ex_seurat.rds", package = "glacier")
#' seurat <- readRDS(seu_path)
#' 
#' results <- score_seurat(seurat, "grp", c("APOE", "CTSZ"))
#' scores <- plot_scores(results$scores, "seurat_clusters", "pca", "grp", "box")
plot_scores <- function(scores, x, y, color, mode) {
  if (mode == "whiskers") {
    se <- function(xs) sqrt(stats::var(xs) / length(xs))
    summary <- scores %>%
      dplyr::group_by(dplyr::across(c(!!x, !!color))) %>%
      dplyr::summarise(mean = mean(.data[[y]]), se = se(.data[[y]]))
    
    return(
      ggplot2::ggplot(summary, ggplot2::aes(.data[[x]], mean, color = .data[[color]])) +
        ggplot2::geom_point() +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se), width = 0.1)
      )
  }
  
  if (mode == "box") func <- ggplot2::geom_boxplot
  else if (mode == "violin") func <- ggplot2::geom_violin
  else stop("plot_scores: mode should be 'whiskers', 'box' or 'violin'")
  
  ggplot2::ggplot(scores, ggplot2::aes(.data[[x]], .data[[y]])) +
    func(ggplot2::aes(color = .data[[color]]))
}

#' Plot value of AUC
#'
#' @param aucs value from \code{\link{score_seurat}} or
#'   \code{\link{score_expr}}
#' @param variable value to plot e.g. "pca"
#'
#' @return ggplot2: AUC
#' @export
#'
#' @examples
#' seu_path <- system.file("extdata", "ex_seurat.rds", package = "glacier")
#' seurat <- readRDS(seu_path)
#' 
#' results <- score_seurat(seurat, "grp", c("APOE", "CTSZ"))
#' auc <- plot_auc(results$aucs, "pca")
plot_auc <- function(aucs, variable) {
  aucs$cluster <- factor(aucs$cluster, levels = aucs$cluster)
  ggplot2::ggplot(aucs, ggplot2::aes(.data$cluster, .data[[variable]])) +
    ggplot2::geom_col() + ggplot2::ylim(0, 1) + ggplot2::ylab("AUC")
}
