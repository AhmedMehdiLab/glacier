#' Plot gene-annotation matches as heatmap
#'
#' @param input value from \code{\link{process.input}}
#' @param stats value from \code{\link{calculate}}
#' @param matches value from \code{\link{calculate}}
#' @param value option for tile color
#' @param value.trans optional value transformation
#'
#' @return ggplot heatmap
#' @export
#'
#' @examples
#' \dontrun{
#' display.overlap(input, stats, matches, 'Gene Value')
#' }
display.overlap <- function(input, stats, matches, value, value.trans = "identity") {
  genes <- input$gene %>% factor(.data, levels = .data)
  annos <- stats$Annotation %>% factor(.data, levels = .data) %>% forcats::fct_rev
  
  # helper function to find gene values
  get_value <- function(gene, anno) {
    if (value == "Gene Value") {
      input %>% dplyr::filter(.data$gene == gene) %>% dplyr::pull(.data$gene)
    } else {
      stats %>% dplyr::filter(.data$Annotation == anno) %>% dplyr::select(!!value) %>% 
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
  ggplot2::ggplot(data, ggplot2::aes(.data$Gene, .data$Anno, fill = .data$Value)) + 
    ggplot2::geom_raster() + ggplot2::labs(x = NULL, y = NULL) + ggplot2::scale_fill_gradientn(colours = grDevices::hcl.colors(3, 
    palette = "Blue-Red 2"), na.value = "transparent", trans = value.trans) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0)) + 
    if (nrow(data) != 0) {
      ggplot2::scale_x_discrete(position = "top")
    }
}

#' Plot overlap statistics on a bar graph
#'
#' @param stats value from \code{\link{calculate}}
#' @param value option for bar length
#' @param color option for bar color
#' @param value.trans optional value transformation
#' @param color.trans optional color transformation
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
#' display.stats(stats, value, color)
#' }
display.stats <- function(stats, value, color, value.trans = "identity", color.trans = "identity", 
  sort = F) {
  # prepare axes
  value <- rlang::sym(value)
  color <- rlang::sym(color)
  
  # order Y axis
  if (sort) {
    stats %<>% dplyr::arrange(dplyr::desc(!!value))
  }
  stats$Annotation %<>% factor(.data, levels = .data) %>% forcats::fct_rev
  
  # plot data
  ggplot2::ggplot(stats, ggplot2::aes(!!value, .data$Annotation, fill = !!color)) + 
    ggplot2::geom_col(colour = "black") + ggplot2::labs(y = NULL) + ggplot2::scale_x_continuous(trans = value.trans) + 
    ggplot2::scale_fill_gradient(low = "black", high = "white", trans = color.trans) + 
    ggplot2::theme_minimal()
}
