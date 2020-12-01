#' Calculate expression scores
#'
#' @param expression expression data
#' @param group grouping variable from expression data
#' @param genes genes of interest
#'
#' @return expression pca, pls and plsda scores with AUC
#' @export
#'
#' @examples
#' expr_path <- system.file("extdata", "ex_expr.rds", package = "glacier")
#' expr <- readRDS(expr_path)
#'
#' results <- score_expr(expr, "grp", c("APOE", "CTSZ"))
score_expr <- function(expression, group, genes) {
  uses(c("mixOmics", "pROC"), stop, "'mixOmics' and 'pROC' are required")
  expression$lvl <- expression[[group]] %>% factor() %>% as.integer()
  data <- expression[, genes]
  
  expression$pca <- mixOmics::pca(t(data))$loadings$X[, 1]
  expression$pls <- mixOmics::pls(data, expression$lvl,
                                  mode = "regression")$variates$X[, 1]
  expression$plsda <- mixOmics::plsda(data, expression$lvl,
                                      mode = "regression")$variates$X[, 1]
  
  aucs <- tibble::tibble(
    cluster = NA,
    pca = as.numeric(pROC::multiclass.roc(lvl ~ pca, expression)$auc),
    pls = as.numeric(pROC::multiclass.roc(lvl ~ pls, expression)$auc),
    plsda = as.numeric(pROC::multiclass.roc(lvl ~ plsda, expression)$auc)
  )
  
  return(list(scores = expression, aucs = aucs))
}

#' Calculate Seurat scores
#'
#' @param seurat Seurat object
#' @param group grouping variable from Seurat data
#' @param genes genes of interest
#'
#' @return Seurat pca, pls and plsda scores with AUC
#' @export
#'
#' @importFrom magrittr %>%
#' @examples
#' seu_path <- system.file("extdata", "ex_seurat.rds", package = "glacier")
#' seurat <- readRDS(seu_path)
#' 
#' results <- score_seurat(seurat, "grp", c("APOE", "CTSZ"))
score_seurat <- function(seurat, group, genes) {
  uses(c("mixOmics", "pROC", "Seurat"), stop,
       "'mixOmics', 'pROC' and 'Seurat' are required")
  scores <- tibble::tibble()
  aucs <- tibble::tibble(cluster = character(), pca = numeric(),
                         pls = numeric(), plsda = numeric())
  
  # score each cluster
  for (i in levels(seurat$seurat_clusters)) {
    cluster <- subset(seurat, idents = i)
    cl_data <- cluster@assays$integrated@data %>% data.frame() %>% t() %>%
      data.frame() %>% dplyr::select(dplyr::all_of(genes))
    cl_meta <- cluster@meta.data
    cl_meta$grp_lvl <- cl_meta[[group]] %>% factor() %>% as.integer()
    
    # calculate scores
    cl_meta$pca <-
      tryCatch(mixOmics::pca(t(cl_data))$loadings$X[, 1],
               error = function(e) NA)
    cl_meta$pls <- 
      tryCatch(mixOmics::pls(cl_data, cl_meta$grp_lvl, mode = "regression")$variates$X[, 1],
               error = function(e) NA)
    cl_meta$plsda <-
      tryCatch(mixOmics::plsda(cl_data, cl_meta$grp_lvl)$variates$X[, 1],
               error = function(e) NA)
    
    # calculate AUC
    scores <- dplyr::bind_rows(scores, cl_meta)
    aucs <- aucs %>% dplyr::add_row(
      cluster = i,
      pca = tryCatch(pROC::multiclass.roc(grp_lvl ~ pca, cl_meta)$auc %>%
                       as.numeric(), error = function(e) NA),
      pls = tryCatch(pROC::multiclass.roc(grp_lvl ~ pls, cl_meta)$auc %>%
                       as.numeric(), error = function(e) NA),
      plsda = tryCatch(pROC::multiclass.roc(grp_lvl ~ plsda, cl_meta)$auc %>%
                         as.numeric(), error = function(e) NA)
    )
  }
  
  return(list(scores = scores, aucs = aucs))
}

