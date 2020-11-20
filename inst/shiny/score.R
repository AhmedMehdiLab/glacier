library(mixOmics)
library(pROC)

score_seurat <- function(seurat, genes) {
  scores <- tibble()
  rocs <- tibble(cluster = character(), pca = numeric(), pls = numeric(), plsda = numeric())
  
  for (i in levels(seurat$seurat_clusters)) {
    cluster <- subset(seurat, idents = i)
    cl_data <- cluster@assays$integrated@data %>% data.frame
    cl_meta <- cluster@meta.data
    cl_meta$grp_lvl <- cl_meta$grp %>% factor %>% as.integer
    
    cl_meta$pca <- pca(cl_data[genes, ], ncomp = 3)$loadings$X[, 1]
    cl_meta$pls <- pls(t(cl_data[genes, ]), cl_meta$grp_lvl, ncomp = 10, mode = "regression")$variates$X[, 1]
    cl_meta$plsda <- plsda(t(cl_data[genes, ]), cl_meta$grp_lvl, ncomp = 2)$variates$X[, 1]
    
    scores <- bind_rows(scores, cl_meta)
    rocs <- rocs %>%
      add_row(cluster = i,
              pca = multiclass.roc(grp_lvl ~ pca, data = cl_meta)$auc %>% as.numeric,
              pls = multiclass.roc(grp_lvl ~ pls, data = cl_meta)$auc %>% as.numeric,
              plsda = multiclass.roc(grp_lvl ~ plsda, data = cl_meta)$auc %>% as.numeric)
  }
  
  return(list(scores = scores, rocs = rocs))
}

score_exp <- function(expression, genes) {
  rocs <- tibble()
  
  expression$grp_lvl <- expression$grp %>% factor %>% as.integer
  expression$pca <- pca(t(expression[, genes]), ncomp = 2)$loadings$X[, 1]
  expression$pls <- pls(expression[, genes], expression$grp_lvl, ncomp = 2, mode = "regression")$variates$X[, 1]
  expression$plsda <- plsda(expression[, genes], expression$grp_lvl, ncomp = 2, mode = "regression")$variates$X[, 1]
  
  rocs <- tibble(pca = multiclass.roc(grp_lvl ~ pca, data = expression)$auc %>% as.numeric,
                 pls = multiclass.roc(grp_lvl ~ pls, data = expression)$auc %>% as.numeric,
                 plsda = multiclass.roc(grp_lvl ~ plsda, data = expression)$auc %>% as.numeric)
  
  return(list(scores = expression, rocs = rocs))
}

score_summary <- function(scores, variable, groups) {
  se <- function(x) sqrt(var(x, na.rm = T) / length(x))
  
  scores %>%
    group_by(across(!!groups)) %>%
    summarise(mean = mean(.data[[variable]], na.rm = T),
              se = se(.data[[variable]]))
}

show_summary <- function(summary, annotation, xaxis) {
  ggplot(summary, aes(x = .data[[xaxis]], y = mean, color = grp)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +
    xlab("Seurat clusters") +
    ylab(str_c(annotation, "scores", sep = " ")) +
    theme_classic() +
    theme(legend.position = "top")
}

show_rocs <- function(rocs, variable) {
  rocs$cluster <- rocs$cluster %>% factor(., levels = .)
  
  ggplot(rocs, aes(x = cluster, y = .data[[variable]])) +
    geom_col() +
    ylim(0, 1) +
    theme_classic()
}
