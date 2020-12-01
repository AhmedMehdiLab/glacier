context("Additional tests")
library(magrittr)
library(tibble)

anno_1 <- system.file("extdata", "ex_anno.csv", package = "glacier")
data_1 <- system.file("extdata", "ex_data.csv", package = "glacier")

test_that("function 'explore_annotation' supports gene filtering", {
  anno <- import_annotations(anno_1, ",", T, c(2, 4), 5)
  data <- import_database(data_1, ",", F, c(2, 4), 0)
  input <- process_input_text("FCN1, FTL")

  info <- anno[c("name", "info")]
  anno_proc <- glacier:::process_annotations(anno, info, "file")
  data_proc <- glacier:::process_database(data)

  expect_equal(
    explore_annotation("Carcinogen", anno_proc$gs_annos,
                       data_proc$gs_genes)$genes,
    c("FTL", "SOX2", "FCN1")
  )
  expect_equal(
    explore_annotation("Carcinogen", anno_proc$gs_annos, data_proc$gs_genes,
                       input$gene)$genes,
    c("FTL", "FCN1")
  )
})

test_that("function 'compute' does not crash", {
  anno <- import_annotations(anno_1, ",", T, c(2, 4), 5)
  data <- import_database(data_1, ",", F, c(2, 4), 0)
  input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")

  expect_error(stats <- compute(input, anno, data), NA)
  expect_error(stats <- compute(input, anno, data, 10000), NA)
  expect_error(stats <- compute(input, anno, data, save = tempfile()), NA)
})

test_that("function 'plot_overlap' does not crash", {
  anno <- import_annotations(anno_1, ",", T, c(2, 4), 5)
  data <- import_database(data_1, ",", F, c(2, 4), 0)
  input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
  res <- compute(input, anno, data, 10000)

  expect_error(plot_overlap(res$matches, "Gene Value", input, res$stats), NA)
  expect_error(plot_overlap(res$matches, "P-value", input, res$stats), NA)
})

test_that("function 'plot_stats' does not crash", {
  anno <- import_annotations(anno_1, ",", T, c(2, 4), 5)
  data <- import_database(data_1, ",", F, c(2, 4), 0)
  input <- process_input_text("FCN1 0.1 FTL 0.8 CLU 0.05")
  stats <- compute(input, anno, data, 10000)

  expect_error(plot_stats(stats$stats, "# genes", "P-value"), NA)
  expect_error(plot_stats(stats$stats, "# genes", "P-value", sort_y = T), NA)
})

test_that("function 'plot_scores' does not crash", {
  seu_path <- system.file("extdata", "ex_seurat.rds", package = "glacier")
  exp_path <- system.file("extdata", "ex_expr.rds", package = "glacier")
  seurat <- readRDS(seu_path)
  exp <- readRDS(exp_path)
  
  seu_results <- score_seurat(seurat, "grp", c("APOE", "CTSZ"))
  exp_results <- score_expr(exp, "grp", c("APOE", "CTSZ"))
  
  expect_error(plot_scores(seu_results$scores, "seurat_clusters", "pca", "grp", "whiskers"), NA)
  expect_error(plot_scores(exp_results$scores, "bin", "pca", "grp", "whiskers"), NA)
  
  expect_error(plot_scores(seu_results$scores, "seurat_clusters", "pca", "grp", "box"), NA)
  expect_error(plot_scores(exp_results$scores, "bin", "pca", "grp", "box"), NA)
  
  expect_error(plot_scores(seu_results$scores, "seurat_clusters", "pca", "grp", "violin"), NA)
  expect_error(plot_scores(exp_results$scores, "bin", "pca", "grp", "violin"), NA)
})

test_that("function 'plot_auc' does not crash", {
  seu_path <- system.file("extdata", "ex_seurat.rds", package = "glacier")
  exp_path <- system.file("extdata", "ex_expr.rds", package = "glacier")
  seurat <- readRDS(seu_path)
  exp <- readRDS(exp_path)
  
  seu_results <- score_seurat(seurat, "grp", c("APOE", "CTSZ"))
  exp_results <- score_expr(exp, "grp", c("APOE", "CTSZ"))
  
  expect_error(plot_auc(seu_results$aucs, "pca"), NA)
  expect_error(plot_auc(exp_results$aucs, "pca"), NA)
})
