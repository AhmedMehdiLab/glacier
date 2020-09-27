context("Seurat object processing")
library(mockery)

seu_path <- system.file("extdata", "ex_seurat.rds", package = "glacier")
seurat <- readRDS(seu_path)

test_that("requires Seurat to function", {
  mockery::stub(process_input_seurat, "requireNamespace", FALSE)
  expect_error(process_input_seurat(seurat, 0)$gene,
               "Package 'Seurat' is required for this feature")
})

test_that("differentially expressed genes can be found", {
  skip_if_not_installed("Seurat")
  expect_equal(process_input_seurat(seurat, 0)$gene,
               c("FCN1", "FTL", "CLU"))
})

test_that("correctly handles case when two identical clusters specified", {
  skip_if_not_installed("Seurat")
  expect_equal(process_input_seurat(seu, 0, 0)$gene, character())
})
