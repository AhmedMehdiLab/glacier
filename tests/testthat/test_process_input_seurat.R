context("Seurat object processing")
library(mockery)

seu <- readRDS(file.path("files", "seurat_pico.rds"))

test_that("requires Seurat to function", {
  mockery::stub(process_input_seurat, "requireNamespace", FALSE)
  expect_error(process_input_seurat(seu, 0)$gene,
               "Package 'Seurat' is required for this feature")
})

test_that("differentially expressed genes can be found", {
  skip_if_not_installed("Seurat")
  expect_equal(process_input_seurat(seu, 0)$gene,
               c("APOC1", "APOE", "CCL18", "NR1H3", "PLA2G7", "FTL", "MGLL",
                 "GM2A", "LGMN", "CYP27A1", "NRP2", "PLD3", "CTSL"))
})

test_that("correctly handles case when two identical clusters specified", {
  skip_if_not_installed("Seurat")
  expect_equal(process_input_seurat(seu, 0, 0)$gene, character())
})
