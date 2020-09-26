context("Seurat object processing")

seu <- readRDS(file.path("files", "seurat_pico.rds"))

test_that("differentially expressed genes can be found", {
  skip_if_not_installed("Seurat")
  expect_equal(process_input_seurat(seu, 0)$gene,
               c("APOC1", "APOE", "CCL18", "NR1H3", "PLA2G7", "FTL", "MGLL",
                 "GM2A", "LGMN", "CYP27A1", "NRP2", "PLD3", "CTSL"))
})
