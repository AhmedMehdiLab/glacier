context("Additional tests")
library(magrittr)
library(tibble)

anno_file <- file.path("files", "anno_two.csv")
data_file <- file.path("files", "data_two.csv")

test_that("function 'explore_annotation' supports gene filtering", {
  anno <- import_annotations(anno_file, ",", F, c(2, 3), 4)
  data <- import_database(data_file, ",", F, c(2, 3), 4)
  input <- process_input_text("GENE1 0.1 GENE2 0.2 GENE3 0.3")

  info <- anno[c("name", "info")]
  anno_proc <- glacier:::process_annotations(anno, info, "file")
  data_proc <- glacier:::process_database(data)

  expect_equal(
    explore_annotation("anno_1", anno_proc$gs_annos,
                       data_proc$gs_genes)$genes,
    c("GENE1", "GENE2", "GENE4")
  )
  expect_equal(
    explore_annotation("anno_1", anno_proc$gs_annos, data_proc$gs_genes,
                       "GENE1")$genes,
    c("GENE1")
  )
})

test_that("function 'compute' does not crash", {
  anno <- import_annotations(anno_file, ",", F, c(2, 3), 4)
  data <- import_database(data_file, ",", F, c(2, 3), 4)
  input <- process_input_text("GENE1 0.1 GENE2 0.2 GENE3 0.3")

  expect_error(stats <- compute(input, anno, data), NA)
  expect_error(stats <- compute(input, anno, data, 10000), NA)
  expect_error(stats <- compute(input, anno, data, save = tempfile()), NA)
})

test_that("function 'plot_overlap' does not crash", {
  anno <- import_annotations(anno_file, ",", F, c(2, 3), 4)
  data <- import_database(data_file, ",", F, c(2, 3), 4)
  input <- process_input_text("GENE1 0.1 GENE2 0.2 GENE3 0.3")
  res <- compute(input, anno, data, 10000)

  expect_error(plot_overlap(res$matches, "Gene Value", input, res$stats), NA)
  expect_error(plot_overlap(res$matches, "P-value", input, res$stats), NA)
})

test_that("function 'plot_stats' does not crash", {
  anno <- import_annotations(anno_file, ",", F, c(2, 3), 4)
  data <- import_database(data_file, ",", F, c(2, 3), 4)
  input <- process_input_text("GENE1 0.1 GENE2 0.2 GENE3 0.3")
  stats <- compute(input, anno, data, 10000)

  expect_error(plot_stats(stats$stats, "# genes", "P-value"), NA)
  expect_error(plot_stats(stats$stats, "# genes", "P-value", sort_y = T), NA)
})
