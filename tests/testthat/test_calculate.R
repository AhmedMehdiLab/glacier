context("Calculate statistics")
library(magrittr)
library(tibble)

apt_1 <- file.path("files", "anno_one.csv")
apt_2 <- file.path("files", "anno_two.csv")
dpt_1 <- file.path("files", "data_one.csv")
dpt_2 <- file.path("files", "data_two.csv")

test_that("calculations can be performed on simple inputs", {
  anno <- import_annotations(apt_1, ",", F, c(2, 2), 3)
  data <- import_database(dpt_1, ",", F, c(2, 2), 3)
  info <- anno[c("name", "info")]

  anno_proc <- process_annotations(anno, info, "file")
  data_proc <- process_database(data)
  input <- process_input_text("")
  expect_error(calculate(input, anno_proc$annos, anno_proc$gs_annos,
                         data_proc$gs_genes, 10000), NA)

  input <- process_input_text("a")
  expect_error(calculate(input, anno_proc$annos, anno_proc$gs_annos,
                         data_proc$gs_genes, 10000), NA)

  input <- process_input_text("a 0.1 b 0.3 c 0.8 d e f g 0.9")
  expect_error(calculate(input, anno_proc$annos, anno_proc$gs_annos,
                         data_proc$gs_genes, 10000), NA)

  input <- process_input_text("GENE1 GENE2 GENE3 0.2")
  expect_error(calculate(input, anno_proc$annos, anno_proc$gs_annos,
                         data_proc$gs_genes, 10000), NA)
})

test_that("calculations can be performed on complex inputs", {
  anno <- import_annotations(apt_2, ",", F, c(2, 3), 4)
  data <- import_database(dpt_2, ",", F, c(2, 3), 4)
  info <- anno[c("name", "info")]

  anno_proc <- process_annotations(anno, info, "file")
  data_proc <- process_database(data)
  input <- process_input_text("")
  expect_error(calculate(input, anno_proc$annos, anno_proc$gs_annos,
                         data_proc$gs_genes, 10000), NA)

  input <- process_input_text("a")
  expect_error(calculate(input, anno_proc$annos, anno_proc$gs_annos,
                         data_proc$gs_genes, 10000), NA)

  input <- process_input_text("a 0.1 b 0.3 c 0.8 d e f g 0.9")
  expect_error(calculate(input, anno_proc$annos, anno_proc$gs_annos,
                         data_proc$gs_genes, 10000), NA)

  input <- process_input_text("GENE1 GENE2 GENE3 0.2")
  expect_error(calculate(input, anno_proc$annos, anno_proc$gs_annos,
                         data_proc$gs_genes, 10000), NA)
})
