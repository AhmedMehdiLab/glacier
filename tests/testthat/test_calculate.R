context("Calculate statistics")
library(magrittr)
library(tibble)

anno_1 <- file.path("files", "anno_one.csv")
anno_2 <- file.path("files", "anno_two.csv")
data_1 <- file.path("files", "data_one.csv")
data_2 <- file.path("files", "data_two.csv")

test_that("calculations can be performed on simple inputs", {
  anno_pre <- import_annotations(read_common(anno_1, ",", F), c(2, 2), 3)
  data_pre <- import_database(read_common(data_1, ",", F), c(2, 2), 3)
  info <- anno_pre[c("name", "info")]

  anno <- process_annotations(anno_pre, info, "file")
  data <- process_database(data_pre)
  input <- process_input("")$input
  expect_error(
    calculate(input, anno$annos, anno$gs_annos, data$gs_genes, 10000),
    NA
  )

  input <- process_input("a")$input
  expect_error(
    calculate(input, anno$annos, anno$gs_annos, data$gs_genes, 10000),
    NA
  )

  input <- process_input("a 0.1 b 0.3 c 0.8 d e f g 0.9")$input
  expect_error(
    calculate(input, anno$annos, anno$gs_annos, data$gs_genes, 10000),
    NA
  )

  input <- process_input("GENE1 GENE2 GENE3 0.2")$input
  expect_error(
    calculate(input, anno$annos, anno$gs_annos, data$gs_genes, 10000),
    NA
  )
})

test_that("calculations can be performed on complex inputs", {
  anno_pre <- import_annotations(read_common(anno_2, ",", F), c(2, 3), 4)
  data_pre <- import_database(read_common(data_2, ",", F), c(2, 3), 4)
  info <- anno_pre[c("name", "info")]

  anno <- process_annotations(anno_pre, info, "file")
  data <- process_database(data_pre)
  input <- process_input("")$input
  expect_error(
    calculate(input, anno$annos, anno$gs_annos, data$gs_genes, 10000),
    NA
  )

  input <- process_input("a")$input
  expect_error(
    calculate(input, anno$annos, anno$gs_annos, data$gs_genes, 10000),
    NA
  )

  input <- process_input("a 0.1 b 0.3 c 0.8 d e f g 0.9")$input
  expect_error(
    calculate(input, anno$annos, anno$gs_annos, data$gs_genes, 10000),
    NA
  )

  input <- process_input("GENE1 GENE2 GENE3 0.2")$input
  expect_error(
    calculate(input, anno$annos, anno$gs_annos, data$gs_genes, 10000),
    NA
  )
})
