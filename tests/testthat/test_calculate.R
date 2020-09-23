context("Calculate statistics")
library(magrittr)
library(tibble)

ANNO_ONE <- file.path("files", "anno_one.csv")
ANNO_TWO <- file.path("files", "anno_two.csv")
DATA_ONE <- file.path("files", "data_one.csv")
DATA_TWO <- file.path("files", "data_two.csv")

test_that("calculations can be performed on simple inputs", {
  anno.pre <- import.annotations(read.common(ANNO_ONE, ",", F), c(2, 2), 3)
  data.pre <- import.database(read.common(DATA_ONE, ",", F), c(2, 2), 3)
  info <- anno.pre[c('name', 'info')]
  
  anno <- process.annotations(anno.pre, info, "file")
  data <- process.database(data.pre)
  input <- process.input("")$input
  expect_error(
    calculate(input, anno$annos, anno$gs_annos, data$gs_genes, 10000),
    NA
  )
  
  input <- process.input("a")$input
  expect_error(
    calculate(input, anno$annos, anno$gs_annos, data$gs_genes, 10000),
    NA
  )
  
  input <- process.input("a 0.1 b 0.3 c 0.8 d e f g 0.9")$input
  expect_error(
    calculate(input, anno$annos, anno$gs_annos, data$gs_genes, 10000),
    NA
  )
  
  input <- process.input("GENE1 GENE2 GENE3 0.2")$input
  expect_error(
    calculate(input, anno$annos, anno$gs_annos, data$gs_genes, 10000),
    NA
  )
})

test_that("calculations can be performed on complex inputs", {
  anno.pre <- import.annotations(read.common(ANNO_TWO, ",", F), c(2, 3), 4)
  data.pre <- import.database(read.common(DATA_TWO, ",", F), c(2, 3), 4)
  info <- anno.pre[c('name', 'info')]
  
  anno <- process.annotations(anno.pre, info, "file")
  data <- process.database(data.pre)
  input <- process.input("")$input
  expect_error(
    calculate(input, anno$annos, anno$gs_annos, data$gs_genes, 10000),
    NA
  )
  
  input <- process.input("a")$input
  expect_error(
    calculate(input, anno$annos, anno$gs_annos, data$gs_genes, 10000),
    NA
  )
  
  input <- process.input("a 0.1 b 0.3 c 0.8 d e f g 0.9")$input
  expect_error(
    calculate(input, anno$annos, anno$gs_annos, data$gs_genes, 10000),
    NA
  )
  
  input <- process.input("GENE1 GENE2 GENE3 0.2")$input
  expect_error(
    calculate(input, anno$annos, anno$gs_annos, data$gs_genes, 10000),
    NA
  )
})