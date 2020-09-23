context("Additional tests")
library(magrittr)
library(tibble)

ANNO <- file.path("files", "anno_two.csv")
DATA <- file.path("files", "data_two.csv")

test_that("function 'display.overlap' does not crash", {
  anno.pre <-
    import.annotations(read.common(ANNO, ",", F), c(2, 3), 4)
  data.pre <- import.database(read.common(DATA, ",", F), c(2, 3), 4)
  input <- process.input('GENE1 0.1 GENE2 0.2 GENE3 0.3')$input
  stats <-
    workflow('GENE1 0.1 GENE2 0.2 GENE3 0.3',
             anno.pre,
             data.pre = data.pre,
             universe = 10000)
  
  expect_error(display.overlap(input, stats$stats, stats$matches, "Gene Value"),
               NA)
  expect_error(display.overlap(input, stats$stats, stats$matches, "Adjusted P-value"),
               NA)
})

test_that("function 'display.stats' does not crash", {
  anno.pre <-
    import.annotations(read.common(ANNO, ",", F), c(2, 3), 4)
  data.pre <- import.database(read.common(DATA, ",", F), c(2, 3), 4)
  stats <-
    workflow('GENE1 0.1 GENE2 0.2 GENE3 0.3',
             anno.pre,
             data.pre = data.pre,
             universe = 10000)
  
  expect_error(display.stats(stats$stats, "Fold Enrichment", "Adjusted P-value"),
               NA)
  expect_error(display.stats(stats$stats, "Fold Enrichment", "Adjusted P-value", sort = T),
               NA)
})

test_that("function 'workflow' does not crash", {
  anno.pre <-
    import.annotations(read.common(ANNO, ",", F), c(2, 3), 4)
  data.pre <- import.database(read.common(DATA, ",", F), c(2, 3), 4)
  
  expect_error(
    workflow(
      'GENE1 0.1 GENE2 0.2 GENE3 0.3',
      anno.pre,
      data.pre = data.pre,
      universe = 10000
    ),
    NA
  )
})