context("Additional tests")
library(magrittr)
library(tibble)

anno_file <- file.path("files", "anno_two.csv")
data_file <- file.path("files", "data_two.csv")

test_that("function 'display_overlap' does not crash", {
  anno_pre <-
    import_annotations(read_common(anno_file, ",", F), c(2, 3), 4)
  data_pre <- import_database(read_common(data_file, ",", F), c(2, 3), 4)
  input <- process_input("GENE1 0.1 GENE2 0.2 GENE3 0.3")$input
  stats <-
    workflow("GENE1 0.1 GENE2 0.2 GENE3 0.3",
             anno_pre,
             data_pre = data_pre,
             universe = 10000)

  expect_error(display_overlap(input, stats$stats, stats$matches, "Gene Value"),
               NA)
  expect_error(display_overlap(input, stats$stats, stats$matches,
                               "Adjusted P-value"), NA)
})

test_that("function 'display_stats' does not crash", {
  anno_pre <-
    import_annotations(read_common(anno_file, ",", F), c(2, 3), 4)
  data_pre <- import_database(read_common(data_file, ",", F), c(2, 3), 4)
  stats <-
    workflow("GENE1 0.1 GENE2 0.2 GENE3 0.3",
             anno_pre,
             data_pre = data_pre,
             universe = 10000)

  expect_error(display_stats(stats$stats, "Fold Enrichment",
                             "Adjusted P-value"), NA)
  expect_error(display_stats(stats$stats, "Fold Enrichment", "Adjusted P-value",
                             sort = T), NA)
})

test_that("function 'workflow' does not crash", {
  anno_pre <-
    import_annotations(read_common(anno_file, ",", F), c(2, 3), 4)
  data_pre <- import_database(read_common(data_file, ",", F), c(2, 3), 4)

  expect_error(
    workflow(
      "GENE1 0.1 GENE2 0.2 GENE3 0.3",
      anno_pre,
      data_pre = data_pre,
      universe = 10000
    ),
    NA
  )
})
