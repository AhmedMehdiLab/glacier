context("Process database")
library(dplyr)
library(magrittr)
library(purrr)
library(tibble)

data_1 <- file.path("files", "data_one.csv")
data_2 <- file.path("files", "data_two.csv")

test_that("simple databases can be processed", {
  data_pre <- import_database(read_common(data_1, ",", F), c(2, 2), 3)
  none <- factor("Not assigned")

  expect_equal(
    process_database(data_pre),
    list(
      gs_genes = list() %>% set_names(),
      gs_info = tibble(name = character(), info = character(), category = none,
                       organism = none),
      genes = NULL
    )
  )
  expect_equal(
    process_database(data_pre, "Not assigned", "Not assigned"),
    list(
      gs_genes = list(SET_1 = "GENE1"),
      gs_info = tibble(name = "SET_1", info = "INFO1", category = none,
                       organism = none),
      genes = "GENE1"
    )
  )
})

test_that("complex databases can be processed", {
  data_pre <- import_database(read_common(data_2, ",", F), c(2, 3), 4)
  none <- factor("Not assigned")

  expect_equal(
    process_database(data_pre),
    list(
      gs_genes = list() %>% set_names(),
      gs_info = tibble(name = character(), info = character(), category = none,
                       organism = none),
      genes = NULL
    )
  )
  expect_equal(
    process_database(data_pre, "Not assigned", "Not assigned"),
    list(
      gs_genes = list(SET_1 = c("GENE1", "GENE2"), SET_2 = "GENE3",
                      SET_3 = "GENE4"),
      gs_info = tibble(
        name = c("SET_1", "SET_2", "SET_3"),
        info = c("INFO1", "INFO2", ""),
        category = none,
        organism = none
      ),
      genes = c("GENE1", "GENE2", "GENE3", "GENE4")
    )
  )
})
