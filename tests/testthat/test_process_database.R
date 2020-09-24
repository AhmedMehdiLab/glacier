context("Process database")
library(magrittr)
library(purrr)
library(tibble)

dpt_1 <- file.path("files", "data_one.csv")
dpt_2 <- file.path("files", "data_two.csv")

test_that("simple databases can be processed", {
  data <- import_database(dpt_1, ",", F, c(2, 2), 3)
  none <- factor("Not assigned")

  expect_equal(
    process_database(data),
    list(
      gs_genes = list() %>% set_names(),
      gs_info = tibble(name = character(), info = character(), category = none,
                       organism = none),
      genes = NULL
    )
  )
  expect_equal(
    process_database(data, "Not assigned", "Not assigned"),
    list(
      gs_genes = list(SET_1 = "GENE1"),
      gs_info = tibble(name = "SET_1", info = "INFO1", category = none,
                       organism = none),
      genes = "GENE1"
    )
  )
})

test_that("complex databases can be processed", {
  data <- import_database(dpt_2, ",", F, c(2, 3), 4)
  none <- factor("Not assigned")

  expect_equal(
    process_database(data),
    list(
      gs_genes = list() %>% set_names(),
      gs_info = tibble(name = character(), info = character(), category = none,
                       organism = none),
      genes = NULL
    )
  )
  expect_equal(
    process_database(data, "Not assigned", "Not assigned"),
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
