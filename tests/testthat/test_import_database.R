context("Import database files")
library(tibble)

dpt_0 <- file.path("files", "blank.csv")
dpt_1 <- file.path("files", "data_one.csv")
dpt_2 <- file.path("files", "data_two.csv")

test_that("blank databases cannot be imported", {
  expect_error(suppressWarnings(import_database(dpt_0, ",", F)),
               "File is empty")
  expect_error(suppressWarnings(import_database(dpt_0, ",", F, c(1, 2), 3)),
               "File is empty")
})

test_that("simple databases can be imported", {
  none <- factor("Not assigned")
  expect_equal(
    import_database(dpt_1, ",", F, c(2, 2), 3),
    list(
      gs_genes = list(SET_1 = "GENE1"),
      gs_info = tibble(name = "SET_1", info = "INFO1", category = none,
                       organism = none)
    )
  )
  expect_equal(
    import_database(dpt_1, ",", F, c(2, 2), 0),
    list(
      gs_genes = list(SET_1 = "GENE1"),
      gs_info = tibble(name = "SET_1", info = "", category = none,
                       organism = none)
    )
  )
  expect_equal(
    import_database(dpt_1, ",", F),
    list(
      gs_genes = list(SET_1 = c("GENE1", "INFO1")),
      gs_info = tibble(name = "SET_1", info = "", category = none,
                       organism = none)
    )
  )
})

test_that("complex databases can be imported", {
  none <- factor("Not assigned")
  expect_equal(
    import_database(dpt_2, ",", F, c(2, 3), 4),
    list(
      gs_genes = list(SET_1 = c("GENE1", "GENE2"), SET_2 = "GENE3",
                      SET_3 = "GENE4"),
      gs_info = tibble(
        name = c("SET_1", "SET_2", "SET_3"),
        info = c("INFO1", "INFO2", ""),
        category = none,
        organism = none
      )
    )
  )
  expect_equal(
    import_database(dpt_2, ",", F, c(2, 3), 0),
    list(
      gs_genes = list(SET_1 = c("GENE1", "GENE2"), SET_2 = "GENE3",
                      SET_3 = "GENE4"),
      gs_info = tibble(
        name = c("SET_1", "SET_2", "SET_3"),
        info = c("", "", ""),
        category = none,
        organism = none
      )
    )
  )
  expect_equal(
    import_database(dpt_2, ",", F),
    list(
      gs_genes = list(SET_1 = c("GENE1", "GENE2", "INFO1"),
                      SET_2 = c("GENE3", "INFO2"),
                      SET_3 = "GENE4"),
      gs_info = tibble(
        name = c("SET_1", "SET_2", "SET_3"),
        info = c("", "", ""),
        category = none,
        organism = none
      )
    )
  )
})
