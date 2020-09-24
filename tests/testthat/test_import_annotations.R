context("Import annotation files")
library(tibble)

apt_0 <- file.path("files", "blank.csv")
apt_1 <- file.path("files", "anno_one.csv")
apt_2 <- file.path("files", "anno_two.csv")

test_that("blank annotation files cannot be imported", {
  expect_error(suppressWarnings(import_annotations(apt_0, ",", F)),
               "File is empty")
  expect_error(suppressWarnings(import_annotations(apt_0, ",", F, c(1, 2), 3)),
               "File is empty")
})

test_that("simple annotation files can be imported", {
  expect_equal(
    import_annotations(apt_1, ",", F, c(2, 2), 3),
    tibble(name = "SET_1", anno_ = "anno_1", info = "INFO1")
  )
  expect_equal(
    import_annotations(apt_1, ",", F, c(2, 2), 0),
    tibble(name = "SET_1", anno_ = "anno_1", info = "")
  )
  expect_equal(
    import_annotations(apt_1, ",", F),
    tibble(name = "SET_1", anno_1 = "anno_1", anno_2 = "INFO1", info = "")
  )
})

test_that("complex annotation files can be imported", {
  expect_equal(
    import_annotations(apt_2, ",", F, c(2, 3), 4),
    tibble(
      name = c("SET_1", "SET_2", "SET_3"),
      anno_1 = c("anno_1", "anno_2", "anno_1"),
      anno_2 = c("anno_2", NA, NA),
      info = c("INFO1", "INFO2", "")
    )
  )
  expect_equal(
    import_annotations(apt_2, ",", F, c(2, 3), 0),
    tibble(
      name = c("SET_1", "SET_2", "SET_3"),
      anno_1 = c("anno_1", "anno_2", "anno_1"),
      anno_2 = c("anno_2", NA, NA),
      info = c("", "", "")
    )
  )
  expect_equal(
    import_annotations(apt_2, ",", F),
    tibble(
      name = c("SET_1", "SET_2", "SET_3"),
      anno_1 = c("anno_1", "anno_2", "anno_1"),
      anno_2 = c("anno_2", NA, NA),
      anno_3 = c("INFO1", "INFO2", NA),
      info = c("", "", "")
    )
  )
})
