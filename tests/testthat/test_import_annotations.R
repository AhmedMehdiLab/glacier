context("Import annotation files")
library(magrittr)
library(tibble)

ONE <- file.path("files", "anno_one.csv")
TWO <- file.path("files", "anno_two.csv")

test_that("blank annotation files cannot be imported", {
  expect_error(import.annotations(tibble()), "anno.raw is empty")
  expect_error(import.annotations(tibble(), c(2, 3), 4), "anno.raw is empty")
})

test_that("simple annotation files can be imported", {
  expect_equal(
    import.annotations(read.common(ONE, ",", F), c(2, 2), 3),
    tibble(name = "SET_1", anno. = "anno_1", info = "INFO1")
  )
  expect_equal(
    import.annotations(read.common(ONE, ",", F), c(2, 2), 0),
    tibble(name = "SET_1", anno. = "anno_1", info = "")
  )
  expect_equal(
    import.annotations(read.common(ONE, ",", F)),
    tibble(name = "SET_1", anno.1 = "anno_1", anno.2 = "INFO1", info = "")
  )
})

test_that("complex annotation files can be imported", {
  expect_equal(
    import.annotations(read.common(TWO, ",", F), c(2, 3), 4),
    tibble(name = c("SET_1", "SET_2", "SET_3"), anno.1 = c("anno_1", "anno_2", "anno_1"), anno.2 = c("anno_2", NA, NA), info = c("INFO1", "INFO2", ""))
  )
  expect_equal(
    import.annotations(read.common(TWO, ",", F), c(2, 3), 0),
    tibble(name = c("SET_1", "SET_2", "SET_3"), anno.1 = c("anno_1", "anno_2", "anno_1"), anno.2 = c("anno_2", NA, NA), info = c("", "", ""))
  )
  expect_equal(
    import.annotations(read.common(TWO, ",", F)),
    tibble(name = c("SET_1", "SET_2", "SET_3"), anno.1 = c("anno_1", "anno_2", "anno_1"), anno.2 = c("anno_2", NA, NA), anno.3 = c("INFO1", "INFO2", NA), info = c("", "", ""))
  )
})