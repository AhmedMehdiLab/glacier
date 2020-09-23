context("Process annotations")
library(dplyr)
library(magrittr)
library(tibble)

ONE <- file.path("files", "anno_one.csv")
TWO <- file.path("files", "anno_two.csv")

test_that("simple annotations can be processed", {
  anno.pre <- import.annotations(read.common(ONE, ",", F), c(2, 2), 3)
  info <- anno.pre[c('name', 'info')]
  
  expect_equal(
    process.annotations(anno.pre, info, ""),
    list(
      gs_annos = anno.pre %>% select("name"),
      annos = NULL
    )
  )
  expect_equal(
    process.annotations(anno.pre, info, "name"),
    list(
      gs_annos = anno.pre %>% select("name") %>% add_column(anno.name = "SET_1"),
      annos = "SET_1"
    )
  )
  expect_equal(
    process.annotations(anno.pre, info, "syms"),
    list(
      gs_annos = anno.pre %>% select("name") %>% add_column(anno.syms = "1"),
      annos = "1"
    )
  )
  expect_equal(
    process.annotations(anno.pre, info, "info"),
    list(
      gs_annos = anno.pre %>% select("name") %>% add_column(anno.info = "INFO1"),
      annos = "INFO1"
    )
  )
  expect_equal(
    process.annotations(anno.pre, info, "auto"),
    list(
      gs_annos = anno.pre %>% select("name") %>% add_column(anno.auto = NA_character_),
      annos = character()
    )
  )
  expect_equal(
    process.annotations(anno.pre, info, "file"),
    list(
      gs_annos = anno.pre %>% select("name") %>% add_column(anno. = "anno_1"),
      annos = "anno_1"
    )
  )
  expect_equal(
    process.annotations(anno.pre, info, c("name", "syms", "info", "auto", "file")),
    list(
      gs_annos = anno.pre %>% select("name") %>% add_column(anno.name = "SET_1", anno.syms = "1", anno.info = "INFO1", anno.auto = NA_character_, anno. = "anno_1"),
      annos = c("SET_1", "1", "INFO1", "anno_1")
    )
  )
})

test_that("complex annotations can be processed", {
  anno.pre <- import.annotations(read.common(TWO, ",", F), c(2, 3), 4)
  info <- anno.pre[c('name', 'info')]
  
  expect_equal(
    process.annotations(anno.pre, info, ""),
    list(
      gs_annos = anno.pre %>% select("name"),
      annos = NULL
    )
  )
  expect_equal(
    process.annotations(anno.pre, info, "name"),
    list(
      gs_annos = anno.pre %>% select("name") %>% add_column(anno.name = c("SET_1", "SET_2", "SET_3")),
      annos = c("SET_1", "SET_2", "SET_3")
    )
  )
  expect_equal(
    process.annotations(anno.pre, info, "syms"),
    list(
      gs_annos = anno.pre %>% select("name") %>% add_column(anno.syms = c("1", "2", "3")),
      annos = c("1", "2", "3")
    )
  )
  expect_equal(
    process.annotations(anno.pre, info, "info"),
    list(
      gs_annos = anno.pre %>% select("name") %>% add_column(anno.info = c("INFO1", "INFO2", "")),
      annos = c("INFO1", "INFO2")
    )
  )
  expect_equal(
    process.annotations(anno.pre, info, "auto"),
    list(
      gs_annos = anno.pre %>% select("name") %>% add_column(anno.auto = NA_character_),
      annos = character()
    )
  )
  expect_equal(
    process.annotations(anno.pre, info, "file"),
    list(
      gs_annos = anno.pre %>% select("name") %>% add_column(anno.1 = c("anno_1", "anno_2", "anno_1"), anno.2 = c("anno_2", NA, NA)),
      annos = c("anno_1", "anno_2")
    )
  )
  expect_equal(
    process.annotations(anno.pre, info, c("name", "syms", "info", "auto", "file")),
    list(
      gs_annos = anno.pre %>% select("name") %>% add_column(
        anno.name = c("SET_1", "SET_2", "SET_3"), anno.syms = c("1", "2", "3"),
        anno.info = c("INFO1", "INFO2", ""), anno.auto = NA_character_,
        anno.1 = c("anno_1", "anno_2", "anno_1"), anno.2 = c("anno_2", NA, NA)),
      annos = c("SET_1", "SET_2", "SET_3", "1", "2", "3", "INFO1", "INFO2", "anno_1", "anno_2")
    )
  )
})