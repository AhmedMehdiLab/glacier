context("Process annotations")
library(dplyr)
library(magrittr)
library(tibble)

anno_1 <- file.path("files", "anno_one.csv")
anno_2 <- file.path("files", "anno_two.csv")

test_that("simple annotations can be processed", {
  anno_pre <- import_annotations(read_common(anno_1, ",", F), c(2, 2), 3)
  info <- anno_pre[c("name", "info")]

  expect_equal(
    process_annotations(anno_pre, info, ""),
    list(
      gs_annos = anno_pre %>% select("name"),
      annos = NULL
    )
  )
  expect_equal(
    process_annotations(anno_pre, info, "name"),
    list(
      gs_annos = anno_pre %>%
        select("name") %>%
        add_column(anno_name = "SET_1"),
      annos = "SET_1"
    )
  )
  expect_equal(
    process_annotations(anno_pre, info, "syms"),
    list(
      gs_annos = anno_pre %>% select("name") %>% add_column(anno_syms = "1"),
      annos = "1"
    )
  )
  expect_equal(
    process_annotations(anno_pre, info, "info"),
    list(
      gs_annos = anno_pre %>%
        select("name") %>%
        add_column(anno_info = "INFO1"),
      annos = "INFO1"
    )
  )
  expect_equal(
    process_annotations(anno_pre, info, "auto"),
    list(
      gs_annos = anno_pre %>%
        select("name") %>%
        add_column(anno_auto = NA_character_),
      annos = character()
    )
  )
  expect_equal(
    process_annotations(anno_pre, info, "file"),
    list(
      gs_annos = anno_pre %>% select("name") %>% add_column(anno_ = "anno_1"),
      annos = "anno_1"
    )
  )
  expect_equal(
    process_annotations(anno_pre, info,
                        c("name", "syms", "info", "auto", "file")),
    list(
      gs_annos = anno_pre %>%
        select("name") %>%
        add_column(
          anno_name = "SET_1",
          anno_syms = "1",
          anno_info = "INFO1",
          anno_auto = NA_character_,
          anno_ = "anno_1"
        ),
      annos = c("SET_1", "1", "INFO1", "anno_1")
    )
  )
})

test_that("complex annotations can be processed", {
  anno_pre <- import_annotations(read_common(anno_2, ",", F), c(2, 3), 4)
  info <- anno_pre[c("name", "info")]

  expect_equal(
    process_annotations(anno_pre, info, ""),
    list(
      gs_annos = anno_pre %>% select("name"),
      annos = NULL
    )
  )
  expect_equal(
    process_annotations(anno_pre, info, "name"),
    list(
      gs_annos = anno_pre %>%
        select("name") %>%
        add_column(anno_name = c("SET_1", "SET_2", "SET_3")),
      annos = c("SET_1", "SET_2", "SET_3")
    )
  )
  expect_equal(
    process_annotations(anno_pre, info, "syms"),
    list(
      gs_annos = anno_pre %>%
        select("name") %>%
        add_column(anno_syms = c("1", "2", "3")),
      annos = c("1", "2", "3")
    )
  )
  expect_equal(
    process_annotations(anno_pre, info, "info"),
    list(
      gs_annos = anno_pre %>%
        select("name") %>%
        add_column(anno_info = c("INFO1", "INFO2", "")),
      annos = c("INFO1", "INFO2")
    )
  )
  expect_equal(
    process_annotations(anno_pre, info, "auto"),
    list(
      gs_annos = anno_pre %>%
        select("name") %>%
        add_column(anno_auto = NA_character_),
      annos = character()
    )
  )
  expect_equal(
    process_annotations(anno_pre, info, "file"),
    list(
      gs_annos = anno_pre %>%
        select("name") %>%
        add_column(anno_1 = c("anno_1", "anno_2", "anno_1"),
                   anno_2 = c("anno_2", NA, NA)),
      annos = c("anno_1", "anno_2")
    )
  )
  expect_equal(
    process_annotations(anno_pre, info,
                        c("name", "syms", "info", "auto", "file")),
    list(
      gs_annos = anno_pre %>%
        select("name") %>%
        add_column(
          anno_name = c("SET_1", "SET_2", "SET_3"),
          anno_syms = c("1", "2", "3"),
          anno_info = c("INFO1", "INFO2", ""),
          anno_auto = NA_character_,
          anno_1 = c("anno_1", "anno_2", "anno_1"),
          anno_2 = c("anno_2", NA, NA)
        ),
      annos = c("SET_1", "SET_2", "SET_3", "1", "2", "3", "INFO1", "INFO2",
                "anno_1", "anno_2")
    )
  )
})
