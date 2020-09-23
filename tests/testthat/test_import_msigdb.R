context("Import MSigDB XML files")
library(magrittr)
library(purrr)
library(tibble)

ZERO <- file.path("files", "msig_zero.xml")
ONE <- file.path("files", "msig_one.xml")
TWO <- file.path("files", "msig_two.xml")

test_that("blank MSigDB XML files can be imported", {
  expect_equal(
    import.msigdb_xml(ZERO),
    list(
      gs_genes = list() %>% set_names(character()),
      gs_info = tibble(
        name = character(), info = character(), desc = factor(),
        category = factor(), organism = factor()
      )
    )
  )
})

test_that("simple MSigDB XML files can be imported", {
  expect_equal(
    import.msigdb_xml(ONE),
    list(
      gs_genes = list(SN1 = "GENE1"),
      gs_info = tibble(
        name = "SN1", info = "example 1", desc = factor("example one"),
        category = as.factor("1"), organism = as.factor("Homo sapiens")
      )
    )
  )
})

test_that("complex MSigDB XML files can be imported", {
  expect_equal(
    import.msigdb_xml(TWO),
    list(
      gs_genes = list(SN1 = "GENE1", SN2 = c("GENE1", "GENE2")),
      gs_info = tibble(
        name = c("SN1", "SN2"), info = c("example 1", "example 2"),
        desc = factor(c("example one", "example two")), category = as.factor("1"),
        organism = as.factor("Homo sapiens")
      )
    )
  )
})