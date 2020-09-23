context("Import MSigDB XML files")
library(magrittr)
library(purrr)
library(tibble)

msig_0 <- file.path("files", "msig_zero.xml")
msig_1 <- file.path("files", "msig_one.xml")
msig_2 <- file.path("files", "msig_two.xml")

test_that("blank MSigDB XML files can be imported", {
  expect_equal(
    import_msigdb_xml(msig_0),
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
    import_msigdb_xml(msig_1),
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
    import_msigdb_xml(msig_2),
    list(
      gs_genes = list(SN1 = "GENE1", SN2 = c("GENE1", "GENE2")),
      gs_info = tibble(
        name = c("SN1", "SN2"),
        info = c("example 1", "example 2"),
        desc = factor(c("example one", "example two")),
        category = as.factor("1"),
        organism = as.factor("Homo sapiens")
      )
    )
  )
})
