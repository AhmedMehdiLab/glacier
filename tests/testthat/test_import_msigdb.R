context("Import MSigDB XML files")
library(magrittr)
library(tibble)

msig_1 <- system.file("extdata", "ex_msig.xml", package = "glacier")

test_that("MSigDB XML files can be imported", {
  expect_equal(
    import_msigdb(msig_1),
    list(
      gs_genes = list(SET_1 = c("FCN1", "FTL", "CLU"),
                      SET_2 = c("PDK4", "PFKFB3"),
                      SET_3 = c("Bim", "CLU", "FoxO3a"),
                      SET_4 = c("FTL", "SOX2", "FCN1")),
      gs_info = tibble(
        name = c("SET_1", "SET_2", "SET_3", "SET_4"),
        info = character(1),
        desc = factor(""),
        category = factor("Not assigned"),
        organism = factor("Not assigned")
      )
    )
  )
})
