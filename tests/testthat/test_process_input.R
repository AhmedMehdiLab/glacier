context("Gene input pre-processing")
library(magrittr)
library(tibble)

BASE <- tibble(gene = character(), value = numeric())

test_that("works for well-formed single gene input", {
  expect_equal(process.input.pre(""), BASE)
  expect_equal(process.input.pre("a"), BASE %>% add_row(gene = "a"))
  expect_equal(process.input.pre("0.1"), BASE)
  expect_equal(process.input.pre("a 0.1"), BASE %>% add_row(gene = "a", value = 0.1))
  
  expect_equal(process.input.pre(" "), BASE)
  expect_equal(process.input.pre("a "), BASE %>% add_row(gene = "a"))
  expect_equal(process.input.pre("0.1 "), BASE)
  expect_equal(process.input.pre("a 0.1 "), BASE %>% add_row(gene = "a", value = 0.1))
})

test_that("works for well-formed multiple gene input", {
  expect_equal(process.input.pre("a b"), BASE %>% add_row(gene = "a") %>% add_row(gene = "b"))
  expect_equal(process.input.pre("a b c"), BASE %>% add_row(gene = "a") %>% add_row(gene = "b") %>% 
    add_row(gene = "c"))
  expect_equal(process.input.pre("a 0.1 b 0.2"), BASE %>% add_row(gene = "a", value = 0.1) %>% 
    add_row(gene = "b", value = 0.2))
  expect_equal(process.input.pre("a 0.1 b 0.2 c 0.1"), BASE %>% add_row(gene = "a", 
    value = 0.1) %>% add_row(gene = "b", value = 0.2) %>% add_row(gene = "c", 
    value = 0.1))
  
  expect_equal(process.input.pre("a b "), BASE %>% add_row(gene = "a") %>% add_row(gene = "b"))
  expect_equal(process.input.pre("a b c "), BASE %>% add_row(gene = "a") %>% add_row(gene = "b") %>% 
    add_row(gene = "c"))
  expect_equal(process.input.pre("a 0.1 b 0.1"), BASE %>% add_row(gene = "a", value = 0.1) %>% 
    add_row(gene = "b", value = 0.1))
  expect_equal(process.input.pre("a 0.1 b 0.1 c 0.1 "), BASE %>% add_row(gene = "a", 
    value = 0.1) %>% add_row(gene = "b", value = 0.1) %>% add_row(gene = "c", 
    value = 0.1))
})

test_that("eliminates duplicate genes", {
  expect_equal(process.input.pre(" "), BASE)
  expect_equal(process.input.pre("  "), BASE)
  expect_equal(process.input.pre("\n\r\t \n\r\t "), BASE)
  
  expect_equal(process.input.pre("a a"), BASE %>% add_row(gene = "a"))
  expect_equal(process.input.pre("a b a "), BASE %>% add_row(gene = "a") %>% add_row(gene = "b"))
  expect_equal(process.input.pre("a a b b "), BASE %>% add_row(gene = "a") %>% 
    add_row(gene = "b"))
  expect_equal(process.input.pre("a b a b a b "), BASE %>% add_row(gene = "a") %>% 
    add_row(gene = "b"))
  expect_equal(process.input.pre("a b c b a a a b b c c c a b c "), BASE %>% add_row(gene = "a") %>% 
    add_row(gene = "b") %>% add_row(gene = "c"))
})

test_that("eliminates duplicate values", {
  expect_equal(process.input.pre("a 0.1 0.2 0.3"), BASE %>% add_row(gene = "a", 
    value = 0.1))
  expect_equal(process.input.pre("a 0.1 b 0.2 a 0.2 b 0.4 a 0.2 0.4"), BASE %>% 
    add_row(gene = "a", value = 0.1) %>% add_row(gene = "b", value = 0.2))
  expect_equal(process.input.pre("a 0.1 b c 0.1 b 0.1 c 0.4 a 0.8"), BASE %>% add_row(gene = "a", 
    value = 0.1) %>% add_row(gene = "b") %>% add_row(gene = "c", value = 0.1))
})

test_that("post-processing works", {
  expect_equal(process.input.post(process.input.pre("")), list(input = BASE, missing = logical()))
  expect_equal(process.input.post(process.input.pre("a 0.1 b")), list(input = BASE %>% 
    add_row(gene = "a", value = 0.1) %>% add_row(gene = "b", value = 0), missing = c(F, 
    T)))
  expect_equal(process.input.post(process.input.pre("a 0.1 b c d e 0.4")), list(input = BASE %>% 
    add_row(gene = "a", value = 0.1) %>% add_row(gene = "b", value = 0) %>% add_row(gene = "c", 
    value = 0) %>% add_row(gene = "d", value = 0) %>% add_row(gene = "e", value = 0.4), 
    missing = c(F, T, T, T, F))
  )
})

test_that("combined processing works", {
  expect_equal(process.input(""), list(input = BASE, missing = logical()))
  expect_equal(process.input("a 0.1 b"), list(
    input = BASE %>% add_row(gene = "a", value = 0.1) %>% add_row(gene = "b", value = 0), missing = c(F, T)))
  expect_equal(
    process.input("a 0.1 b c d e 0.4"),
    list(input = BASE %>% 
           add_row(gene = "a", value = 0.1) %>%
           add_row(gene = "b", value = 0) %>%
           add_row(gene = "c",  value = 0) %>%
           add_row(gene = "d", value = 0) %>%
           add_row(gene = "e", value = 0.4), 
         missing = c(F, T, T, T, F))
  )
})
