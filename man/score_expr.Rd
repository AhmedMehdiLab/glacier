% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score.R
\name{score_expr}
\alias{score_expr}
\title{Calculate expression scores}
\usage{
score_expr(expression, group, genes)
}
\arguments{
\item{expression}{expression data}

\item{group}{grouping variable from expression data}

\item{genes}{genes of interest}
}
\value{
expression pca, pls and plsda scores with AUC
}
\description{
Calculate expression scores
}
\examples{
expr_path <- system.file("extdata", "ex_expr.rds", package = "glacier")
expr <- readRDS(expr_path)

results <- score_expr(expr, "grp", c("APOE", "CTSZ"))
}
