% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/display.R
\name{plot_auc}
\alias{plot_auc}
\title{Plot value of AUC}
\usage{
plot_auc(aucs, variable)
}
\arguments{
\item{aucs}{value from \code{\link{score_seurat}} or
\code{\link{score_expr}}}

\item{variable}{value to plot e.g. "pca"}
}
\value{
ggplot2: AUC
}
\description{
Plot value of AUC
}
\examples{
seu_path <- system.file("extdata", "ex_seurat.rds", package = "glacier")
seurat <- readRDS(seu_path)

results <- score_seurat(seurat, "grp", c("APOE", "CTSZ"))
auc <- plot_auc(results$aucs, "pca")
}
