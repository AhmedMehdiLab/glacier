% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import.R
\name{import_msigdb}
\alias{import_msigdb}
\title{Import MSigDB XML database file into glacier-specific format}
\usage{
import_msigdb(path)
}
\arguments{
\item{path}{path to file}
}
\value{
glacier-specific imported database
}
\description{
Import MSigDB XML database file into glacier-specific format
}
\examples{
msig_path <- system.file("extdata", "ex_msig.xml", package = "glacier")
data <- import_msigdb(msig_path)
}
