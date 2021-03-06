% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsetting flowsets.R
\name{matrix_to_flowframe}
\alias{matrix_to_flowframe}
\title{A wrapper to flowFrame constructor that also updates the "desc" parameters based on fs_src}
\usage{
matrix_to_flowframe(matrix)
}
\arguments{
\item{matrix}{A numeric matrix}
}
\value{
a flowFrame whose expression correspond to the matrix provided. If some channels are common with fs_src, the function updates the flowframe with the corresponding description (ie protein targets)
}
\description{
A wrapper to flowFrame constructor that also updates the "desc" parameters based on fs_src
}
\details{
If replace is true, each element ff will have n events. If it is false, ff will have min(n,nrow(ff))
}
