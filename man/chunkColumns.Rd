% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare.R
\name{chunkColumns}
\alias{chunkColumns}
\title{Assign columns into chunks}
\usage{
chunkColumns(x, nrows, chunkGiB)
}
\arguments{
\item{x}{Vector of elements}

\item{nrows}{Numeric: number of rows}

\item{chunkGiB}{Numeric: size (in gibibytes) of chunks to load
\code{reference} file; only if argument \code{reference} is a file path}
}
\value{
List of chunks with equally distributed columns
}
\description{
Assign columns into chunks
}
\keyword{internal}
