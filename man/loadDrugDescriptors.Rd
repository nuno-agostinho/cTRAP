% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drugSetEnrichment.R
\name{loadDrugDescriptors}
\alias{loadDrugDescriptors}
\title{Load table with drug descriptors}
\usage{
loadDrugDescriptors(
  source = c("NCI60", "CMap"),
  type = c("2D", "3D"),
  file = NULL
)
}
\arguments{
\item{source}{Character: molecular descriptors for compounds in \code{NCI60}
or \code{CMap}}

\item{type}{Character: load \code{2D} or \code{3D} molecular descriptors}

\item{file}{Character: filepath to drug descriptors (automatically downloaded
if file does not exist)}
}
\value{
Data table with drug descriptors
}
\description{
Load table with drug descriptors
}
\examples{
loadDrugDescriptors()
}
\seealso{
Other functions for drug set enrichment analysis: 
\code{\link{analyseDrugSetEnrichment}()},
\code{\link{plotDrugSetEnrichment}()},
\code{\link{prepareDrugSets}()}
}
\concept{functions for drug set enrichment analysis}
