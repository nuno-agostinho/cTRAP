% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare.R
\name{as.table.referenceComparison}
\alias{as.table.referenceComparison}
\title{Cross Tabulation and Table Creation}
\usage{
\method{as.table}{referenceComparison}(x, ..., clean = TRUE)
}
\arguments{
\item{x}{\code{referenceComparison} object}

\item{...}{Extra parameters not currently used}

\item{clean}{Boolean: only show certain columns (to avoid redundancy)?}
}
\value{
Complete table with metadata based on a \code{targetingDrugs} object
}
\description{
Cross Tabulation and Table Creation
}
\seealso{
Other functions related with the ranking of CMap perturbations: 
\code{\link{filterCMapMetadata}()},
\code{\link{getCMapConditions}()},
\code{\link{getCMapPerturbationTypes}()},
\code{\link{loadCMapData}()},
\code{\link{loadCMapZscores}()},
\code{\link{parseCMapID}()},
\code{\link{plot.perturbationChanges}()},
\code{\link{plot.referenceComparison}()},
\code{\link{plotTargetingDrugsVSsimilarPerturbations}()},
\code{\link{prepareCMapPerturbations}()},
\code{\link{print.similarPerturbations}()},
\code{\link{rankSimilarPerturbations}()}

Other functions related with the prediction of targeting drugs: 
\code{\link{listExpressionDrugSensitivityAssociation}()},
\code{\link{loadExpressionDrugSensitivityAssociation}()},
\code{\link{plot.referenceComparison}()},
\code{\link{plotTargetingDrugsVSsimilarPerturbations}()},
\code{\link{predictTargetingDrugs}()}
}
\concept{functions related with the prediction of targeting drugs}
\concept{functions related with the ranking of CMap perturbations}
