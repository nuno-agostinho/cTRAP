#' L1000 package
#'
#' Compare differential gene expression results with those from big datasets
#' (e.g. L1000), allowing to infer which types of perturbations may explain the
#' observed difference in gene expression.
#'
#' \strong{Input:} To use this package, a named vector of differentially
#' expressed gene metric is needed, where its values represent the significance
#' and magnitude of the differentially expressed genes (e.g. t-statistic) and
#' its names are gene symbols.
#'
#' \strong{Workflow:} The differentially expressed genes will be compared
#' against selected perturbation conditions by:
#' \itemize{
#'     \item{Spearman or Pearson correlation with z-scores of differentially
#'     expressed genes after perturbations from L1000. Use function
#'     \code{compareAgainstL1000} with \code{method = "spearman"} or
#'     \code{method = "pearson"}}
#'     \item{Gene set enrichment analysis (GSEA) using the (around) 12 000 genes
#'     from L1000. Use function \code{compareAgainstL1000} with
#'     \code{method = gsea}.}
#' }
#'
#' Available perturbation conditions for L1000 include:
#' \itemize{
#'     \item{Cell line(s).}
#'     \item{Perturbation type (gene knockdown, gene upregulation or drug
#'     intake).}
#'     \item{Drug concentration.}
#'     \item{Time points.}
#' }
#'
#' Values for each perturbation type can be listed with
#' \code{getL1000PerturbationTypes()}
#'
#' \strong{Output:} The output includes a data frame of ranked perturbations
#' based on the associated statistical values and respective p-values.
#'
#' @name l1000
#' @docType package
NULL
