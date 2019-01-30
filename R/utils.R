#' Perform differential gene expression based on ENCODE data
#'
#' @param counts Data frame: gene expression
#'
#' @importFrom stats model.matrix aggregate
#' @importFrom limma voom lmFit eBayes topTable
#'
#' @return Data frame with differential gene expression results between
#' knockdown and control
#' @export
#'
#' @examples
#' data("ENCODEsamples")
#'
#' ## Download ENCODE metadata for a specific cell line and gene
#' # cellLine <- "HepG2"
#' # gene <- "EIF4G1"
#' # ENCODEmetadata <- downloadENCODEknockdownMetadata(cellLine, gene)
#'
#' ## Download samples based on filtered ENCODE metadata
#' # ENCODEsamples <- downloadENCODEsamples(ENCODEmetadata)
#'
#' counts <- prepareENCODEgeneExpression(ENCODEsamples)
#'
#' # Remove low coverage (at least 10 counts shared across two samples)
#' minReads   <- 10
#' minSamples <- 2
#' filter <- rowSums(counts[ , -c(1, 2)] >= minReads) >= minSamples
#' counts <- counts[filter, ]
#'
#' ## Convert ENSEMBL identifier to gene symbol
#' # library(biomaRt)
#' # mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#' # genes <- sapply(strsplit(counts$gene_id, "\\."), `[`, 1)
#' # geneConversion <- getBM(filters="ensembl_gene_id", values=genes, mart=mart,
#' #                         attributes=c("ensembl_gene_id", "hgnc_symbol"))
#' # counts$gene_id <- geneConversion$hgnc_symbol[
#' #     match(genes, geneConversion$ensembl_gene_id)]
#'
#' ## Perform differential gene expression analysis
#' # diffExpr <- performDifferentialExpression(counts)
performDifferentialExpression <- function(counts) {
    counts <- data.frame(counts)

    # Design matrix
    Sample_info <- data.frame(
        sample = c("shRNA1", "shRNA2", "control1", "control2"),
        condition = c("shRNA", "shRNA", "control", "control"))
    design <- model.matrix(~ condition, Sample_info)
    rownames(design) <- Sample_info$sample

    # Check: identical(names(counts[ , 3:6]), rownames(design_matrix))
    voom <- voom(counts[ , -c(1, 2)], design=design, plot=FALSE,
                 normalize.method="quantile")

    # Fit linear model
    fit <- lmFit(voom[ , colnames(voom$E) %in% rownames(design)], design=design)
    ebayes <- eBayes(fit)
    results <- topTable(ebayes, coef=2, number=nrow(ebayes),
                        genelist=counts$gene_id)

    # Mean-aggregation per gene symbol to compare unique gene knockdowns
    meanAggr <- aggregate(results[ , -1], data=results, FUN=mean,
                          by=list(Gene_symbol=results$ID))

    # Remove non-matching genes (if any)
    meanAggr           <- meanAggr[meanAggr$Gene_symbol != "", ]
    rownames(meanAggr) <- meanAggr$Gene_symbol
    return(meanAggr)
}

#' Download data if a file does not exist
#'
#' @param file Character: filepath
#' @param link Character: link to download file
#' @param gz Boolean: is downloaded file compressed?
#'
#' @importFrom utils download.file
#'
#' @return Download file if a file does not exist
#' @keywords internal
downloadIfNeeded <- function(file, link, gz=TRUE) {
    if (!file.exists(file)) {
        if (gz) {
            file <- paste0(file, ".gz")
            mode <- "wb"
        } else {
            mode <- "w"
        }
        download.file(link, file, mode="wb")
        if (gz) gunzip(file)
    }
}

#' Collapse duplicated rows based on a column of a data frame
#'
#' @param df Data frame
#' @param column Character: name of the column with elements to indicate which
#' rows to collapse
#'
#' @importFrom plyr ldply
#'
#' @return Data frame with duplicated rows collapsed
#' @keywords internal
collapseDuplicatedRows <- function(df, column) {
    cluster <- split(seq(df[[column]]), df[[column]])
    collapseRowInfo <- function(thisCluster, df) {
        tmp <- df[thisCluster, ]
        apply(tmp, 2, function(x) paste(unique(x), collapse=", "))
    }
    collapsed <- pblapply(cluster, collapseRowInfo, df)
    collapsed <- ldply(collapsed, .id=column)
    return(collapsed)
}

# l1000perturbations object ----------------------------------------------------

#' Subset an \code{l1000perturbations} object
#'
#' @param x \code{l1000perturbations} object
#' @param ... Extra parameters passed to \code{`[`}
#'
#' @return \code{l1000perturbations} object with subset data
#' @export
`[.l1000perturbations` <- function(x, ...) {
    out <- unclass(x)
    out <- `[`(out, ..., drop=FALSE)

    # Inherit the same attributes
    attrs <- attributes(x)
    attrs$dim <- NULL
    attrs$dimnames <- NULL
    attrs$names <- NULL

    # Trim metadata to contain subset information
    if (ncol(x) != ncol(out) && !is.null(attrs$metadata)) {
        samples <- attrs$metadata$sig_id %in% colnames(out)
        attrs$metadata <- attrs$metadata[samples, ]
    }
    attributes(out) <- c(attributes(out), attrs)
    return(out)
}

#' @inherit base::as.data.frame
#' @export
as.data.frame.l1000perturbations <- function(x, ...) {
    as.data.frame(unclass(x), ...)
}

#' @inherit utils::head
#' @export
head.l1000perturbations <- function(x, ...) head(unclass(x), ...)

#' @inherit utils::tail
#' @export
tail.l1000perturbations <- function(x, ...) tail(unclass(x), ...)
