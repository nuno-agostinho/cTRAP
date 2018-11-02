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
    rownames(counts) <- counts$gene_id

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
    results <- topTable(ebayes, coef=2, number=nrow(ebayes), sort.by="logFC",
                        resort.by="p")

    # Mean-aggregation per gene symbol to compare unique gene knockdowns
    results2 <- aggregate(results[ , 1:6], data=results, FUN=mean,
                          by=list(Gene_symbol=rownames(results)))

    # Remove non-matching genes (if any)
    results2 <- results2[rownames(results2) != "", ]
    return(results2)
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
