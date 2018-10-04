#' Perform differential gene expression based on ENCODE data
#'
#' @param counts Data frame: gene expression
#' @param geneAnnot Data frame: gene annotation
#'
#' @importFrom stats model.matrix aggregate
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom AnnotationDbi select keys
#' @importFrom org.Hs.eg.db org.Hs.eg.db org.Hs.egENSEMBL2EG
#'
#' @export
#' @return Data frame with differential gene expression results between
#' knockdown and control
performDifferentialExpression <- function(counts, geneAnnot) {
    counts <- data.frame(counts)
    rownames(counts) <- counts$gene_id

    # Design matrix
    Sample_info <- data.frame(
        sample = c("shRNA1", "shRNA2", "control1", "control2"),
        condition = c("shRNA", "shRNA", "control", "control"))
    design <- model.matrix(~ condition, Sample_info)
    rownames(design) <- Sample_info$sample

    # Check: identical(names(counts[ , 3:6]), rownames(design_matrix))
    voom <- voom(counts[ , 3:6], design=design, plot=FALSE,
                 normalize.method="quantile")

    # Fit linear model
    fit     <- lmFit(voom[ , colnames(voom$E) %in% rownames(design)],
                     design=design)
    ebayes  <- eBayes(fit)
    results <- topTable(ebayes, coef = 2, number = nrow(ebayes),
                        sort.by = "logFC", resort.by = "p")

    # Convert to gene symbol
    geneConversion <- suppressMessages(
        select(org.Hs.eg.db, keys=keys(org.Hs.egENSEMBL2EG),
               columns=c("ENSEMBL", "SYMBOL"), keytype="ENSEMBL"))
    results$Gene_symbol <- geneConversion$SYMBOL[
        match(gsub("\\..*", "", rownames(results)), geneConversion$ENSEMBL)]

    # Mean-aggregation per gene symbol to compare unique gene knockdowns
    results2 <- aggregate(results[ , 1:6], data=results, FUN=mean,
                          by=list(Gene_symbol=results$Gene_symbol))
    return(results2)
}

#' Download data if a file does not exist
#'
#' @param file Character: filepath
#' @param link Character: link to download file
#' @param gz Boolean: is downloaded file compressed?
downloadIfNeeded <- function(file, link, gz=TRUE) {
    if (!file.exists(file)) {
        if (gz) file <- paste0(file, ".gz")
        download.file(link, file)
        if (gz) gunzip(file)
    }
}
