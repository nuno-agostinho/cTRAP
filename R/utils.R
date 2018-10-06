#' Perform differential gene expression based on ENCODE data
#'
#' @param counts Data frame: gene expression
#'
#' @importFrom stats model.matrix aggregate
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom biomaRt useDataset getBM useMart
#'
#' @return Data frame with differential gene expression results between
#' knockdown and control
#' @export
#'
#' @examples
#' gene <- "EIF4G1"
#' cellLine <- "HepG2"
#'
#' ENCODEmetadata <- downloadENCODEknockdownMetadata(cellLine, gene)
#' counts <- loadENCODEgeneExpression(ENCODEmetadata)
#'
#' # Remove low coverage (at least 10 counts shared across two samples)
#' minReads   <- 10
#' minSamples <- 2
#' filter <- rowSums(counts[ , 3:6] >= minReads) >= minSamples
#' counts <- counts[filter, ]
#'
#' # Perform differential gene expression analysis
#' diffExpr <- performDifferentialExpression(counts)
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
    voom <- voom(counts[ , 3:6], design=design, plot=FALSE,
                 normalize.method="quantile")

    # Fit linear model
    fit     <- lmFit(voom[ , colnames(voom$E) %in% rownames(design)],
                     design=design)
    ebayes  <- eBayes(fit)
    results <- topTable(ebayes, coef = 2, number = nrow(ebayes),
                        sort.by = "logFC", resort.by = "p")

    # Convert to gene symbol
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    genes <- sapply(strsplit(rownames(results), "\\."), `[`, 1)
    geneConversion <- getBM(filters="ensembl_gene_id", values=genes, mart=mart,
                            attributes=c("ensembl_gene_id", "hgnc_symbol"))
    results$Gene_symbol <- geneConversion$hgnc_symbol[
        match(genes, geneConversion$ensembl_gene_id)]

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
#'
#' @importFrom utils download.file
#'
#' @return Download file if a file does not exist
downloadIfNeeded <- function(file, link, gz=TRUE) {
    if (!file.exists(file)) {
        if (gz) file <- paste0(file, ".gz")
        download.file(link, file)
        if (gz) gunzip(file)
    }
}

#' Load internal data for use in vignettes and examples
#' 
#' @param x Character: name of object to load
#' 
#' @importFrom utils getFromNamespace
#' 
#' @return Object
#' @export
#' 
#' @examples 
#' loadInternalData("l1000perturbationsKnockdown")
loadInternalData <- function(x) {
    res <- getFromNamespace(x, pos="package:cTRAP")
    
    if (x == "compareKnockdown")
        l1000perturbationsVar <- "l1000perturbationsKnockdown"
    else if (x == "compareSmallMolecule")
        l1000perturbationsVar <- "l1000perturbationsSmallMolecules"
    
    if (x %in% paste0("compare", c("Knockdown", "SmallMolecule"))) {
        perturbations <- getFromNamespace(l1000perturbationsVar, 
                                          pos="package:cTRAP")
        attr(res$spearman, "perturbations") <- perturbations
        attr(res$pearson,  "perturbations") <- perturbations
        attr(res$gsea,     "perturbations") <- perturbations
    }
    
    return(res)
}