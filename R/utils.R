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
#' # ENCODEsamples <- loadENCODEsamples(ENCODEmetadata)[[1]]
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

# cmapPerturbations object -----------------------------------------------------

#' Subset a \code{cmapPerturbations} object
#'
#' @param x \code{cmapPerturbations} object
#' @param ... Extra parameters passed to \code{`[`}
#'
#' @return \code{cmapPerturbations} object with subset data
#' @export
`[.cmapPerturbations` <- function(x, ...) {
    out <- NextMethod("[", drop=FALSE)

    # Inherit input's attributes
    attrs <- attributes(x)
    attrs$dim <- NULL
    attrs$dimnames <- NULL
    attrs$names <- NULL

    # Trim metadata to only contain subset information
    if (!is.null(ncol(out)) && ncol(x) != ncol(out) &&
        !is.null(attrs$metadata)) {
        samples <- attrs$metadata$sig_id %in% colnames(out)
        attrs$metadata <- attrs$metadata[samples, , drop=FALSE]
    }
    attributes(out) <- c(attributes(out), attrs)
    return(out)
}

#' @inherit base::as.data.frame
#' @export
as.data.frame.cmapPerturbations <- function(x, ...) NextMethod("as.data.frame")

#' @inherit utils::head
#' @export
head.cmapPerturbations <- function(x, ...) NextMethod("head")

#' @inherit utils::tail
#' @export
tail.cmapPerturbations <- function(x, ...) NextMethod("tail", ...)

# cmapComparison object --------------------------------------------------------

#' Print a \code{cmapComparison} object
#'
#' @param x \code{cmapComparison} object
#' @param perturbation Character (perturbation identifier) or numeric
#'   (perturbation index)
#' @param ... Extra parameters passed to \code{print}
#'
#' @return Information on \code{cmapPerturbations} object or on specific
#'   perturbations (if \code{perturbation} is set)
#' @export
print.cmapComparison <- function(x, perturbation=NULL, ...) {
    if (is.null(perturbation)) {
        NextMethod("print")
    } else {
        if (is.numeric(perturbation)) perturbation <- x[[1]][perturbation]

        metadata <- attr(x, "metadata")
        if (!is.null(metadata)) {
            selectMetadata <- metadata[metadata$sig_id %in% perturbation]
            if (nrow(selectMetadata) == 0) {
                # Check to see if using identifiers referring to summary stats
                summaryID <- gsub("\\_[A-Z].*\\_", "\\_", metadata$sig_id)
                selectMetadata <- metadata[summaryID %in% perturbation]
            }
        }

        compoundInfo <- attr(x, "compoundInfo")
        if (!is.null(compoundInfo)) {
            compound <- selectMetadata$pert_iname
            selectCompounds <- compoundInfo[
                compoundInfo$pert_iname %in% compound]
            res <- list(metadata=selectMetadata, compoundInfo=selectCompounds)
        } else {
            selectCompounds <- NULL
            res <- list(metadata=selectMetadata)
        }
        return(res)
    }
}

#' Cross Tabulation and Table Creation
#'
#' @param x \code{cmapComparison} object
#' @param ... Extra parameters passed to \code{table}
#' @param clean Boolean: only show certain columns (to avoid redundancy)?
#'
#' @return Complete table with metadata based on a \code{cmapComparison} object
as.table.cmapComparison <- function(x, ..., clean=TRUE) {
    metadata <- attr(x, "metadata")
    if (!is.null(metadata)) {
        nonCellID <- "non_cell_id"

        summaryID <- gsub("\\_[A-Z].*\\_", "\\_", metadata$sig_id)
        metadata[[nonCellID]] <- summaryID
        metadataSubset <- metadata[unique(match(summaryID, summaryID)), ]
        metadataSubset[ , c("cell_id", "sig_id", "distil_id")] <- NULL

        x[[nonCellID]] <- gsub("\\_[A-Z].*\\_", "\\_", x[[1]])
        res <- merge(x, metadataSubset, all.x=TRUE, by=nonCellID)
        res[[nonCellID]] <- NULL

        compoundInfo <- attr(x, "compoundInfo")
        if (!is.null(compoundInfo)) {
            res <- merge(res, compoundInfo, by="pert_iname", all.x=TRUE)
            # Place "pert_iname" column after "pert_id" one
            m <- match("pert_id", colnames(res))
            res <- res[ , c(2:m, 1, (m+1):(ncol(res))), with=FALSE]
        }
    } else {
        res <- x
    }

    if (clean) {
        hideCols <- c(colnames(res)[endsWith(colnames(res), "value") |
                                        endsWith(colnames(res), "value_rank")],
                      "pert_dose", "pert_dose_unit",
                      "pert_time", "pert_time_unit")
        hideCols <- hideCols[hideCols %in% colnames(res)]
        if (length(hideCols) > 0) res <- res[ , -hideCols, with=FALSE]
    }

    return(res)
}
