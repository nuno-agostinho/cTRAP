#' Assign vector elements into chunks
#'
#' @param x Vector of elements
#' @param nElems Numeric: number of chunks
#'
#' @keywords internal
#'
#' @return List of chunks with the original vector elements divided
chunkVector <- function(x, nElems) {
    groups <- ceiling(length(x)/nElems)
    split(x, factor(sort(rank(x) %% groups)))
}

processChunk <- function(chunk, data, FUN, ..., progress) {
    zscores <- loadCMapZscores(data[ , chunk], verbose=FALSE)
    processCall(zscores, chunk, FUN, ..., progress=progress)
}

processCall <- function(data, cols, calledFUN, ..., progress) {
    data <- unclass(data)
    gc()

    setFUNprogress <- function(col, data, ..., progress, calledFUN) {
        res <- calledFUN(col, data, ...)
        setpb(progress, getpb(progress) + 1)
        return(res)
    }

    # Suppress warnings to avoid "Cannot compute exact p-value with ties"
    res <- suppressWarnings(lapply(cols, setFUNprogress, data, ...,
                                   progress=progress, calledFUN=calledFUN))
    # res <- plyr::compact(res) # in case of NULL elements for GSEA
    return(res)
}

#' Process data column by chunks
#'
#' Columns will be processed per chunk if argument \code{data} is a
#' \code{perturbationChanges} object containing a file path instead of a data
#' matrix. Otherwise, the data will be processed as a single chunk.
#'
#' For instance, loading a chunk of 10000 CMap pertubations requires ~1GB of RAM
#' compared to loading the whole dataset.
#'
#' @param data Data matrix or \code{perturbationChanges} object
#' @param FUN Function: function to run for each chunk
#' @param ... Arguments passed to \code{FUN}
#' @param chunkSize Integer: number of columns to load on-demand (a higher value
#'   increases RAM usage, but decreases running time)
#'
#' @importFrom pbapply startpb getpb setpb closepb
#'
#' @return Results of running \code{FUN}
#' @keywords internal
processByChunks <- function(data, FUN, ..., chunkSize=10000) {
    pb <- startpb(max=ncol(data))
    loadFromFile <- is.character(data)
    if (loadFromFile && !file.exists(data)) {
        msg <- "%s not found: has the CMap z-scores file been moved or deleted?"
        stop(sprintf(msg, data))
    }

    if (loadFromFile) {
        chunks <- chunkVector(colnames(data), chunkSize)
        res <- lapply(chunks, processChunk, data, FUN, ..., progress=pb)
        res <- unlist(res, recursive=FALSE, use.names=FALSE)
    } else {
        res <- processCall(data, colnames(data), FUN, ..., progress=pb)
    }
    closepb(pb)
    return(res)
}

#' Rank columns in a dataset
#'
#' @details The rank product's rank is calculated if more than one method is
#'   ranked.
#'
#' @note The first column of \code{data} and \code{rankingInfo} must contain
#'   common identifiers.
#'
#' @param table Data table: data; first column must be identifiers
#' @param rankingInfo Data table: boolean values of which rows to rank based on
#'   columns (column names to be ranked must exactly match those available in
#'   argument \code{table}); first column must be identifiers
#' @param sort Boolean: sort data based on rank product's rank (if multiple
#'   methods are available) or by available ranks
#' @inheritParams compareAgainstReference
#'
#' @importFrom data.table setkeyv
#'
#' @return Data table with the contents of \code{table} and extra columns with
#'   respective rankings
#' @keywords internal
rankColumns <- function(table, rankingInfo, rankByAscending=TRUE, sort=FALSE) {
    setkeyv(table, colnames(table)[[1]])
    colsToRank <- colnames(rankingInfo)[-1]
    rankedCols <- NULL
    for(col in colsToRank) {
        toRank     <- rankingInfo[[col]]
        rowsToRank <- rankingInfo[[1]][toRank]
        dataToRank <- table[rowsToRank][[col]]
        if (rankByAscending) dataToRank <- -dataToRank
        ranked     <- rank(dataToRank, na.last="keep")
        newCol     <- paste0(gsub("(.*)_.*$", "\\1", col), "_rank")
        table[rowsToRank, newCol] <- ranked
        rankedCols <- c(rankedCols, newCol)
    }
    if (length(rankedCols) > 1) {
        # Calculate rank product's rank
        ranks    <- table[rowsToRank, rankedCols, with=FALSE]
        rankProd <- apply(ranks, 1, prod) ^ (1 / ncol(ranks))
        table[rowsToRank, "rankProduct_rank"] <- rank(rankProd, na.last="keep")
        sortingCol <- "rankProduct_rank"
    } else {
        sortingCol <- rankedCols
    }
    if (sort) table <- table[order(table[[sortingCol]])]
    return(table)
}

#' Correlate against data columns
#'
#' @inheritParams rankSimilarPerturbations
#' @param pAdjust Character: method to use for p-value adjustment
#'
#' @importFrom stats p.adjust cor.test
#' @importFrom data.table data.table
#'
#' @return Data frame with correlation results per data column
#'
#' @keywords internal
correlateAgainstReference <- function(diffExprGenes, reference, method,
                                      pAdjust="BH") {
    # Subset based on intersecting genes
    genes <- intersect(names(diffExprGenes), rownames(reference))
    diffExprGenes <- diffExprGenes[genes]
    reference     <- reference[genes, ]

    # Correlate per data column
    corPerColumn <- function(k, data, diffExprGenes, method) {
        cor.test(data[ , k], diffExprGenes, method=method)
    }
    cors <- processByChunks(reference, corPerColumn,
                            diffExprGenes=diffExprGenes, method=method)
    cor  <- sapply(cors, "[[", "estimate")
    pval <- sapply(cors, "[[", "p.value")
    qval <- p.adjust(pval, pAdjust)
    names(cor) <- names(pval) <- names(qval) <- colnames(reference)

    res  <- data.table(names(cor), cor, pval, qval)
    cols <- sprintf("%s_%s", method, c("coef", "pvalue", "qvalue"))
    names(res) <- c("identifier", cols)
    attr(res, "colsToRank") <- cols[[1]]
    return(res)
}

#' Prepare GSEA gene sets
#'
#' @inheritParams compareAgainstReferencePerMethod
#' @return List of gene sets
#' @keywords internal
prepareGSEAgenesets <- function(input, geneSize) {
    input <- prepareGeneInput(input)
    isGeneset <- isTRUE(attr(input, "isGeneset"))
    if(isGeneset) {
        geneset <- list("custom"=input)
        geneSize <- c(length(input), 0)
    } else {
        # Check if length of input is too small for geneSize
        genes <- names(input)
        if (length(geneSize) == 1) {
            isGeneNumberLow <- floor(length(genes)/2) < geneSize
        } else if (length(geneSize) == 2) {
            isGeneNumberLow <- length(genes) < sum(geneSize)
        }

        if (isGeneNumberLow) {
            geneSize <- floor(length(genes)/2)
            msg <- paste(
                "'input' contains a smaller number of genes than available for",
                "'geneSize'; 'geneSize' parameter was reduced to", geneSize)
            warning(msg, immediate.=TRUE)
        }
        if (length(geneSize) == 1) geneSize <- rep(geneSize, 2)

        ordered  <- order(input, decreasing=TRUE)
        topGenes <- genes[head(ordered, geneSize[[1]])]
        if (geneSize[[1]] <= 0) topGenes <- NULL

        bottomGenes <- genes[tail(ordered, geneSize[[2]])]
        if (geneSize[[2]] <= 0) bottomGenes <- NULL

        geneset <- list("top"=topGenes, "bottom"=bottomGenes)
    }
    attr(geneset, "geneSize") <- unique(geneSize)
    return(geneset)
}

#' Perform gene set enrichment (GSA) against data columns
#'
#' @inheritParams rankSimilarPerturbations
#' @inheritParams fgsea::fgsea
#'
#' @importFrom fgsea fgsea
#' @importFrom data.table data.table
#' @importFrom dplyr bind_rows
#'
#' @return Data frame containing gene set enrichment analysis (GSEA) results per
#' data column
#' @keywords internal
performGSEAagainstReference <- function(reference, geneset) {
    # Calculate GSEA per data column
    gseaPerColumn <- function(k, data, geneset) {
        signature        <- data[ , k]
        names(signature) <- rownames(data)
        signature        <- sort(signature)
        score            <- fgsea(pathways=geneset, stats=signature,
                                  minSize=15, maxSize=500, nperm=1)
        return(score)
    }
    gsa <- processByChunks(reference, gseaPerColumn, geneset)
    gsaRes <- bind_rows(gsa)

    calcWTCS <- all(sort(unique(gsaRes$pathway)) == c("bottom", "top"))
    if (calcWTCS) {
        # Weighted connectivity score (WTCS) as per CMap paper (page e8)
        isTop  <- gsaRes$pathway == "top"
        top    <- gsaRes[["ES"]][isTop]
        bottom <- gsaRes[["ES"]][!isTop]
        wtcs   <- ifelse(sign(top) != sign(bottom), (top - bottom) / 2, 0)
        score  <- wtcs
    } else {
        score  <- gsaRes[["ES"]]
    }
    if (length(score) == 0) score <- NA
    results <- data.table("identifier"=colnames(reference), "GSEA"=score)
    attr(results, "colsToRank") <- "GSEA"
    return(results)
}

#' Compare single method
#'
#' @param input \code{Named numeric vector} of differentially expressed genes
#'   whose names are gene identifiers and respective values are a statistic that
#'   represents significance and magnitude of differentially expressed genes
#'   (e.g. t-statistics); or \code{character} of gene symbols composing a gene
#'   set that is tested for enrichment in reference data (only used if
#'   \code{method} includes \code{gsea})
#' @param geneSize Numeric: number of top up-/down-regulated genes to use as
#'   gene sets to test for enrichment in reference data; if a 2-length numeric
#'   vector, the first index is the number of top up-regulated genes and the
#'   second index is the number of down-regulated genes used to create gene
#'   sets; only used if \code{method} includes \code{gsea} and if \code{input}
#'   is not a gene set
#' @param method Character: one or more methods to compare data
#'   (\code{spearman}, \code{pearson} or \code{gsea})
#' @param reference Data matrix or \code{perturbationChanges} object (CMap
#'   perturbations; see \code{\link{prepareCMapPerturbations}()})
#' @param cellLines Integer: number of unique cell lines
#' @param cellLineMean Boolean: add a column with the mean score across cell
#'   lines? If \code{cellLineMean = "auto"} (default), the mean score will be
#'   added when data for more than one cell line is available.
#' @param rankByAscending Boolean: rank values based on their ascending
#'   (\code{TRUE}) or descending (\code{FALSE}) order?
#' @param rankPerCellLine Boolean: rank results based on both individual cell
#'   lines and mean scores across cell lines (\code{TRUE}) or based on mean
#'   scores alone (\code{FALSE})? If \code{cellLineMean = FALSE}, individual
#'   cell line conditions are always ranked.
#'
#' @importFrom utils head tail
#' @importFrom R.utils capitalize
#'
#' @keywords internal
#' @return Data frame containing the results per method of comparison
compareAgainstReferencePerMethod <- function(method, input, reference,
                                             geneSize=150, cellLines=NULL,
                                             cellLineMean="auto",
                                             rankPerCellLine=FALSE) {
    startTime <- Sys.time()

    # Check immediately if there is something wrong with the geneset(s)
    if (method == "gsea") geneset <- prepareGSEAgenesets(input, geneSize)

    type <- attr(reference, "type")
    if (is.null(type)) type <- "comparisons"
    compareMsg <- paste(c(ncol(reference),
                          attr(reference, "source"), type), collapse=" ")

    if (method %in% c("spearman", "pearson")) {
        methodStr <- paste0(capitalize(method), "'s correlation")
        if (is.null(cellLines) || cellLines == 0) {
            msg <- methodStr
        } else {
            msg <- sprintf("%s cell line%s; %s", cellLines,
                           ifelse(cellLines == 1, "", "s"), methodStr)
        }
        message(sprintf("Correlating against %s (%s)...", compareMsg, msg))
        geneset  <- NULL
        rankedRef <- correlateAgainstReference(input, reference, method)
    } else if (method == "gsea") {
        if (is.null(cellLines) || cellLines == 0) {
            msg <- compareMsg
        } else {
            cellLinesMsg <- sprintf("(%s cell line%s)...", cellLines,
                                    ifelse(cellLines == 1, "", "s"))
            msg <- paste(compareMsg, cellLinesMsg)
        }
        message(sprintf("Performing GSEA against %s...", msg))
        geneSize  <- attr(geneset, "geneSize")
        rankedRef <- performGSEAagainstReference(reference, geneset)
    }
    colsToRank <- attr(rankedRef, "colsToRank")

    # Set whether to calculate the mean value across cell lines
    if (cellLineMean == "auto") cellLineMean <- cellLines > 1

    # Retrieve ranking information
    if (is(reference, "perturbationChanges")) {
        cellLine <- parseCMapID(rankedRef[["identifier"]], cellLine=TRUE)
        names(cellLine) <- rankedRef[["identifier"]]
    } else {
        cellLine <- rankedRef[["identifier"]]
        names(cellLine) <- cellLine
    }

    metadata <- attr(reference, "metadata")
    if (cellLineMean) {
        aggregated  <- calculateCellLineMean(rankedRef, cellLine, metadata,
                                             rankPerCellLine)
        rankedRef   <- aggregated$reference
        rankingInfo <- aggregated$rankingInfo
        metadata    <- aggregated$metadata
    } else {
        rankingInfo <- data.table(names(cellLine), TRUE)
    }
    names(rankingInfo) <- c("cTRAP_id", colsToRank)

    # Inherit information of interest
    attr(rankedRef, "rankingInfo") <- rankingInfo
    attr(rankedRef, "metadata")    <- metadata
    attr(rankedRef, "geneset")     <- geneset

    # Report run settings and time
    diffTime <- format(round(Sys.time() - startTime, 2))
    msg      <- "Comparison against %s using '%s' %sperformed in %s\n"
    extra    <- ifelse(method == "gsea",
                       sprintf("(gene size of %s) ", geneSize), "")
    message(sprintf(msg, compareMsg, method, extra, diffTime))
    return(rankedRef)
}

prepareGeneInput <- function(input) {
    if (is.null(names(input)) && is.character(input)) {
        isGeneset   <- TRUE
        geneSymbols <- input
    } else if (!is.null(names(input)) && is.numeric(input)) {
        isGeneset   <- FALSE
        geneSymbols <- names(input)
    } else {
        stop("argument 'input' must be a named numeric vector or a character",
             " vector with gene symbols.")
    }
    attr(input, "isGeneset")   <- isGeneset
    attr(input, "geneSymbols") <- geneSymbols
    return(input)
}

#' Compare multiple methods and rank reference accordingly
#'
#' @inheritParams compareAgainstReferencePerMethod
#'
#' @importFrom data.table :=
#'
#' @keywords internal
#' @return List of data frame containing the results per methods of comparison
compareAgainstReference <- function(input, reference,
                                    method=c("spearman", "pearson", "gsea"),
                                    geneSize=150, cellLines=NULL,
                                    cellLineMean="auto", rankByAscending=TRUE,
                                    rankPerCellLine=FALSE) {
    startTime <- Sys.time()
    # Check if any of supplied methods are supported
    supported <- c("spearman", "pearson", "gsea")
    method <- unique(method)
    method <- method[method %in% supported]
    if (length(method) == 0) {
        stop(paste(
            "Argument 'method' must contain one of the following supported",
            "comparison methods:", paste(supported, collapse=", ")))
    }
    input <- prepareGeneInput(input)
    if (attr(input, "isGeneset")) {
        if (!"gsea" %in% method) {
            msg <- paste("Method 'gsea' is automatically performed if argument",
                         "'input' is a gene set.")
            warning(msg)
        }
        method <- "gsea"
    }

    # Summary of intersecting genes
    genes       <- intersect(attr(input, "geneSymbols"), rownames(reference))
    intersected <- length(genes)
    total       <- length(input)
    message(sprintf(paste(
        "Subsetting data based on %s intersecting genes",
        "(%s%% of the %s input genes)..."),
        intersected, round(intersected / total * 100, 0), total))
    if (intersected == 0) {
        stop("No intersecting genes found. Check if argument 'input' is a",
             " named numeric vector or a character vector with gene symbols.")
    }

    names(method) <- method
    res <- lapply(method, compareAgainstReferencePerMethod, input=input,
                  reference=reference, geneSize=geneSize, cellLines=cellLines,
                  cellLineMean=cellLineMean, rankPerCellLine=rankPerCellLine)

    # Rank columns
    rankingInfo <- Reduce(merge, lapply(res, attr, "rankingInfo"))
    replaceNAsWithFALSE <- function(DT) {
        for (i in names(DT)) DT[is.na(get(i)), (i):=FALSE]
        return(DT)
    }
    rankingInfo <- replaceNAsWithFALSE(rankingInfo)
    merged <- Reduce(merge, res)
    ranked <- rankColumns(merged, rankingInfo, rankByAscending=rankByAscending,
                          sort=TRUE)

    # Inherit metadata
    attr(ranked, "metadata") <- Reduce(merge, lapply(res, attr, "metadata"))
    attr(ranked, "geneInfo")        <- attr(reference, "geneInfo")
    attr(ranked, "compoundInfo")    <- attr(reference, "compoundInfo")
    attr(ranked, "zscoresFilename") <- attr(reference, "zscoresFilename")

    # Inherit from input
    attr(ranked, "input")         <- input
    attr(ranked, "rankingInfo")   <- rankingInfo
    attr(ranked, "geneset")       <- c(attr(res[["gsea"]], "pathways"), # legacy
                                       attr(res[["gsea"]], "geneset"))
    attr(ranked, "runtime")       <- Sys.time() - startTime

    class(ranked) <- c("referenceComparison", class(ranked))
    return(ranked)
}

# referenceComparison object ---------------------------------------------------

convertToTable <- function(x, clean=TRUE) {
    metadata <- attr(x, "metadata")
    isMetadataUseful <- !is.null(metadata) && nrow(metadata) > 0
    if (isMetadataUseful) {
        nonCellID <- "non_cell_id"

        summaryID <- parseCMapID(metadata$sig_id, cellLine=FALSE)
        metadata[[nonCellID]] <- summaryID
        metadataSubset <- metadata[unique(match(summaryID, summaryID)), ]
        metadataSubset[ , c("sig_id", "distil_id")] <- NULL

        x[[nonCellID]] <- parseCMapID(x[[1]], cellLine=FALSE)
        res <- merge(x, metadataSubset, all.x=TRUE, by=nonCellID)
        res[[nonCellID]] <- NULL
    } else {
        res <- x
    }

    geneInfo <- attr(x, "geneInfo")
    isGeneInfoUseful <- !is.null(geneInfo) && any(
        res[["pert_iname"]] %in% geneInfo[["pr_gene_symbol"]])
    if (isGeneInfoUseful) {
        res <- merge(res, geneInfo, by.x="pert_iname",
                     by.y="pr_gene_symbol", all.x=TRUE)
        # Place "pert_iname" column after "pert_id" one
        m <- match("pert_id", colnames(res))
        res <- res[ , c(2:m, 1, (m+1):(ncol(res))), with=FALSE]
    }

    compoundInfo <- attr(x, "compoundInfo")
    if (is(x, "similarPerturbations")) {
        compoundCol1 <- compoundCol2 <- "pert_iname"
    } else {
        compoundCol1 <- "compound"
        compoundCol2 <- names(compoundInfo)[[1]]
    }
    isCompoundInfoUseful <- !is.null(compoundInfo) && any(
        res[[compoundCol1]] %in% compoundInfo[[compoundCol2]])
    if (isCompoundInfoUseful) {
        res[[compoundCol1]] <- as.character(res[[compoundCol1]])
        compoundInfo[[compoundCol2]] <- as.character(
            compoundInfo[[compoundCol2]])
        res <- merge(res, compoundInfo, by.x=compoundCol1, by.y=compoundCol2,
                     all.x=TRUE)
        pos <- match("pert_id", colnames(res))
        if (!is.na(pos)) {
            # Place "pert_iname" column after "pert_id"
            res <- res[ , c(2:pos, 1, (pos + 1):(ncol(res))), with=FALSE]
        }
    }
    if (clean) {
        hideCols <- c(colnames(res)[endsWith(colnames(res), "value") |
                                        endsWith(colnames(res), "value_rank")],
                      "pert_dose", "pert_dose_unit",
                      "pert_time", "pert_time_unit")
        hideCols <- hideCols[hideCols %in% colnames(res)]
        if (length(hideCols) > 0) res <- res[ , -hideCols, with=FALSE]
    }
    res <- data.table(res)

    # Reorder table based on original ordering
    if (all(x[[1]] %in% res[[1]])) {
        res <- res[match(x[[1]], res[[1]]), , drop=FALSE]
    }
    return(res)
}

#' Cross Tabulation and Table Creation
#'
#' @param x \code{referenceComparison} object
#' @param ... Extra parameters not currently used
#' @param clean Boolean: only show certain columns (to avoid redundancy)?
#'
#' @family functions related with the ranking of CMap perturbations
#' @family functions related with the prediction of targeting drugs
#' @return Complete table with metadata based on a \code{targetingDrugs} object
#' @export
as.table.referenceComparison <- function(x, ..., clean=TRUE) {
    return(convertToTable(x, clean=clean))
}
