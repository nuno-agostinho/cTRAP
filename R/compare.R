#' Divide columns into chunks to be loaded into memory at a time
#' @keywords internal
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
        msg <- "%s not found: have you moved or deleted the CMap z-scores file?"
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
#'
#' @inheritParams compareAgainstReference
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

#' Prepare GSEA pathways
#'
#' @param diffExprGenes Numeric: named vector of differentially expressed genes
#'   whose names are gene identifiers and respective values are a statistic
#'   that represents significance and magnitude of differentially expressed
#'   genes (e.g. t-statistics)
#' @param geneSize Number: top and bottom number of differentially expressed
#'   genes for gene set enrichment (only used if \code{method = gsea})
#'
#' @return List of top and bottom differentially expressed genes
#' @keywords internal
prepareGSEApathways <- function(diffExprGenes, geneSize) {
    ordered         <- order(diffExprGenes, decreasing=TRUE)
    topGenes        <- names(diffExprGenes)[head(ordered, geneSize)]
    bottomGenes     <- names(diffExprGenes)[tail(ordered, geneSize)]
    pathways        <- list(topGenes, bottomGenes)
    names(pathways) <- c("top", "bottom")
    return(pathways)
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
performGSEAagainstReference <- function(diffExprGenes, reference, pathways) {
    # Calculate GSEA per data column
    gseaPerColumn <- function(k, data, pathways) {
        signature        <- data[ , k]
        names(signature) <- rownames(data)
        signature        <- sort(signature)
        score            <- fgsea(pathways=pathways, stats=signature,
                                  minSize=15, maxSize=500, nperm=1)
        return(score)
    }
    gsa <- processByChunks(reference, gseaPerColumn, pathways=pathways)

    # Calculate weighted connectivity score (WTCS) based on CMap paper (page e8)
    gsaRes <- bind_rows(gsa)
    isTop  <- gsaRes$pathway == "top"
    top    <- gsaRes[["ES"]][isTop]
    bottom <- gsaRes[["ES"]][!isTop]
    wtcs   <- ifelse(sign(top) != sign(bottom), (top - bottom) / 2, 0)

    results <- data.table("identifier"=colnames(reference), "GSEA"=wtcs)
    attr(results, "colsToRank") <- "GSEA"
    return(results)
}

#' Compare single method
#'
#' @inheritParams prepareGSEApathways
#' @param method Character: comparison methods to run (\code{spearman},
#'   \code{pearson} or \code{gsea}); multiple methods can be selected
#' @param reference Data matrix or \code{perturbationChanges} object (CMap
#'   perturbations; read \code{\link{prepareCMapPerturbations}})
#' @param cellLines Integer: number of unique cell lines
#' @param cellLineMean Boolean: add a column with the mean score across cell
#'   lines? If \code{cellLineMean = "auto"} (default) the mean score will be
#'   added if more than one cell line is available
#' @param rankByAscending Boolean: rank values based on their ascending (TRUE)
#'   or descending (FALSE) order?
#' @param rankPerCellLine Boolean: when ranking results, also rank them based on
#'   individual cell lines instead of only focusing on the mean score across
#'   cell lines; if \code{cellLineMean = FALSE}, individual cell line conditions
#'   are always ranked
#'
#' @importFrom utils head tail
#' @importFrom tidyr gather
#' @importFrom dplyr bind_rows
#' @importFrom R.utils capitalize
#'
#' @keywords internal
#' @return Data frame containing the results per method of comparison
compareAgainstReferencePerMethod <- function(method, diffExprGenes, reference,
                                             geneSize=150, cellLines=NULL,
                                             cellLineMean="auto",
                                             rankPerCellLine=FALSE) {
    startTime <- Sys.time()

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
        pathways  <- NULL
        rankedRef <- correlateAgainstReference(diffExprGenes=diffExprGenes,
                                               reference=reference,
                                               method=method)
    } else if (method == "gsea") {
        if (is.null(cellLines) || cellLines == 0) {
            msg <- compareMsg
        } else {
            cellLinesMsg <- sprintf("(%s cell line%s)...", cellLines,
                                    ifelse(cellLines == 1, "", "s"))
            msg <- paste(compareMsg, cellLinesMsg)
        }
        message(sprintf("Performing GSEA against %s...", msg))

        pathways  <- prepareGSEApathways(diffExprGenes, geneSize)
        rankedRef <- performGSEAagainstReference(diffExprGenes=diffExprGenes,
                                                 reference=reference,
                                                 pathways=pathways)
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
    attr(rankedRef, "pathways")    <- pathways

    # Report run settings and time
    diffTime <- format(round(Sys.time() - startTime, 2))
    msg      <- "Comparison against %s using '%s' %sperformed in %s\n"
    extra    <- ifelse(method == "gsea",
                       sprintf("(gene size of %s) ", geneSize), "")
    message(sprintf(msg, compareMsg, method, extra, diffTime))
    return(rankedRef)
}

#' Compare multiple methods and rank reference accordingly
#' @inheritParams compareAgainstReferencePerMethod
compareAgainstReference <- function(diffExprGenes, reference,
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

    # Summary of intersecting genes
    genes       <- intersect(names(diffExprGenes), rownames(reference))
    intersected <- length(genes)
    total       <- length(diffExprGenes)
    message(sprintf(paste(
        "Subsetting data based on %s intersecting genes",
        "(%s%% of the %s input genes)..."),
        intersected, round(intersected / total * 100, 0), total))
    if (intersected == 0) {
        stop(paste("No intersecting genes found. Check if argument",
                   "'diffExprGenes' is a named vector with gene symbols."))
    }

    names(method) <- method
    res <- lapply(method, compareAgainstReferencePerMethod,
                  diffExprGenes=diffExprGenes, reference=reference,
                  geneSize=geneSize, cellLines=cellLines,
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
    attr(ranked, "geneInfo")      <- attr(reference, "geneInfo")
    attr(ranked, "compoundInfo")  <- attr(reference, "compoundInfo")

    # Inherit from input differential expression profile
    attr(ranked, "diffExprGenes") <- diffExprGenes
    attr(ranked, "rankingInfo")   <- rankingInfo
    attr(ranked, "pathways")      <- attr(res[["gsea"]], "pathways")
    attr(ranked, "runtime")       <- Sys.time() - startTime

    class(ranked) <- c("referenceComparison", class(ranked))
    return(ranked)
}