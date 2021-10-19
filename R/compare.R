# Process a function by chunks of loaded CMap z-scores -------------------------

#' Assign columns into chunks
#'
#' @param x Vector of elements
#' @param nrows Numeric: number of rows
#' @inheritParams processByChunks
#'
#' @return List of chunks with equally distributed columns
#' @keywords internal
chunkColumns <- function(x, nrows, chunkGiB) {
    ncolsPerChunk <- ceiling( 1024^3 / (nrows * 8) )
    nchunks <- ceiling( length(x) / ncolsPerChunk )
    return( split(x, factor(sort(rank(x) %% nchunks))) )
}

processChunk <- function(chunk, data, FUN, ..., progress, verbose=FALSE) {
    if (verbose) message(Sys.time(), " - starting to load data chunk")
    if (is(data, "perturbationChanges")) {
        loaded <- loadCMapZscores(data[ , chunk], verbose=FALSE)
    } else if (is(data, "expressionDrugSensitivityAssociation")) {
        loaded <- readExpressionDrugSensitivityCorHDF5(
            data, cols=chunk, loadValues=TRUE, verbose=FALSE)
    }
    if (verbose) message(Sys.time(), " - loaded data chunk")
    res <- FUN(loaded, chunk, ..., progress=progress)
    return(res)
}

#' Process data by chunks
#'
#' @note All rows from file are currently loaded when processing chunks.
#'
#' @param data Character containing a HDF5 file path (allowing partial loading)
#'   or data matrix (processed as single chunk if data matrix)
#' @param FUN Function: function to run for each chunk
#' @param num Numeric: numbers of methods to run per chunk
#' @param ... Arguments passed to \code{FUN}
#' @inheritParams compareWithAllMethods
#'
#' @return Results of running \code{FUN}
#' @keywords internal
processByChunks <- function(data, FUN, num, ..., threads=1, chunkGiB=1,
                            verbose=FALSE) {
    loadFromFile <- is.character(data)
    if (loadFromFile && !file.exists(data)) {
        if (is(data, "perturbationChanges")) {
            type <- "z-scores"
        } else if (is(data, "expressionDrugSensitivityAssociation")) {
            type <- "gene expression and drug sensitivity association"
        }
        source <- attr(data, "source")
        msg <- "%s not found: has the %s %s file been moved or deleted?"
        stop(sprintf(msg, data, source, type))
    }

    # Display progress per chunk (if multi-threaded and on-demand file loading)
    # or per perturbation (if single-threaded)
    pb <- NULL
    if (threads == 1) pb <- startpb(max=ncol(data) * num)

    if (loadFromFile) {
        chunks <- chunkColumns(colnames(data), nrow(data), chunkGiB)
        if (threads > 1) pb <- startpb(max=length(chunks))
        resTmp <- lapply(chunks, processChunk, data, FUN, ..., threads=threads,
                         progress=pb, verbose=verbose)
        names(resTmp) <- NULL

        # Organise lists by the results of each method
        unlisted        <- unlist(resTmp, recursive=FALSE)
        len             <- sapply(unlisted, length)
        methods         <- names(unlisted)
        names(unlisted) <- NULL
        groups          <- factor(rep(methods, len), unique(methods))
        res             <- split(unlist(unlisted, recursive=FALSE), groups)
    } else {
        res <- FUN(data, colnames(data), ..., threads=threads, progress=pb)
    }
    if (!is.null(pb)) closepb(pb)
    return(res)
}

# Compare similarity of data against reference ---------------------------------

#' @importFrom stats cor.test
correlateAgainstReference <- function(k, data, diffExprGenes, method,
                                      progress) {
    res <- cor.test(data[ , k], diffExprGenes, method=method)
    if (!is.null(progress)) setpb(progress, getpb(progress) + 1)
    return(res)
}

#' @importFrom stats p.adjust
prepareCorrelationResults <- function(cors, method, pAdjust="BH") {
    cor  <- sapply(cors, "[[", "estimate")
    pval <- sapply(cors, "[[", "p.value")
    qval <- p.adjust(pval, pAdjust)
    names(cor) <- names(pval) <- names(qval) <- names(cors)

    res  <- data.table(names(cor), cor, pval, qval)
    cols <- sprintf("%s_%s", method, c("coef", "pvalue", "qvalue"))
    names(res) <- c("identifier", cols)
    attr(res, "colsToRank") <- cols[[1]]
    return(res)
}

#' Prepare GSEA gene sets
#'
#' @inheritParams compareWithAllMethods
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

#' @importFrom fastmatch fmatch
performGSEAagainstReference <- function(k, data, geneset, progress) {
    signature        <- data[ , k]
    names(signature) <- rownames(data)
    gseaParam        <- 1
    signature        <- sort(signature, decreasing=TRUE) ^ gseaParam

    filterPathways  <- function(p, stats) na.omit(fmatch(p, names(stats)))
    genesetFiltered <- lapply(geneset, filterPathways, signature)
    score <- sapply(genesetFiltered, calcGseaStat, stats=signature)
    if (!is.null(progress)) setpb(progress, getpb(progress) + 1)
    return(score)
}

prepareGSEAresults <- function(gsa) {
    gsaRes <- do.call(rbind, gsa)
    calcWTCS <- all(sort(unique(colnames(gsaRes))) == c("bottom", "top"))
    if (calcWTCS) {
        # Weighted connectivity score (WTCS) as per CMap paper (page e8)
        isTop  <- colnames(gsaRes) == "top"
        top    <- gsaRes[ , which(isTop)]
        bottom <- gsaRes[ , which(!isTop)]
        wtcs   <- ifelse(sign(top) != sign(bottom), (top - bottom) / 2, 0)
        score  <- wtcs
    } else {
        score  <- gsaRes[[1]]
    }
    if (length(score) == 0) score <- NA
    results <- data.table("identifier"=rownames(gsaRes), "GSEA"=score)
    attr(results, "colsToRank") <- "GSEA"
    return(results)
}

prepareRankedResults <- function(rankedRef, cellLineMean, cellLines, reference,
                                 rankPerCellLine, geneset) {
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
    return(rankedRef)
}

#' @importFrom parallel mclapply
comparePerMethod <- function(m, cols, reference, geneset, referenceSubset,
                             diffExprGenes, progress, threads=1) {
    if (any(m %in% c("spearman", "pearson"))) {
        cmp <- suppressWarnings(mclapply(cols, correlateAgainstReference,
                                         referenceSubset, progress=progress,
                                         diffExprGenes=diffExprGenes, method=m,
                                         mc.cores=threads))
    } else if (m %in% "gsea") {
        cmp          <- NULL
        nulls        <- sapply(geneset, is.null) # Genesets with size of 0
        validGeneset <- geneset[!nulls]
        if (length(validGeneset) > 0) {
            cmp <- mclapply(cols, performGSEAagainstReference, reference,
                            validGeneset, progress=progress, mc.cores=threads)
        }

        if (any(nulls)) {
            ns   <- names(geneset)[nulls]
            gsea <- rep(NA, length(ns))
            names(gsea) <- ns

            if (!is.null(cmp)) {
                cmp <- lapply(cmp, c, gsea)
            } else {
                cmp <- rep(list(gsea), length(cols))
                names(cmp) <- cols
            }
        }
    }
    return(cmp)
}

messageComparisonStats <- function(reference, method, cellLines, geneSize) {
    type <- attr(reference, "type")
    if (is.null(type)) type <- "comparisons"
    compareMsg <- paste(c(ncol(reference), attr(reference, "source"), type),
                        collapse=" ")

    if (is.null(cellLines) || cellLines == 0) {
        msg <- compareMsg
    } else {
        cellLinesMsg <- sprintf("(%s cell line%s)", cellLines,
                                ifelse(cellLines == 1, "", "s"))
        msg <- paste(compareMsg, cellLinesMsg)
    }
    geneSizeMsg <- ifelse("gsea" %in% method,
                          sprintf(" (gene size of %s)", geneSize), "")
    msg <- sprintf("Comparing against %s using '%s'%s...", msg,
                   paste(method, collapse=", "), geneSizeMsg)
    message(msg)
}

compareChunk <- function(reference, cols, method, diffExprGenes, genes, geneset,
                         progress, threads=1) {
    names(cols) <- cols
    reference   <- unclass(reference)
    gc()
    if (any(c("spearman", "pearson") %in% method)) {
        referenceSubset <- reference[genes, ]
    }
    names(method) <- method
    res <- lapply(method, comparePerMethod, cols, reference, geneset,
                  referenceSubset, diffExprGenes, progress, threads=threads)
    if (threads > 1 && !is.null(progress)) setpb(progress, getpb(progress) + 1)
    return(res)
}

#' Compare reference using all methods
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
#' @param method Character: comparison method (\code{spearman}, \code{pearson}
#'   or \code{gsea}; multiple methods may be selected at once)
#' @param reference Data matrix or \code{character} object with file path to
#'   CMap perturbations (see \code{\link{prepareCMapPerturbations}()}) or gene
#'   expression and drug sensitivity association (see
#'   \code{\link{loadExpressionDrugSensitivityAssociation}()})
#' @param cellLines Integer: number of unique cell lines
#' @param cellLineMean Boolean: add rows with the mean of \code{method} across
#'   cell lines? If \code{cellLineMean = "auto"} (default), rows will be added
#'   when data for more than one cell line is available.
#' @param rankByAscending Boolean: rank values based on their ascending
#'   (\code{TRUE}) or descending (\code{FALSE}) order?
#' @param rankPerCellLine Boolean: rank results based on both individual cell
#'   lines and mean scores across cell lines (\code{TRUE}) or based on mean
#'   scores alone (\code{FALSE})? If \code{cellLineMean = FALSE}, individual
#'   cell line conditions are always ranked.
#' @param threads Integer: number of parallel threads
#' @param verbose Boolean: print additional details?
#' @param chunkGiB Numeric: size (in gibibytes) of chunks to load
#'   \code{reference} file; only if argument \code{reference} is a file path
#'
#' @section GSEA score:
#'   When \code{method = "gsea"}, weighted connectivity scores (WTCS) are
#'   calculated (\url{https://clue.io/connectopedia/cmap_algorithms}).
#'
#' @importFrom utils head tail
#' @importFrom R.utils capitalize
#' @importFrom pbapply startpb getpb setpb closepb
#'
#' @return List of data tables with correlation and/or GSEA score results
#' @keywords internal
compareWithAllMethods <- function(method, input, reference, geneSize=150,
                                  cellLines=NULL, cellLineMean="auto",
                                  rankPerCellLine=FALSE, threads=1, chunkGiB=1,
                                  verbose=FALSE) {
    startTime <- Sys.time()
    geneset <- NULL
    if ("gsea" %in% method) {
        # Check if there is something wrong with the geneset(s)
        geneset  <- prepareGSEAgenesets(input, geneSize)
        geneSize <- attr(geneset, "geneSize")
    }
    genes   <- NULL
    if (any(c("spearman", "pearson") %in% method)) {
        # Subset based on intersecting genes
        genes <- intersect(names(input), rownames(reference))
        input <- input[genes]
    }
    messageComparisonStats(reference, method, cellLines, geneSize)

    rankedRef <- processByChunks(
        reference, compareChunk, length(method), method=method, verbose=verbose,
        diffExprGenes=input, genes=genes, geneset=geneset, threads=threads,
        chunkGiB=chunkGiB)
    for (m in method) {
        if (m == "gsea") {
            rankedRef[[m]] <- prepareGSEAresults(rankedRef[[m]])
        } else {
            rankedRef[[m]] <- prepareCorrelationResults(rankedRef[[m]], m)
        }
    }
    rankedRef <- lapply(rankedRef, prepareRankedResults, cellLineMean,
                        cellLines, reference, rankPerCellLine, geneset)
    # Report run settings and time
    diffTime <- format(round(Sys.time() - startTime, 2))
    message(sprintf("Comparison performed in %s\n", diffTime))
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
#' @inheritParams rankAgainstReference
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

#' Compare multiple methods and rank against reference accordingly
#' @inherit compareWithAllMethods
#' @param chunkGiB Numeric: if second argument is a path to an HDF5 file
#'   (\code{.h5} extension), that file is loaded and processed in chunks of a
#'   given size in gibibytes (GiB); lower values decrease peak RAM usage (see
#'   details below)
#'
#' @section Process data by chunks:
#'   If a file path to a valid HDF5 (\code{.h5}) file is provided instead of a
#'   data matrix, that file can be loaded and processed in chunks of size
#'   \code{chunkGiB}, resulting in decreased peak memory usage.
#'
#'   The default value of 1 GiB (1 GiB = 1024^3 bytes) allows loading chunks of ~10000 columns and
#'   14000 rows (\code{10000 * 14000 * 8 bytes / 1024^3 = 1.04 GiB}).
#'
#' @importFrom data.table :=
#' @return Data table with correlation and/or GSEA score results
#' @keywords internal
rankAgainstReference <- function(input, reference,
                                 method=c("spearman", "pearson", "gsea"),
                                 geneSize=150, cellLines=NULL,
                                 cellLineMean="auto", rankByAscending=TRUE,
                                 rankPerCellLine=FALSE, threads=1, chunkGiB=1,
                                 verbose=FALSE) {
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
    res <- compareWithAllMethods(
        method=method, input=input, reference=reference, geneSize=geneSize,
        cellLines=cellLines, cellLineMean=cellLineMean,
        rankPerCellLine=rankPerCellLine, threads=threads, chunkGiB=chunkGiB,
        verbose=verbose)

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

# Match compound identifiers between datasets ----------------------------------

filterKeys <- function(keys, cols, keyList) {
    if (is.null(keys)) keys <- keyList[keyList %in% cols]
    if (length(keys) == 0) keys <- cols[[1]]
    names(keys) <- keys
    return(keys)
}

compareDatasetIds <- function(data1, data2, key1, key2) {
    values1 <- stripStr(tolower(data1[[key1]]))
    values2 <- stripStr(tolower(data2[[key2]]))
    matches <- which(values1 %in% na.omit(values2)) # Avoid matching NAs
    return(data1[[key1]][matches])
}

getCompoundIntersectingKeyList <- function() {
    keyList       <- list()
    keyList$cmap  <- c("compound_perturbation", "pert_iname", "pert_id",
                       "smiles", "InChIKey", "pubchem_cid")
    keyList$nci60 <- c("compound", "PubChem SID", "PubChem CID", "SMILES")
    keyList$ctrp  <- c("compound", "name", "broad id", "SMILES")
    keyList$gdsc  <- c("compound", "name")
    keyList       <- unique(unlist(keyList))
    return(keyList)
}

#' Check for intersecting compounds across specific columns on both datasets
#'
#' @return List containing three elements: matching compounds
#'   \code{commonCompounds} between column \code{key 1} and \code{key 2} from
#'   the first and second datasets, respectively
#' @keywords internal
findIntersectingCompounds <- function(data1, data2, keys1=NULL, keys2=NULL) {
    showSelectedCols <- is.null(keys1) || is.null(keys2)

    # Filter keys based on dataset columns
    keyList <- getCompoundIntersectingKeyList()
    keys1 <- filterKeys(keys1, colnames(data1), keyList)
    keys2 <- filterKeys(keys2, colnames(data2), keyList)

    # Compare dataset key columns
    res <- list(key1=NULL, key2=NULL, commonCompounds=NULL)
    for (col1 in keys1) {
        for (col2 in keys2) {
            cmp <- compareDatasetIds(data1, data2, col1, col2)
            if (length(cmp) >= length(res$commonCompounds)) {
                # Save params if number of matching compounds is same or larger
                res$key1 <- col1
                res$key2 <- col2
                res$commonCompounds <- cmp
            }
        }
    }

    if (showSelectedCols) {
        message(sprintf(paste(
            "Columns '%s' and '%s' were matched based on %s common values; to",
            "manually select columns to compare, please set arguments starting",
            "with 'keyCol'"),
            res$key1, res$key2, length(res$commonCompounds)))
    }
    return(res)
}

mergeDatasets <- function(data2, data1, key2=NULL, key1=NULL,
                          suffixes=paste0(".", 1:2), ...,
                          removeKey2ColNAs=FALSE) {
    keys <- findIntersectingCompounds(data1, data2, key1, key2)
    key1 <- keys$key1
    key2 <- keys$key2

    # Convert key columns to same class if needed
    key1val <- data1[[key1]]
    key2val <- data2[[key2]]
    areNotClass <- function(key1val, key2val, cmp) cmp(key1val) && !cmp(key2val)

    FUN <- NULL
    if (length(keys$commonCompounds) > 0) {
        if (areNotClass(key1val, key2val, is.character)) {
            FUN <- as.character
        } else if (areNotClass(key1val, key2val, is.integer)) {
            FUN <- as.integer
        } else if (areNotClass(key1val, key2val, is.numeric)) {
            FUN <- as.numeric
        } else if (areNotClass(key1val, key2val, is.factor)) {
            FUN <- as.factor
        } else if (areNotClass(key1val, key2val, is.logical)) {
            FUN <- as.logical
        }
    }
    if (!is.null(FUN)) data2[[key2]] <- FUN(data2[[key2]])

    # Avoid matching NAs from key2 column of data2
    if (removeKey2ColNAs) data2 <- data2[!is.na(data2[[key2]]), ]

    # Merge data based on intersecting compounds
    data1[["matched_terms"]] <- stripStr(tolower(data1[[key1]]))
    data2[["matched_terms"]] <- stripStr(tolower(data2[[key2]]))
    df <- merge(data2, data1, by="matched_terms", suffixes=rev(suffixes), ...)
    
    id <- key1
    if (!key1 %in% colnames(df)) id <- paste0(key1, suffixes[[2]])
    df[["matched_terms"]] <- df[[id]]
    
    attr(df, "keys") <- keys
    return(df)
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
