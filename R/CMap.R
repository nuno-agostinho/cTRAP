loadCMapMetadata <- function(file, nas) {
    link <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&",
                   "format=file&", "file=GSE92742_Broad_LINCS_sig_info.txt.gz")
    downloadIfNotFound(file, link)
    message(sprintf("Loading data from %s...", file))
    data <- fread(file, sep="\t", na.strings=nas)

    # Fix issues with specific metadata values
    data$pert_dose[data$pert_dose == "300.0|300.000000"] <- 300
    data$pert_dose <- as.numeric(data$pert_dose)
    data$pert_idose[data$pert_idose == "300 ng|300 ng"] <- "300 ng"
    return(data)
}

prepareCMapZscores <- function(file, zscoresID=NULL) {
    link <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&",
                   "format=file&file=GSE92742_Broad_LINCS_Level5_COMPZ.",
                   "MODZ_n473647x12328.gctx.gz")
    downloadIfNotFound(file, link, ask=TRUE)
    data <- normalizePath(file)
    attr(data, "genes")         <- readGctxIds(data, dimension="row")
    attr(data, "perturbations") <- processIds(
        zscoresID, readGctxIds(data, dimension="col"), type="cid")$ids
    return(data)
}

#' Load matrix of CMap zscores
#'
#' @param data \code{perturbationChanges} object
#' @param perturbationChanges Boolean: convert to \code{perturbationChanges}
#'   object?
#' @param verbose Boolean: print messages?
#'
#' @return Matrix containing CMap perturbation z-scores (genes as rows,
#'   perturbations as columns)
#' @export
#'
#'
#' @examples
#' \donttest{
#'   metadata <- loadCMapData("cmapMetadata.txt", "metadata")
#'   metadata <- filterCMapMetadata(metadata, cellLine="HepG2")
#'   perts <- prepareCMapPerturbations(metadata, "cmapZscores.gctx",
#'                                     "cmapGeneInfo.txt")
#'   zscores <- loadCMapZscores(perts[ , 1:10])
#' }
loadCMapZscores <- function(data, perturbationChanges=FALSE, verbose=TRUE) {
    if (verbose) message(sprintf("Loading data from %s...", data))
    zscores  <- new("GCT", src=data, cid=colnames(data), verbose=verbose)@mat
    geneInfo <- attr(data, "geneInfo")
    if (!is.null(geneInfo)) {
        rownames(zscores) <- geneInfo$pr_gene_symbol[
            match(rownames(zscores), geneInfo$pr_gene_id)]
        if (!setequal(attr(data, "genes"), rownames(zscores)))
            zscores <- zscores[attr(data, "genes"), , drop=FALSE]
    }

    if (perturbationChanges) {
        class(zscores) <- c("perturbationChanges", class(zscores))
        # Inherit input's attributes
        attrs <- attributes(data)
        attrs <- attrs[!names(attrs) %in% c(names(attributes(zscores)),
                                            "genes", "perturbations")]
        attributes(zscores) <- c(attributes(zscores), attrs)
    }
    return(zscores)
}

loadCMapCompoundInfo <- function(file, nas) {
    file <- gsub("\\_drugs|\\_samples", "", file)
    file <- sprintf("%s%s.%s", file_path_sans_ext(file),
                    c("_drugs", "_samples"), file_ext(file))
    names(file) <- c("drugs", "samples")

    readAfterComments <- function(file, comment.char="!") {
        # Ignore first rows starting with a comment character
        firstRows  <- fread(file, sep="\t", na.strings=nas, nrows=20)
        ignoreExpr <- paste0("^\\", comment.char)
        skipRows   <- min(grep(ignoreExpr, firstRows[[1]], invert=TRUE)) - 1
        data       <- fread(file, sep="\t", na.strings=nas, skip=skipRows)
        return(data)
    }

    # Process drug data
    link <- paste0(
        "https://s3.amazonaws.com/data.clue.io/repurposing/downloads/",
        "repurposing_drugs_20180907.txt")
    downloadIfNotFound(file[["drugs"]], link)

    # Replace separation symbols for targets
    message(sprintf("Loading compound data from %s...", file[["drugs"]]))
    drugData <- readAfterComments(file[["drugs"]])
    drugData$target <- gsub("|", ", ", drugData$target, fixed=TRUE)

    # Process perturbation data
    link <- paste0(
        "https://s3.amazonaws.com/data.clue.io/repurposing/downloads/",
        "repurposing_samples_20180907.txt")
    downloadIfNotFound(file[["samples"]], link)

    message(sprintf("Loading compound data from %s...", file[["samples"]]))
    pertData <- readAfterComments(file[["samples"]])
    pertData <- pertData[ , c("pert_iname", "expected_mass", "smiles",
                              "InChIKey", "pubchem_cid")]
    pertData <- unique(pertData)
    pertData <- aggregate(pertData, by=list(pertData$pert_iname),
                          function(x) paste(unique(na.omit(x)), collapse=", "))
    data <- merge(drugData, pertData, all=TRUE)
    data[data == ""] <- NA # Fix missing values
    return(data)
}

loadCMapGeneInfo <- function(file, nas) {
    link <- paste0(
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&",
        "format=file&", "file=GSE92742_Broad_LINCS_gene_info.txt.gz")
    downloadIfNotFound(file, link)
    message(sprintf("Loading data from %s...", file))
    data <- fread(file, sep="\t", na.strings=nas)
    return(data)
}

#' Load CMap data
#'
#' Load CMap data (if not found, \code{file} will be automatically downloaded)
#'
#' @note If \code{type = "compoundInfo"}, two files from
#' \strong{The Drug Repurposing Hub} will be downloaded containing information
#' about drugs and perturbations. The files will be named \code{file} with
#' \code{_drugs} and \code{_samples} before their extension, respectively.
#'
#' @param file Character: path to file
#' @param type Character: type of data to load (\code{metadata},
#' \code{geneInfo}, \code{zscores} or \code{compoundInfo})
#' @param zscoresID Character: identifiers to partially load z-scores file
#' (for performance reasons)
#'
#' @importFrom data.table fread
#' @importFrom tools file_ext file_path_sans_ext
#'
#' @return Metadata as a data table
#' @export
#'
#' @examples
#' # Load CMap metadata (data is automatically downloaded if not available)
#' cmapMetadata <- loadCMapData("cmapMetadata.txt", "metadata")
#'
#' # Load CMap gene info
#' loadCMapData("cmapGeneInfo.txt", "geneInfo")
#'
#' # Load CMap zscores based on filtered metadata
#' cmapMetadataKnockdown <- filterCMapMetadata(
#'   cmapMetadata, cellLine="HepG2",
#'   perturbationType="Consensus signature from shRNAs targeting the same gene")
#'
#' \donttest{
#' loadCMapData("cmapZscores.gctx.gz", "zscores",
#'              cmapMetadataKnockdown$sig_id)
#' }
loadCMapData <- function(file, type=c("metadata", "geneInfo", "zscores",
                                      "compoundInfo"),
                         zscoresID=NULL) {
    if (is.null(file)) stop("File cannot be NULL, please provide a filename")

    type <- match.arg(type)
    nas  <- c("", "NA", "na", "-666", "-666.0", "-666 -666",
              "-666 -666|-666 -666", "-666.000000", "-666.0|-666.000000")
    if (type == "metadata") {
        data <- loadCMapMetadata(file, nas)
    } else if (type == "geneInfo") {
        data <- loadCMapGeneInfo(file, nas)
    } else if (type == "zscores") {
        data <- prepareCMapZscores(file, zscoresID)
    } else if (type == "compoundInfo") {
        data <- loadCMapCompoundInfo(file, nas)
    }
    return(data)
}

#' List available conditions in CMap datasets
#'
#' Downloads metadata if not available
#'
#' @inheritParams filterCMapMetadata
#' @param control Boolean: show controls for perturbation types?
#'
#' @return List of conditions in CMap datasets
#' @export
#'
#' @examples
#' data("cmapMetadata")
#' # cmapMetadata <- loadCMapData("cmapMetadata.txt", "metadata")
#' getCMapConditions(cmapMetadata)
getCMapConditions <- function(metadata, cellLine=NULL, timepoint=NULL,
                              dosage=NULL, perturbationType=NULL,
                              control=FALSE) {
    metadata <- filterCMapMetadata(metadata, cellLine=cellLine,
                                   timepoint=timepoint, dosage=dosage,
                                   perturbationType=perturbationType)
    pertTypes <- getCMapPerturbationTypes(control=control)
    pertTypes <- names(pertTypes)[pertTypes %in% unique(metadata$pert_type)]

    # Order categories of value with units
    sortNumericUnitChar <- function(data, levels=NULL) {
        uniq   <- unique(na.omit(data))
        values <- as.numeric(sapply(strsplit(na.omit(uniq), " "), "[[", 1))
        units  <- sapply(strsplit(na.omit(uniq), " "), "[[", 2)
        units  <- factor(units, levels=unique(c(levels, unique(units))))
        sorted <- c(if (any(is.na(data))) NA, uniq[order(units, values)])
        return(sorted)
    }
    dose <- sortNumericUnitChar(
        metadata$pert_idose, c("%", "nM", "\u00B5M", "\u00B5L", "ng",
                               "ng/\u00B5L", "ng/mL"))
    timepoint <- sortNumericUnitChar(metadata$pert_itime)

    list("perturbationType"=pertTypes,
         "cellLine"=sort(unique(metadata$cell_id)),
         "dosage"=dose,
         "timepoint"=timepoint)
}

#' Process perturbation changes per chunk
#'
#' Perturbations will be processed per chunk if argument \code{perturbations} is
#' a file path instead of a data matrix. Otherwise, the data will be processed
#' in full (i.e., as a single chunk).
#'
#' Loading a chunk of 10000 CMap pertubations requires ~1GB of RAM.
#'
#' @inheritParams performGSEAagainstCMap
#' @param FUN Function: function to run for each chunk
#' @param ... Arguments passed to \code{FUN}
#' @param chunkSize Integer: number of perturbations to load on-demand (a higher
#'   value increases RAM usage, but also time-wise performance)
#'
#' @importFrom pbapply startpb getpb setpb closepb
#'
#' @return Results of running \code{FUN}
#' @keywords internal
processPertsPerChunk <- function(perturbations, FUN, ..., chunkSize=10000) {
    processPertCall <- function(data, perts, pertFUN, ..., progress) {
        data <- unclass(data)
        gc()

        pertFUNprogress <- function(pert, data, ..., progress, pertFUN) {
            res <- pertFUN(pert, data, ...)
            setpb(progress, getpb(progress) + 1)
            return(res)
        }

        # Suppress warnings to avoid "Cannot compute exact p-value with ties"
        res <- suppressWarnings(lapply(perts, pertFUNprogress, data, ...,
                                       progress=progress, pertFUN=pertFUN))
        # res <- plyr::compact(res) # in case of NULL elements for GSEA
        return(res)
    }

    pb <- startpb(max=ncol(perturbations))
    loadPerturbationsFromFile <- is.character(perturbations)
    if (loadPerturbationsFromFile) {
        # Divide perturbations into chunks (load into memory a chunk at a time)
        chunkVector <- function(x, nElems) {
            groups <- ceiling(length(x)/nElems)
            split(x, factor(sort(rank(x) %% groups)))
        }
        chunks <- chunkVector(colnames(perturbations), chunkSize)

        processChunk <- function(chunk, perturbations, FUN, ..., progress) {
            data <- loadCMapZscores(perturbations[ , chunk], verbose=FALSE)
            processPertCall(data, chunk, FUN, ..., progress=progress)
        }
        res <- lapply(chunks, processChunk, perturbations, FUN, ...,
                      progress=pb)
        res <- unlist(res, recursive=FALSE, use.names=FALSE)
    } else {
        res <- processPertCall(perturbations, colnames(perturbations), FUN, ...,
                               progress=pb)
    }
    closepb(pb)
    return(res)
}

#' Perform gene set enrichment (GSA) against CMap perturbations
#'
#' @inheritParams rankSimilarPerturbations
#' @inheritParams fgsea::fgsea
#'
#' @importFrom fgsea fgsea
#' @importFrom data.table data.table
#' @importFrom dplyr bind_rows
#'
#' @return Data frame containing gene set enrichment analysis (GSEA) results per
#' perturbation
#' @keywords internal
performGSEAagainstCMap <- function(diffExprGenes, perturbations, pathways,
                                   cellLines) {
    msg <- "Performing GSEA against %s CMap perturbations (%s cell lines)..."
    message(sprintf(msg, ncol(perturbations), cellLines))

    # Calculate GSEA per perturbation
    gseaPert <- function(k, data, pathways) {
        signature        <- data[ , k]
        names(signature) <- rownames(data)
        signature        <- sort(signature)
        score            <- fgsea(pathways=pathways, stats=signature,
                                  minSize=15, maxSize=500, nperm=1)
        return(score)
    }
    gsa <- processPertsPerChunk(perturbations, gseaPert, pathways=pathways)

    # Calculate weighted connectivity score (WTCS) based on CMap paper (page e8)
    gsaRes <- bind_rows(gsa)
    isTop  <- gsaRes$pathway == "top"
    top    <- gsaRes[["ES"]][isTop]
    bottom <- gsaRes[["ES"]][!isTop]
    wtcs   <- ifelse(sign(top) != sign(bottom), (top - bottom) / 2, 0)

    results <- data.table("identifier"=colnames(perturbations), "GSEA"=wtcs)
    return(results)
}

#' Correlate against CMap perturbations
#'
#' @inheritParams rankSimilarPerturbations
#' @param pAdjust Character: method to use for p-value adjustment
#'
#' @importFrom stats p.adjust cor.test
#'
#' @return Data frame with correlation results per perturbation
#'
#' @keywords internal
correlateAgainstCMap <- function(diffExprGenes, perturbations, method,
                                 cellLines, pAdjust="BH") {
    methodStr <- switch(method,
                        "spearman"="Spearman's correlation",
                        "pearson" ="Pearson's correlation",
                        "gsea"    ="GSEA")
    msg <- "Correlating against %s CMap perturbations (%s cell lines; %s)..."
    message(sprintf(msg, ncol(perturbations), cellLines, methodStr))

    # Subset based on intersecting genes
    genes <- intersect(names(diffExprGenes), rownames(perturbations))
    diffExprGenes <- diffExprGenes[genes]
    perturbations <- perturbations[genes, ]

    # Correlate per perturbation
    corPert <- function(k, data, diffExprGenes, method) {
        cor.test(data[ , k], diffExprGenes, method=method)
    }
    cors <- processPertsPerChunk(perturbations, corPert,
                                 diffExprGenes=diffExprGenes, method=method)
    cor  <- sapply(cors, "[[", "estimate")
    pval <- sapply(cors, "[[", "p.value")
    qval <- p.adjust(pval, pAdjust)
    names(cor) <- names(pval) <- names(qval) <- colnames(perturbations)

    res <- data.table(names(cor), cor, pval, qval)
    names(res) <- c("identifier", sprintf("%s_%s", method,
                                          c("coef", "pvalue", "qvalue")))
    return(res)
}

#' Prepare GSEA pathways
#'
#' @param diffExprGenes Numeric: named vector of differentially expressed genes
#'   where the name of the vector are gene names and the values are a statistic
#'   that represents significance and magnitude of differentially expressed
#'   genes (e.g. t-statistics)
#' @param geneSize Number: top and bottom number of differentially expressed
#'   genes to use for gene set enrichment (GSE) (only used if
#'   \code{method = gsea})
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

#' Calculate cell line mean
#'
#' @param data Data table: comparison against CMap data
#' @param cellLine Character: perturbation identifiers as names and respective
#' cell lines as values
#' @param metadata Data table: \code{data} metadata
#' @inheritParams compareAgainstCMapPerMethod
#'
#' @return A list with two items:
#' \describe{
#' \item{\code{data}}{input \code{data} with extra rows containing cell line
#'   average scores (if calculated)}
#' \item{\code{cellLineInfo}}{data table with cell line information}
#' \item{\code{metadata}{metadata associated with output \code{data}, including
#'   for identifiers regarding mean cell line scores}}
#' }
#' @keywords internal
calculateCellLineMean <- function(data, cellLine, metadata,
                                  rankCellLinePerturbations) {
    scoreCol <- 2
    # Remove cell line information from the identifier
    allIDs <- parseCMapID(data$identifier, cellLine=FALSE)
    idsFromMultipleCellLines <- names(table(allIDs)[table(allIDs) > 1])
    names(idsFromMultipleCellLines) <- idsFromMultipleCellLines

    # Calculate mean scores across cell lines
    calcMeanScores <- function(id, allIDs, score, cellLine) {
        matches <- id == allIDs
        list(cellLines=paste(cellLine[matches], collapse=", "),
             mean=mean(score[matches]))
    }
    res <- lapply(idsFromMultipleCellLines, calcMeanScores, allIDs=allIDs,
                  score=data[[scoreCol]], cellLine=cellLine)

    if (length(idsFromMultipleCellLines) > 0) {
        # Prepare data including for mean perturbation scores
        avg <- sapply(res, "[[", "mean")
        avgDF <- data.frame(names(avg), avg, stringsAsFactors=FALSE)
        colnames(avgDF) <- colnames(data)[c(1, scoreCol)]
        dataJoint <- bind_rows(list(data, avgDF))

        # Prepare cell line information for mean perturbation scores
        isSummarised <- allIDs %in% idsFromMultipleCellLines
        toRank <- rankCellLinePerturbations | !isSummarised
        avgCellLines <- sapply(res, "[[", "cellLines")

        cellLineInfo <- data.table(
            c(names(cellLine), names(avgCellLines)),
            c(cellLine, avgCellLines),
            c(toRank, rep(TRUE, length(avgCellLines))))

        # Append metadata associated with mean perturbation scores
        avgCellLinesMetadata <- metadata[
            match(names(avgCellLines), parseCMapID(metadata$sig_id)), ]
        avgCellLinesMetadata$sig_id <- names(avgCellLines)
        avgCellLinesMetadata$distil_id <- NA
        metadataJoint <- rbind(avgCellLinesMetadata, metadata)
    } else {
        cellLineInfo  <- data.table(names(cellLine), cellLine, TRUE)
        dataJoint     <- data
        metadataJoint <- metadata
    }
    res <- list("data"=dataJoint, "cellLineInfo"=cellLineInfo,
                "metadata"=metadataJoint)
    return(res)
}

#' Compare single method
#'
#' @inheritParams prepareGSEApathways
#' @param method Character: comparison method (\code{spearman}, \code{pearson}
#'   or \code{gsea}; multiple methods may be selected at once)
#' @param perturbations \code{perturbationChanges} object: CMap perturbations
#'   (check \code{\link{prepareCMapPerturbations}})
#' @param cellLineMean Boolean: add a column with the mean score across cell
#'   lines? If \code{cellLineMean = "auto"} (default) the mean score will be
#'   added if more than one cell line is available
#' @param rankCellLinePerturbations Boolean: when ranking results, also rank
#'   perturbations regarding individual cell lines instead of the mean
#'   perturbation score alone; if \code{cellLineMean = FALSE}, individual cell
#'   line perturbations are always ranked
#'
#' @importFrom utils head tail
#' @importFrom tidyr gather
#' @importFrom dplyr bind_rows
#'
#' @keywords internal
#' @return Data frame containing the results per method of comparison
compareAgainstCMapPerMethod <- function(
    method, diffExprGenes=diffExprGenes, perturbations=perturbations,
    geneSize=geneSize, cellLineMean=cellLineMean,
    rankCellLinePerturbations=FALSE) {

    startTime <- Sys.time()
    metadata  <- attr(perturbations, "metadata")
    cellLine  <- unique(metadata$cell_id)

    # Summary of intersecting genes
    genes <- intersect(names(diffExprGenes), rownames(perturbations))
    intersected <- length(genes)
    total       <- length(diffExprGenes)
    message(sprintf(
        "Subsetting %s intersecting genes (%s%% of the %s input genes)...",
        intersected, round(intersected / total * 100, 0), total))

    pathways <- NULL
    if (method %in% c("spearman", "pearson")) {
        data <- correlateAgainstCMap(
            diffExprGenes=diffExprGenes, perturbations=perturbations,
            method=method, cellLines=length(cellLine))
    } else if (method == "gsea") {
        pathways <- prepareGSEApathways(diffExprGenes, geneSize)
        data <- performGSEAagainstCMap(
            diffExprGenes=diffExprGenes, perturbations=perturbations,
            pathways=pathways, cellLines=length(cellLine))
    }

    # Set whether to calculate the mean value across cell lines
    if (cellLineMean == "auto") cellLineMean <- length(cellLine) > 1

    # Retrieve cell line information
    cellLine <- parseCMapID(data$identifier, cellLine=TRUE)
    names(cellLine) <- data$identifier

    if (cellLineMean) {
        aggregated <- calculateCellLineMean(data, cellLine, metadata,
                                            rankCellLinePerturbations)
        data         <- aggregated$data
        cellLineInfo <- aggregated$cellLineInfo
        metadata     <- aggregated$metadata
    } else {
        cellLineInfo <- data.table(names(cellLine), cellLine, TRUE)
    }
    names(cellLineInfo) <- c("cTRAP_id", "cellLines", paste0(method, "_rank"))
    attr(data, "cellLineInfo") <- cellLineInfo
    attr(data, "metadata") <- metadata

    # Relabel the "identifier" column name to be more descriptive
    pertType <- unique(metadata$pert_type)
    if (length(pertType) == 1) {
        pertTypes <- getCMapPerturbationTypes()
        pertType  <- names(pertTypes[pertTypes == pertType])

        if (pertType == "Compound") {
            id <- "compound_perturbation"
        } else if (grepl("biological agents", pertType)) {
            id <- "biological_agent_perturbation"
        } else {
            id <- "gene_perturbation"
        }
        colnames(data)[colnames(data) == "identifier"] <- id
    }
    rownames(data) <- data$genes

    # Add pathway information based on GSEA run if available
    if (method == "gsea" && !is.null(pathways)) {
        attr(data, "pathways") <- pathways
    }

    # Report run settings and time
    diffTime <- format(round(Sys.time() - startTime, 2))
    msg <- "Comparison against %s perturbations using '%s' %sperformed in %s\n"
    extra <- ifelse(method == "gsea",
                    sprintf("(gene size of %s) ", geneSize), "")
    message(sprintf(msg, ncol(perturbations), method, extra, diffTime))
    return(data)
}

#' Compare differential expression results against CMap perturbations
#'
#' Weighted connectivity scores (WTCS) are calculated when
#' \code{method = "gsea"}. For more information on WTCS, read
#' \url{https://clue.io/connectopedia/cmap_algorithms}.
#'
#' @inheritParams compareAgainstCMapPerMethod
#'
#' @importFrom data.table setkeyv :=
#'
#' @return Data table with correlation or GSEA results comparing differential
#' gene expression values with those associated with CMap perturbations
#' @export
#'
#' @examples
#' data("cmapPerturbationsCompounds")
#' perturbations <- cmapPerturbationsCompounds
#' data("diffExprStat")
#'
#' # Compare differential expression results against CMap perturbations
#' rankSimilarPerturbations(diffExprStat, perturbations)
#'
#' # Compare using only Spearman's correlation
#' rankSimilarPerturbations(diffExprStat, perturbations, method="spearman")
rankSimilarPerturbations <- function(diffExprGenes, perturbations,
                                     method=c("spearman", "pearson", "gsea"),
                                     geneSize=150, cellLineMean="auto",
                                     rankCellLinePerturbations=FALSE) {
    supported <- c("spearman", "pearson", "gsea")
    method <- unique(method)
    method <- method[method %in% supported]

    if (length(method) == 0) {
        stop(paste(
            "Method must contain one of the following supported comparison",
            "methods:", paste(supported, collapse=", ")))
    }

    names(method) <- method
    res <- lapply(method, compareAgainstCMapPerMethod,
                  diffExprGenes=diffExprGenes, perturbations=perturbations,
                  geneSize=geneSize, cellLineMean=cellLineMean,
                  rankCellLinePerturbations=rankCellLinePerturbations)

    # Merge metadata
    metadata <- Reduce(merge, lapply(res, attr, "metadata"))

    # Merge cell line information
    colsPerMethod <- sapply(res, length) - 1
    cellLineInfo  <- Reduce(merge, lapply(res, attr, "cellLine"))
    replaceNAsWithFALSE <- function(DT) {
        for (i in names(DT)) DT[is.na(get(i)), (i):=FALSE]
        return(DT)
    }
    cellLineInfo <- replaceNAsWithFALSE(cellLineInfo)

    pathways <- NULL
    if (!is.null(res$gsea)) pathways <- attr(res$gsea, "pathways")
    merged <- Reduce(merge, res)

    # Rank perturbations
    rankPerturbations <- function(data, cellLineInfo, colsPerMethod) {
        setkeyv(data, colnames(data)[[1]])
        for(k in seq(length(colsPerMethod))) {
            if (k == 1) {
                col <- 2
            } else {
                col <- cumsum(colsPerMethod)[k - 1] + 2
            }
            rankCol    <- paste0(names(colsPerMethod[k]), "_rank")
            toRank     <- cellLineInfo[[rankCol]]
            colsToRank <- cellLineInfo[[1]][toRank]
            ranked     <- rank(-data[colsToRank][[col]], na.last="keep")
            data[colsToRank, rankCol] <- ranked
        }
        rankCols <- grep("_rank", colnames(data))
        if (length(rankCols) > 1) {
            # Calculate rank product's rank
            ranks    <- data[colsToRank, rankCols, with=FALSE]
            rankProd <- apply(ranks, 1, prod) ^ (1 / ncol(ranks))
            data[colsToRank, "rankProduct_rank"] <- rank(rankProd,
                                                         na.last="keep")
        }
        return(data)
    }
    ranked <- rankPerturbations(merged, cellLineInfo, colsPerMethod)

    # Inherit metadata from perturbations and other useful information
    attr(ranked, "metadata")      <- metadata
    attr(ranked, "geneInfo")      <- attr(perturbations, "geneInfo")
    attr(ranked, "compoundInfo")  <- attr(perturbations, "compoundInfo")
    attr(ranked, "diffExprGenes") <- diffExprGenes
    attr(ranked, "cellLineInfo")  <- cellLineInfo
    attr(ranked, "pathways")      <- pathways

    class(ranked) <- c("similarPerturbations", class(ranked))
    return(ranked)
}

#' Filter CMap metadata
#'
#' @param metadata Data frame (CMap metadata) or character (respective filepath)
#' @param cellLine Character: cell line (if \code{NULL}, all values are loaded)
#' @param timepoint Character: timepoint (if \code{NULL}, all values are loaded)
#' @param dosage Character: dosage (if \code{NULL}, all values are loaded)
#' @param perturbationType Character: type of perturbation (if \code{NULL}, all
#' perturbation types are loaded)
#'
#' @return Filtered CMap metadata
#'
#' @export
#' @examples
#' data("cmapMetadata")
#' # cmapMetadata <- loadCMapData("cmapMetadata.txt", "metadata")
#' filterCMapMetadata(cmapMetadata, cellLine="HEPG2", timepoint="2 h",
#'                    dosage="25 ng/mL")
filterCMapMetadata <- function(metadata, cellLine=NULL, timepoint=NULL,
                               dosage=NULL, perturbationType=NULL) {
    if (is.character(metadata)) metadata <- loadCMapData(metadata, "metadata")

    filter <- list()
    if (!is.null(cellLine)) {
        metadata <- metadata[tolower(metadata$cell_id) %in% tolower(cellLine), ]
        filter$cellLine <- cellLine
    }

    if (!is.null(timepoint)) {
        metadata <- metadata[metadata$pert_itime %in% timepoint, ]
        filter$timepoint <- timepoint
    }

    if (!is.null(dosage)) {
        metadata <- metadata[metadata$pert_idose %in% dosage, ]
        filter$dosage <- dosage
    }

    if (!is.null(perturbationType)) {
        filter$perturbationType <- perturbationType
        tmp <- getCMapPerturbationTypes()[perturbationType]
        if (!is.na(tmp)) perturbationType <- tmp
        metadata <- metadata[metadata$pert_type %in% perturbationType, ]
    }
    if (length(filter) > 0) attr(metadata, "filter") <- filter

    return(metadata)
}

#' Prepare CMap perturbation data
#'
#' @param metadata Data frame (CMap metadata) or character (respective filepath
#'   to load data from file)
#' @param zscores Data frame (GCTX z-scores) or character (respective filepath
#'   to load data from file)
#' @param geneInfo Data frame (CMap gene info) or character (respective
#'   filepath to load data from file)
#' @param compoundInfo Data frame (CMap compound info) or character (respective
#'   filepath to load data from file)
#' @param loadZscores Boolean: load perturbation z-scores? Not recommended in
#'   memory-constrained systems
#'
#' @importFrom R.utils gunzip
#' @importFrom methods new
#'
#' @return CMap perturbation data attributes and filename
#' @export
#' @examples
#' \donttest{
#'   metadata <- loadCMapData("cmapMetadata.txt", "metadata")
#'   metadata <- filterCMapMetadata(metadata, cellLine="HepG2")
#'   prepareCMapPerturbations(metadata, "cmapZscores.gctx", "cmapGeneInfo.txt")
#' }
prepareCMapPerturbations <- function(metadata, zscores, geneInfo,
                                     compoundInfo=NULL, loadZscores=FALSE) {
    if (is.character(metadata)) metadata <- loadCMapData(metadata, "metadata")
    if (is.character(geneInfo)) geneInfo <- loadCMapData(geneInfo, "geneInfo")
    if (is.character(zscores)) {
        zscores <- loadCMapData(zscores, "zscores", metadata$sig_id)
    }
    if (is.character(compoundInfo)) {
        compoundInfo <- loadCMapData(compoundInfo, "compoundInfo")
    }

    attr(zscores, "genes") <- geneInfo$pr_gene_symbol[
        match(attr(zscores, "genes"), geneInfo$pr_gene_id)]
    attr(zscores, "metadata") <- metadata
    attr(zscores, "geneInfo") <- geneInfo
    attr(zscores, "compoundInfo") <- compoundInfo
    class(zscores) <- c("perturbationChanges", class(zscores))

    if (loadZscores) zscores <- loadCMapZscores(zscores,
                                                perturbationChanges=TRUE)

    # Display summary message of loaded perturbations
    filters <- attr(metadata, "filter")
    summaryMsg <- sprintf(
        "\nSummary: %s perturbations measured across %s genes",
        ncol(zscores), nrow(zscores))
    if (!is.null(filters)) {
        filterNames <- c("cellLine"="Cell lines",
                         "perturbationType"="Perturbation types",
                         "dosage"="Perturbation doses",
                         "timepoint"="Time points")
        filterNames <- filterNames[names(filters)]

        filterMsg <- paste0("  - ", filterNames, ": ",
                            sapply(filters, paste, collapse=", "),
                            collapse="\n")
        summaryMsg <- sprintf("%s filtered by:\n%s", summaryMsg, filterMsg)
    }
    message(summaryMsg)
    return(zscores)
}

#' Get perturbation types
#'
#' @param control Boolean: return perturbation types used as control?
#'
#' @return Perturbation types and respective codes as used by CMap datasets
#' @export
#'
#' @examples
#' getCMapPerturbationTypes()
getCMapPerturbationTypes <- function (control=FALSE) {
    perts <- c("Compound"="trt_cp",
      "Peptides and other biological agents (e.g. cytokine)"="trt_lig",
      "shRNA for loss of function (LoF) of gene"="trt_sh",
      "Consensus signature from shRNAs targeting the same gene"="trt_sh.cgs",
      "cDNA for overexpression of wild-type gene"="trt_oe",
      "cDNA for overexpression of mutated gene"="trt_oe.mut",
      "CRISPR for LLoF"="trt_xpr")

    if (control) {
        controlPerts <- c("ctl_vehicle", "ctl_vector", "trt_sh.css",
                          "ctl_vehicle.cns", "ctl_vector.cns", "ctl_untrt.cns",
                          "ctl_untrt")
        names(controlPerts) <- c(
            "vehicle for compound treatment (e.g DMSO)",
            "vector for genetic perturbation (e.g empty vector, GFP)",
            "consensus signature from shRNAs that share a common seed sequence",
            "consensus signature of vehicles",
            "consensus signature of vectors",
            "consensus signature of many untreated wells",
            "Untreated cells")
        names(controlPerts) <- paste("Controls -", names(controlPerts))

        res <- c(perts, controlPerts)
    } else {
        res <- perts
    }
    return(res)
}
