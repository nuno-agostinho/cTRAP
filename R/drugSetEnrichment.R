#' Load table with drug descriptors
#'
#' @param source Character: molecular descriptors for compounds in \code{NCI60}
#'   or \code{CMap}
#' @param type Character: load \code{2D} or \code{3D} molecular descriptors
#' @param file Character: filepath to drug descriptors (automatically downloaded
#'   if file does not exist)
#'
#' @family functions for drug set enrichment analysis
#' @return Data table with drug descriptors
#' @export
#'
#' @examples
#' loadDrugDescriptors()
loadDrugDescriptors <- function(source=c("NCI60", "CMap"), type=c("2D", "3D"),
                                file=NULL) {
    source <- match.arg(source)
    type   <- match.arg(type)
    if (source == "NCI60" && type == "2D") {
        link <- "599ok2w9ahysdga/compound_descriptors_NCI60_2D.rds"
    } else if (source == "NCI60" && type == "3D") {
        link <- "c2hbmk8qi3tyrh4/compound_descriptors_NCI60_3D.rds"
    } else if (source == "CMap" && type == "2D") {
        link <- "u1ath10e753x6en/compound_descriptors_CMap_2D.rds"
    } else if (source == "CMap" && type == "3D") {
        link <- "tpu3sq53mpy5fvt/compound_descriptors_CMap_3D.rds"
    } else {
        stop("selected 'source' and 'type' are not supported")
    }

    link  <- sprintf("https://www.dropbox.com/s/%s?raw=1", link)
    if (is.null(file))
        file  <- sprintf("molecular_descriptors_%s_%s.rds", source, type)
    file  <- downloadIfNotFound(link, file)
    table <- readRDS(file)
    return(table)
}

#' Calculate evenly-distributed bins
#'
#' @param numbers Numeric
#' @param maxBins Numeric: maximum number of bins for numeric columns
#' @param k Numeric: constant; the higher the constant, the smaller the bin size
#'   (check \code{minpts})
#' @param minPoints Numeric: minimum number of points in a bin (if \code{NULL},
#'   the minimum number of points is the number of non-missing values divided by
#'   \code{maxBins} divided by \code{k})
#' @inheritDotParams binr::bins -x -target.bins -minpts
#'
#' @importFrom binr bins bins.getvals
#'
#' @return Factor containing the respective group of each element in
#'   \code{numbers}
#' @keywords internal
calculateEvenlyDistributedBins <- function(numbers, maxBins=15, k=5,
                                           minPoints=NULL, ..., ids=NULL) {
    nas     <- is.na(numbers)
    if (all(nas)) return(NULL)

    numbers <- round(numbers[!nas])
    if (max(numbers) - min(numbers) == 0) return(setNames(numbers, ids[!nas]))
    if (is.null(minPoints)) {
        minPoints <- round( length(numbers) / maxBins / k )
        if (minPoints == 0) return(NULL)
    }
    bin <- bins(numbers, target.bins=maxBins, minpts=minPoints, ...)

    # Replace labels of single number intervals
    names(bin$binct) <- gsub("\\[([-]?[0-9]*), \\1\\]", "\\1", names(bin$binct))
    # Suppress warnings to avoid integer64 overflow warnings
    factors <- suppressWarnings(
        cut(numbers, bins.getvals(bin), labels=names(bin$binct)))
    if (!is.null(ids)) names(factors) <- ids[!nas]
    return(factors)
}

#' Prepare drug sets from a table with compound descriptors
#'
#' Create a list of drug sets for each character and numeric column. For each
#' character column, drugs are split across that column's unique values (see
#' argument \code{maxUniqueElems}). For each numeric column, drugs are split
#' across evenly-distributed bins.
#'
#' @param table Data frame: drug descriptors
#' @param id Integer or character: index or name of the identifier column
#' @param maxUniqueElems Numeric: ignore character columns with more unique
#'   elements than \code{maxUniqueElems}
#' @inheritParams calculateEvenlyDistributedBins
#'
#' @importFrom pbapply pblapply
#'
#' @family functions for drug set enrichment analysis
#' @return Named list of characters: named drug sets with respective compound
#'   identifiers as list elements
#' @export
#'
#' @examples
#' descriptors <- loadDrugDescriptors("NCI60")
#' prepareDrugSets(descriptors)
prepareDrugSets <- function(table, id=1, maxUniqueElems=15, maxBins=15, k=5,
                            minPoints=NULL) {
    # Remove elements with no ID
    valid <- !is.na(table[[id]])
    table <- table[valid, ]

    # Prepare sets from character columns if # unique values <= maxUniqueElems
    isCharacter <- sapply(table, is.character)
    uniqueElems <- sapply(lapply(table[ , isCharacter, with=FALSE], unique),
                          length)
    sets <- names(uniqueElems[uniqueElems <= maxUniqueElems])
    subtable <- table[ , sets, with=FALSE]

    res <- lapply(subtable, function(x, ids) split(ids, x, drop=TRUE),
                  ids=table[[id]])

    nonCharacterTable <- table[ , !isCharacter, with=FALSE]
    res2 <- pblapply(nonCharacterTable, calculateEvenlyDistributedBins,
                     ids=table[[id]], maxBins=maxBins, k=k, minPoints=minPoints)
    res2 <- Filter(length, res2)
    res2 <- lapply(res2, function(x) split(names(x), x, drop=TRUE))

    res <- c(res, res2)
    symbol <- ": "
    names(res) <- paste0(names(res), symbol)
    res <- unlist(res, recursive=FALSE)
    names(res) <- gsub(paste0(symbol, "."), symbol, names(res), fixed=TRUE)

    # Inherit input attributes
    attributes(res) <- c(attributes(res),
                         attributes(table)[c("compoundInfo", "source", "type")])
    return(res)
}

#' Prepare stats' compound information
#' @keywords internal
prepareStatsCompoundInfo <- function(stats) {
    statsCompoundInfo <- attr(stats, "compoundInfo")
    if (is.vector(stats)) {
        statsInfo  <- data.table("id"=names(stats), "values"=stats)
        statsIDcol <- colnames(statsInfo)[[1]]
    } else if (is(stats, "referenceComparison")) {
        statsInfo  <- as.table(stats, clean=FALSE)
        statsIDcol <- colnames(statsInfo)[[1]]
    } else if (!is.null(statsCompoundInfo)) {
        if (is(stats, "data.table") && !is(statsCompoundInfo, "data.table")) {
            statsCompoundInfo <- data.table(statsCompoundInfo)
        }
        statsIDcol                <- colnames(stats)[[1]]
        stats[[statsIDcol]]       <- as.character(stats[[statsIDcol]])
        statsCompoundInfo[["id"]] <- as.character(statsCompoundInfo[["id"]])
        statsInfo <- merge(stats, statsCompoundInfo, by.x=statsIDcol, by.y="id",
                           all.x=TRUE)
    } else {
        msg <- paste(
            "Argument 'stats' needs to be a vector or a 'referenceComparison'",
            "object or have an attribute called 'compoundInfo'")
        stop(msg)
    }
    return(list("statsInfo"=statsInfo, "statsIDcol"=statsIDcol))
}

#' Get drug sets' compound info
#' @keywords internal
prepareSetsCompoundInfo <- function(sets) {
    setsCompoundInfo <- attr(sets, "compoundInfo")
    setsIDcol <- "id"
    if (is.null(setsCompoundInfo)) {
        setsCompoundInfo <- data.table("id"=unique(unlist(sets)))
    } else if (!setsIDcol %in% colnames(setsCompoundInfo)) {
        setsIDcol <- colnames(setsCompoundInfo)[[1]]
    }
    return(list("setsCompoundInfo"=setsCompoundInfo, "setsIDcol"=setsIDcol))
}

#' Match identifiers between data and drug sets
#'
#' @inheritParams analyseDrugSetEnrichment
#' @importFrom data.table data.table
#' @keywords internal
#'
#' @return Statistic values from input data and corresponding identifiers as
#'   names (if no match is found, the original identifier from argument
#'   \code{stats} is used)
matchStatsWithDrugSetsID <- function(sets, stats, col="values",
                                     keyColSets=NULL, keyColStats=NULL) {
    res        <- prepareStatsCompoundInfo(stats)
    statsInfo  <- res$statsInfo
    statsIDcol <- res$statsIDcol

    res              <- prepareSetsCompoundInfo(sets)
    setsCompoundInfo <- res$setsCompoundInfo
    setsIDcol        <- res$setsIDcol

    checkIfIDwasReplacedAfterMerging <- function(statsIDcol, data) {
        keys <- attr(data, "keys")
        if (keys$key1 == statsIDcol) statsIDcol <- keys$key2
        return(statsIDcol)
    }

    # Return statistical values with corresponding identifier (or original
    # identifier if no match is found)
    df  <- mergeDatasets(setsCompoundInfo, statsInfo, key1=keyColSets,
                         key2=keyColStats, all.y=TRUE, removeKey2ColNAs=TRUE)
    res <- setNames(df[[col]], df[[setsIDcol]])
    statsIDcol <- checkIfIDwasReplacedAfterMerging(statsIDcol, df)
    names(res)[is.na(names(res))] <- df[[statsIDcol]][is.na(names(res))]
    return(res)
}

processStats <- function(stats, col) {
    isRank <- endsWith(col, "_rank")
    stats  <- sort(stats, decreasing=!isRank)
    stats  <- stats[unique(names(stats))]
    stats  <- ifelse(!isRank, `+`, `-`)(stats)
    return(stats)
}

prepareStatsCol <- function(col, stats) {
    if (is.vector(stats)) {
        col <- "values"
    } else if (is.null(col)) {
        cols <- paste0(c("rankProduct", "spearman", "pearson", "GSEA"), "_rank")
        col  <- cols[cols %in% colnames(stats)][[1]]
        if (!col %in% colnames(stats)) {
            stop("no suitable column to analyse; explicitly set argument 'col'")
        } else {
            msg <- sprintf(
                paste("Ordering results by column '%s'; to manually select",
                      "column to order by, please set argument 'col'"), col)
            message(msg)
        }
    } else if (!col %in% colnames(stats)) {
        stop(sprintf("specified column '%s' not found", col))
    }
    return(col)
}

#' Analyse drug set enrichment
#'
#' @param stats Named numeric vector or either a \code{similarPerturbations} or
#'   a \code{targetingDrugs} object (obtained after running
#'   \code{\link{rankSimilarPerturbations}} or
#'   \code{\link{predictTargetingDrugs}}, respectively)
#' @param sets Named list of characters: named sets containing compound
#'   identifiers (obtain drug sets by running \code{prepareDrugSets()})
#' @param col Character: name of the column to use for statistics (only required
#'   if class of \code{stats} is either \code{similarPerturbations} or
#'   \code{targetingDrugs})
#' @inheritParams fgsea::fgsea
#' @inheritDotParams fgsea::fgsea -pathways -stats -nperm -maxSize
#' @param keyColSets Character: column from \code{sets} to compare with column
#'   \code{keyColStats} from \code{stats}; automatically selected if \code{NULL}
#' @param keyColStats Character: column from \code{stats} to compare with column
#'   \code{keyColSets} from \code{sets}; automatically selected if \code{NULL}
#'
#' @importFrom fgsea fgsea
#'
#' @family functions for drug set enrichment analysis
#' @return Enrichment analysis based on GSEA
#' @export
#'
#' @examples
#' descriptors <- loadDrugDescriptors()
#' drugSets <- prepareDrugSets(descriptors)
#'
#' # Analyse drug set enrichment in ranked targeting drugs for a differential
#' # expression profile
#' data("diffExprStat")
#' gdsc      <- loadExpressionDrugSensitivityAssociation("GDSC")
#' predicted <- predictTargetingDrugs(diffExprStat, gdsc)
#'
#' analyseDrugSetEnrichment(drugSets, predicted)
analyseDrugSetEnrichment <- function(sets, stats, col=NULL, nperm=10000,
                                     maxSize=500, ..., keyColSets=NULL,
                                     keyColStats=NULL) {
    message("Matching compounds with those available in drug sets...")
    col   <- prepareStatsCol(col, stats)
    stats <- matchStatsWithDrugSetsID(
        sets=sets, stats=stats, col=col,
        keyColSets=keyColSets, keyColStats=keyColStats)
    stats <- processStats(stats, col)

    message("Performing enrichment analysis...")
    gseaRes <- suppressWarnings(
        fgsea(sets, stats, nperm=nperm, maxSize=maxSize, ...))
    gseaRes <- gseaRes[order(gseaRes$padj), ]
    colnames(gseaRes)[[1]] <- "descriptor"
    return(gseaRes)
}

#' Plot drug set enrichment
#'
#' @inheritParams analyseDrugSetEnrichment
#' @param selectedSets Character: drug sets to plot (if \code{NULL}, plot all)
#'
#' @family functions for drug set enrichment analysis
#' @return List of GSEA plots per drug set
#' @export
#'
#' @importFrom pbapply pblapply
#'
#' @examples
#' descriptors <- loadDrugDescriptors()
#' drugSets <- prepareDrugSets(descriptors)
#'
#' # Analyse drug set enrichment in ranked targeting drugs for a differential
#' # expression profile
#' data("diffExprStat")
#' gdsc      <- loadExpressionDrugSensitivityAssociation("GDSC")
#' predicted <- predictTargetingDrugs(diffExprStat, gdsc)
#'
#' plotDrugSetEnrichment(drugSets, predicted)
plotDrugSetEnrichment <- function(sets, stats, col="rankProduct_rank",
                                  selectedSets=NULL, keyColSets=NULL,
                                  keyColStats=NULL) {
    message("Matching compounds with those available in drug sets...")
    col   <- prepareStatsCol(col, stats)
    stats <- matchStatsWithDrugSetsID(
        sets=sets, stats=stats, col=col,
        keyColSets=keyColSets, keyColStats=keyColStats)
    stats <- processStats(stats, col)

    message("Plotting enrichment analysis...")
    plotGSEAperSet <- function(k, sets, stats, col) {
        plot <- plotGSEA(stats, sets[[k]], title=names(sets)[[k]])
        return(plot)
    }
    if (!is.null(selectedSets)) sets <- sets[selectedSets]
    plots <- pblapply(seq(sets), plotGSEAperSet, sets, stats, col)
    names(plots) <- names(sets)
    return(plots)
}
