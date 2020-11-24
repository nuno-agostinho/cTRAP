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

#' @importFrom binr bins bins.getvals
calculateEvenlyDistributedBins <- function(numbers, target.bins=15, minpts=NULL,
                                           ..., ids=NULL) {
    nas     <- is.na(numbers)
    numbers <- round(numbers[!nas])
    if (max(numbers) - min(numbers) == 0) return(setNames(numbers, ids))
    if (is.null(minpts)) minpts <- round( length(numbers) / target.bins / 5 )

    bin <- bins(numbers, target.bins=target.bins, minpts=minpts, ...)

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
#' @param table Data frame: drug descriptors
#' @param id Integer or character: index or name of the column containing
#'   identifiers
#' @param maxUniqueElems Numeric: maximum number of unique elements in a
#'   descriptor to consider when creating discrete drug sets
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
prepareDrugSets <- function(table, id=1, maxUniqueElems=15) {
    # Remove elements with no ID
    valid <- !is.na(table[[id]])
    table <- table[valid, ]

    isCharacter <- sapply(table, class) == "character"

    # Prepare sets from character columns if # unique values <= maxUniqueElems
    uniqueElems <- sapply(lapply(table[ , isCharacter, with=FALSE], unique),
                          length)
    sets <- names(uniqueElems[uniqueElems <= maxUniqueElems])
    subtable <- table[ , sets, with=FALSE]

    res <- lapply(subtable, function(x, ids) split(ids, x, drop=TRUE),
                  ids=table[[id]])

    nonCharacterTable <- table[ , !isCharacter, with=FALSE]
    res2 <- pblapply(nonCharacterTable, calculateEvenlyDistributedBins,
                     ids=table[[id]])
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

#' Match identifiers between data and drug sets
#'
#' @inheritParams analyseDrugSetEnrichment
#'
#' @importFrom pbapply pbapply
#' @importFrom data.table data.table
#' @keywords internal
#'
#' @return Statistic values from input data and corresponding identifiers as
#'   names (if no match is found, the original identifier from argument
#'   \code{stats} is used)
matchStatsWithDrugSetsID <- function(sets, stats, col=NULL) {
    statsSuffix <- ".stats"
    setsSuffix  <- ".sets"

    # Prepare stats' compound information
    statsCompoundInfo <- attr(stats, "compoundInfo")
    if (is.vector(stats)) {
        statsInfo  <- data.table("id"=names(stats), "values"=stats)
        info       <- statsInfo
        statsIDcol <- colnames(statsInfo)[[1]]
    } else if (is(stats, "similarPerturbations")) {
        statsIDcol <- colnames(stats)[[1]]
        # Merge both metadata and compound information
        metadata   <- attr(stats, "metadata")
        pertsID    <- unique(metadata[ , c("pert_id", "pert_iname")])
        info       <- merge(pertsID, statsCompoundInfo, by="pert_iname",
                            all.x=TRUE)
        allInfo    <- merge(metadata, statsCompoundInfo, by="pert_iname",
                            all.x=TRUE)
        statsInfo  <- merge(stats, allInfo, by.x=colnames(stats)[[1]],
                            by.y="sig_id", all.x=TRUE)
    } else if (!is.null(statsCompoundInfo)) {
        info      <- statsCompoundInfo
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
            "Argument 'stats' needs to be a vector or a 'similarPerturbations'",
            "object or have an attribute called 'compoundInfo'")
        stop(msg)
    }
    colnames(statsInfo) <- paste0(colnames(statsInfo), statsSuffix)
    colnames(info)      <- paste0(colnames(info), statsSuffix)
    statsIDcol          <- paste0(statsIDcol, statsSuffix)

    # Select column containing values of interest
    if (is.vector(stats)) {
        col <- "values"
    } else if (is.null(col)) {
        cols <- paste0(c("rankProduct", "spearman", "pearson", "GSEA"), "_rank")
        col  <- cols[paste0(cols, statsSuffix) %in% colnames(statsInfo)][[1]]
    }
    col <- paste0(col, statsSuffix)
    if (!col %in% colnames(statsInfo)) {
        stop(sprintf("specified column '%s' not found", col))
    }

    # Get drug sets' compound info
    setsCompoundInfo <- attr(sets, "compoundInfo")
    setsIDcol <- "id"
    if (is.null(setsCompoundInfo)) {
        setsCompoundInfo <- data.table("id"=unique(unlist(sets)))
    } else if (!setsIDcol %in% colnames(setsCompoundInfo)) {
        setsIDcol <- colnames(setsCompoundInfo)[[1]]
    }
    colnames(setsCompoundInfo) <- paste0(colnames(setsCompoundInfo), setsSuffix)
    setsIDcol <- paste0(setsIDcol, setsSuffix)

    # Check matches using available compound information
    intersectWithStatsID <- function(setInfo_k, statsInfo_j) {
        length(intersect(stripStr(setInfo_k), stripStr(statsInfo_j)))
    }
    intersected <- pbapply(info, 2, function(i)
        apply(setsCompoundInfo, 2, intersectWithStatsID, i))

    # Get column with most matches
    if (is.vector(intersected)) {
        setsCol  <- names(setsCompoundInfo)
        statsCol <- names(which.max(intersected))
    } else {
        setsCol  <- names(which.max(apply(intersected, 1, max)))
        statsCol <- names(which.max(apply(intersected, 2, max)))
    }

    # Ignore missing values in sets' compound information
    setsCompoundInfo <- setsCompoundInfo[!is.na(setsCompoundInfo[[setsCol]]), ]
    statsInfo[[statsCol]] <- as.character(statsInfo[[statsCol]])
    setsCompoundInfo[[setsCol]] <- as.character(setsCompoundInfo[[setsCol]])

    if (!is(statsInfo, "data.table")) {
        statsInfo <- data.table(statsInfo)
    }
    if (!is(setsCompoundInfo, "data.table")) {
        setsCompoundInfo <- data.table(setsCompoundInfo)
    }
    merged <- merge(statsInfo, setsCompoundInfo, by.x=statsCol, by.y=setsCol,
                    all.x=TRUE, suffixes=c(statsSuffix, setsSuffix))
    # If column was merged, take new column name
    if (setsIDcol == setsCol) setsIDcol <- statsCol

    # Create vector of statistical results whose names include matches
    matchedStats <- setNames(merged[[col]], merged[[setsIDcol]])
    names(matchedStats)[is.na(names(matchedStats))] <-
        merged[[statsIDcol]][is.na(names(matchedStats))]
    # Save original column used for values
    attr(matchedStats, "valuesCol") <- gsub(paste0(statsSuffix, "$"), "", col)
    return(matchedStats)
}

processStats <- function(stats, col) {
    isRank <- endsWith(col, "_rank")
    stats  <- sort(stats, decreasing=!isRank)
    stats  <- stats[unique(names(stats))]
    stats  <- (if (!isRank) `+` else `-`)(stats)
    return(stats)
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
                                     maxSize=500, ...) {
    message("Matching compounds with those available in drug sets...")
    stats <- matchStatsWithDrugSetsID(sets=sets, stats=stats, col=col)
    stats <- processStats(stats, attr(stats, "valuesCol"))

    message("Performing enrichment analysis...")
    gseaRes <- suppressWarnings(
        fgsea(sets, stats, nperm=nperm, maxSize=maxSize, ...))
    gseaRes <- gseaRes[order(gseaRes$padj), ]
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
                                  selectedSets=NULL) {
    message("Matching compounds with those available in drug sets...")
    stats <- matchStatsWithDrugSetsID(sets=sets, stats=stats, col=col)
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
