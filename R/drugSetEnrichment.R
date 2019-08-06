#' Load table with drug descriptors
#'
#' @importFrom data.table fread
#'
#' @return Data table with drug descriptors (first column, named )
#' @export
#' @examples
#' loadDrugDescriptors()
loadDrugDescriptors <- function() {
    source <- "NCI60"
    type   <- "2D"
    if (source == "NCI60" && type == "2D")
        link <- "2b9b8spodcwkbq6/compound_descriptors_NCI60_2D.txt"
    else
        stop("Unsupported source and type")

    link <- sprintf("https://www.dropbox.com/s/%s?raw=1", link)
    file <- sprintf("molecular_descriptors_%s_%s.txt", source, type)
    file <- downloadIfNotFound(link, file)
    table <- fread(file)
    return(table)
}

#' Prepare drug sets from a table with compound descriptors
#'
#' @param table Data frame: drug descriptors
#' @param id Integer or character: column index or name containing identifiers
#' @param maxUniqueElems Numeric: maximum number of unique elements in a
#'   descriptor to consider when creating discrete drug sets
#'
#' @family functions for drug set enrichment analysis
#' @return Named list of characters: named drug sets with respective compound
#'   identifiers as list elements
#' @export
#' @examples
#' descriptors <- loadDrugDescriptors()
#' prepareDrugSets(descriptors)
prepareDrugSets <- function(table, id=1, maxUniqueElems=15) {
    nas <- !is.na(table[[id]])
    table <- table[nas, ]
    uniqueElems <- sapply(lapply(table, unique), length)
    sets <- names(uniqueElems[uniqueElems <= maxUniqueElems])
    subtable <- table[ , sets, with=FALSE]

    res <- lapply(subtable, function(x, ids) {
        split(ids, x, drop=TRUE)
    }, ids=table[[id]])
    res <- unlist(res, recursive=FALSE)
    return(res)
}

#' Analyse drug set enrichment
#'
#' @param stats Named vector of characters or either a
#'   \code{similarPerturbations} or a \code{targetingDrugs} object (obtained
#'   after running \code{\link{rankSimilarPerturbations}} or
#'   \code{\link{predictTargetingDrugs}}, respectively)
#' @param sets Named list of characters: named sets containing compound
#'   identifiers (obtain drug sets by running \code{prepareDrugSets()})
#' @param col Character: name of the column to use for statistics (only required
#'   if class of \code{stats} is either \code{similarPerturbations} or
#'   \code{targetingDrugs})
#' @inheritParams fgsea::fgsea
#' @inheritDotParams fgsea::fgsea -pathways -stats -nperm -maxSize
#'
#' @family functions for drug set enrichment analysis
#' @return Enrichment analysis based on GSEA
#' @export
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
analyseDrugSetEnrichment <- function(sets, stats, col="rankProduct_rank",
                                     nperm=10000, maxSize=500, ...) {
    if (is(stats, "similarPerturbations")) {
        # Use negative values to invert perturbation relevance
        rank   <- stats[[col]]
        values <- na.omit(setNames(rank, stats[[1]])) # Discard missing values
        metadata     <- attr(stats, "metadata")
        compoundInfo <- attr(stats, "compoundInfo")

        # Convert to perturbation name (discard if non-matching)
        names(values) <- metadata[
            match(names(values), metadata$sig_id), ][["pert_iname"]]
        values <- values[!is.na(names(values))]

        # Convert to PubChem identifier
        names(values) <- compoundInfo[
            match(names(values), compoundInfo$pert_iname), ][["pubchem_cid"]]
        values <- values[!is.na(names(values))]

        dt <- data.table("PubChem"=names(values), "rank"=values)
        dt <- dt[ , list("rank"=max(rank)), by="PubChem"]
        stats <- setNames(-dt[[2]], dt[["PubChem"]])
    } else if (is(stats, "targetingDrugs")) {
        # Use negative values to invert perturbation relevance
        rank   <- stats[[col]]
        values <- na.omit(setNames(rank, stats[[1]])) # Discard missing values
        compoundInfo <- attr(stats, "compoundInfo")
        stats <- setNames(-stats[[col]], stats$compound)
    }
    stats   <- sort(stats)
    gseaRes <- suppressWarnings(
        fgsea(sets, stats, nperm=nperm, maxSize=maxSize, ...))
    gseaRes <- gseaRes[order(gseaRes$padj), ]
    return(gseaRes)
}

#' Plot drug set enrichment
#'
#' @inheritParams analyseDrugSetEnrichment
#'
#' @family functions for drug set enrichment analysis
#' @return List of GSEA plots per drug set
#' @export
plotDrugSetEnrichment <- function(sets, stats, col="rankProduct_rank") {
    plotGSEAperSet <- function(k) {
        ranks <- setNames(stats[[col]], stats$compound)
        plot  <- plotGSEA(ranks, sets[[k]])
        return(plot)
    }
    plots <- lapply(seq(sets), plotGSEAperSet)
    return(plots)
}
