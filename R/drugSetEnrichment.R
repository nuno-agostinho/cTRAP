#' Prepare drug sets from a table with compound descriptors
#'
#' @param table Data frame: drug descriptors
#' @param id Character: column name containing compound identifiers
#' @param maxUniqueElems Numeric: maximum number of unique elements in a
#'   descriptor to consider when creating discrete drug sets
#'
#' @return Named list of characters: named drug sets with respective compound
#'   identifiers as list elements
#' @export
prepareDrugSets <- function(table, id, maxUniqueElems=15) {
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
#' @inheritDotParams fgsea::fgsea -pathways -stats -nperm
#'
#' @return Enrichment analysis based on GSEA
#' @export
analyseDrugSetEnrichment <- function(sets, stats, col="rankProduct_rank",
                                     nperm=10000, maxSize=500, ...) {
    if (is(stats, "similarPerturbations")) {
        # Use negative values to invert perturbation relevance
        rank   <- -stats[[col]]
        values <- na.omit(setNames(rank, stats[[1]])) # Discard missing values
        metadata     <- attr(stats, "metadata")
        compoundInfo <- attr(stats, "compoundInfo")

        # Convert to perturbation name (discard if non-matching)
        names(values) <- metadata[
            metadata$sig_id == names(values), ][["pert_iname"]]
        values <- values[!is.na(names(values))]

        # Convert to PubChem identifier
        names(values) <- compoundInfo[
            match(names(values), compoundInfo$pert_iname), ][["pubchem_cid"]]

        values <- values[!is.na(names(values))]
        dt <- data.table(PubChem=names(values), rank=values)[
            , .("rank"=max(rank)), by=PubChem]
        stats <- setNames(dt$rank, dt$PubChem)
    } else if (is(stats, "targetingDrugs")) {
        # Use negative values to invert perturbation relevance
        rank   <- -stats[[col]]
        values <- na.omit(setNames(rank, stats[[1]])) # Discard missing values
        compoundInfo <- attr(stats, "compoundInfo")
        stats <- setNames(stats[[col]], stats$compound)
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
