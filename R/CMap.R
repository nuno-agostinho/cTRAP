# Retrieve CMap data -----------------------------------------------------------

#' Parse CMap identifier
#'
#' @param id Character: CMap identifier
#' @param cellLine Boolean: if \code{TRUE}, return cell line information from
#'   CMap identifier; else, return the CMap identifier without the cell line
#'
#' @family functions related with the ranking of CMap perturbations
#' @return Character vector with information from CMap identifiers
#' @export
#'
#' @examples
#' id <- c("CVD001_HEPG2_24H:BRD-K94818765-001-01-0:4.8",
#'         "CVD001_HEPG2_24H:BRD-K96188950-001-04-5:4.3967",
#'         "CVD001_HUH7_24H:BRD-A14014306-001-01-1:4.1")
#' parseCMapID(id, cellLine=TRUE)
#' parseCMapID(id, cellLine=FALSE)
parseCMapID <- function(id, cellLine=FALSE) {
    if (cellLine) {
        # Retrieve cell line
        res <- gsub(".*\\_([A-Z].*)\\_.*", "\\1", id)
        # Assign missing values to identifiers of summarised perturbation scores
        res <- ifelse(grepl(":", res), NA, res)
    } else {
        # Remove cell line identifier
        res <- gsub("\\_[A-Z].*\\_", "\\_", id)
    }
    names(res) <- id
    return(res)
}

#' Get CMap perturbation types
#'
#' @param control Boolean: return perturbation types used as control?
#'
#' @family functions related with the ranking of CMap perturbations
#' @return Perturbation types and respective codes as used by CMap datasets
#' @export
#'
#' @examples
#' getCMapPerturbationTypes()
getCMapPerturbationTypes <- function (control=FALSE) {
    perts <- c(
        "Compound"="trt_cp",
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

loadCMapMetadata <- function(file, nas) {
    link <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742",
                   "&format=file&file=GSE92742_Broad_LINCS_sig_info.txt.gz")
    downloadIfNotFound(link, file)
    message(sprintf("Loading CMap metadata from %s...", file))
    data <- fread(file, sep="\t", na.strings=nas)

    # Fix issues with specific metadata values
    data$pert_dose[data$pert_dose == "300.0|300.000000"] <- 300
    data$pert_dose <- as.numeric(data$pert_dose)
    data$pert_idose[data$pert_idose == "300 ng|300 ng"] <- "300 ng"
    return(data)
}

prepareCMapZscores <- function(file, zscoresID=NULL) {
    link <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742",
                   "&format=file&file=GSE92742_Broad_LINCS_Level5_COMPZ.",
                   "MODZ_n473647x12328.gctx.gz")
    downloadIfNotFound(link, file, ask=TRUE)
    data <- normalizePath(file)
    attr(data, "genes")         <- readGctxIds(data, dimension="row")
    attr(data, "perturbations") <- processIds(
        zscoresID, readGctxIds(data, dimension="col"), type="cid")$ids
    return(data)
}

#' Load matrix of CMap perturbation's differential expression z-scores
#' (optional)
#'
#' @param data \code{perturbationChanges} object
#' @param inheritAttrs Boolean: convert to \code{perturbationChanges} object and
#'   inherit attributes from \code{data}?
#' @param verbose Boolean: print additional details?
#'
#' @family functions related with the ranking of CMap perturbations
#' @return Matrix containing CMap perturbation z-scores (genes as rows,
#'   perturbations as columns)
#' @export
#'
#' @examples
#' metadata <- loadCMapData("cmapMetadata.txt", "metadata")
#' metadata <- filterCMapMetadata(metadata, cellLine="HepG2")
#' \dontrun{
#' perts <- prepareCMapPerturbations(metadata, "cmapZscores.gctx",
#'                                   "cmapGeneInfo.txt")
#' zscores <- loadCMapZscores(perts[ , 1:10])
#' }
loadCMapZscores <- function(data, inheritAttrs=FALSE, verbose=TRUE) {
    if (verbose) {
        msg <- paste("Loading CMap perturbation's differential expression",
                     "z-scores from %s...")
        message(sprintf(msg, data))
    }
    zscores  <- new("GCT", src=data, cid=colnames(data), verbose=verbose)@mat
    geneInfo <- attr(data, "geneInfo")
    if (!is.null(geneInfo)) {
        rownames(zscores) <- geneInfo$pr_gene_symbol[
            match(rownames(zscores), geneInfo$pr_gene_id)]
        if (!setequal(attr(data, "genes"), rownames(zscores)))
            zscores <- zscores[attr(data, "genes"), , drop=FALSE]
    }

    if (inheritAttrs) {
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
    downloadIfNotFound(link, file[["drugs"]])

    # Replace separation symbols for targets
    message(sprintf("Loading CMap compound data [1/2] from %s...",
                    file[["drugs"]]))
    drugData <- readAfterComments(file[["drugs"]])
    drugData$target <- gsub("|", ", ", drugData$target, fixed=TRUE)

    # Process perturbation data
    link <- paste0(
        "https://s3.amazonaws.com/data.clue.io/repurposing/downloads/",
        "repurposing_samples_20180907.txt")
    downloadIfNotFound(link, file[["samples"]])

    message(sprintf("Loading CMap compound data [2/2] from %s...",
                    file[["samples"]]))
    pertData <- readAfterComments(file[["samples"]])
    pertData <- pertData[ , c("pert_iname", "expected_mass", "smiles",
                              "InChIKey", "pubchem_cid")]
    pertData <- unique(pertData)
    pertData <- aggregate(pertData, by=list(pertData$pert_iname),
                          function(x) paste(unique(na.omit(x)), collapse=", "))
    data <- merge(drugData, pertData[ , -1], all=TRUE)
    data[data == ""] <- NA # Fix missing values
    return(data)
}

#' @include utils.R
loadCMapGeneInfo <- function(file, nas) {
    link <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742",
                   "&format=file&file=GSE92742_Broad_LINCS_gene_info.txt.gz")
    downloadIfNotFound(link, file)
    message(sprintf("Loading CMap gene information from %s...", file))
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
#' (for performance reasons; if \code{NULL}, all identifiers will be loaded)
#'
#' @importFrom data.table fread
#' @importFrom tools file_ext file_path_sans_ext
#'
#' @family functions related with the ranking of CMap perturbations
#' @return Metadata as a data table
#' @export
#'
#' @examples
#' # Load CMap metadata (data is automatically downloaded if not available)
#' cmapMetadata <- loadCMapData("cmapMetadata.txt", "metadata")
#'
#' # Load CMap gene info
#' loadCMapData("cmapGeneInfo.txt", "geneInfo")
#' \dontrun{
#' # Load CMap zscores based on filtered metadata
#' cmapMetadataKnockdown <- filterCMapMetadata(
#'   cmapMetadata, cellLine="HepG2",
#'   perturbationType="Consensus signature from shRNAs targeting the same gene")
#' loadCMapData("cmapZscores.gctx.gz", "zscores", cmapMetadataKnockdown$sig_id)
#' }
loadCMapData <- function(file, type=c("metadata", "geneInfo", "zscores",
                                      "compoundInfo"),
                         zscoresID=NULL) {
    if (is.null(file)) stop("'file' cannot be NULL, please provide a filename")

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
#' @family functions related with the ranking of CMap perturbations
#' @return List of conditions in CMap datasets
#' @export
#'
#' @examples
#' \dontrun{
#' cmapMetadata <- loadCMapData("cmapMetadata.txt", "metadata")
#' }
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

#' Filter CMap metadata
#'
#' @param metadata Data frame (CMap metadata) or character (respective filepath)
#' @param cellLine Character: cell line (if \code{NULL}, all values are loaded)
#' @param timepoint Character: timepoint (if \code{NULL}, all values are loaded)
#' @param dosage Character: dosage (if \code{NULL}, all values are loaded)
#' @param perturbationType Character: type of perturbation (if \code{NULL}, all
#' perturbation types are loaded)
#'
#' @family functions related with the ranking of CMap perturbations
#' @return Filtered CMap metadata
#' @export
#'
#' @examples
#' cmapMetadata <- loadCMapData("cmapMetadata.txt", "metadata")
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
        tmp <- getCMapPerturbationTypes(control=TRUE)[perturbationType]
        if (!all(is.na(tmp))) perturbationType <- tmp
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
#' @inheritDotParams filterCMapMetadata
#' @param loadZscores Boolean: load matrix of perturbation z-scores? Not
#'   recommended in systems with less than 30GB of RAM; if \code{FALSE},
#'   downstream functions will load and process the file directly chunk by
#'   chunk, resulting in a lower memory footprint
#'
#' @importFrom R.utils gunzip
#' @importFrom methods new
#'
#' @family functions related with the ranking of CMap perturbations
#' @return CMap perturbation data attributes and filename
#' @export
#'
#' @examples
#' metadata <- loadCMapData("cmapMetadata.txt", "metadata")
#' metadata <- filterCMapMetadata(metadata, cellLine="HepG2")
#' \dontrun{
#' prepareCMapPerturbations(metadata, "cmapZscores.gctx", "cmapGeneInfo.txt")
#' }
prepareCMapPerturbations <- function(metadata, zscores, geneInfo,
                                     compoundInfo=NULL, ...,
                                     loadZscores=FALSE) {
    if (is.character(metadata)) metadata <- loadCMapData(metadata, "metadata")
    if (!is.null(list(...))) metadata <- filterCMapMetadata(metadata, ...)

    if (is.character(geneInfo)) geneInfo <- loadCMapData(geneInfo, "geneInfo")
    if (is.character(zscores)) {
        zscores <- loadCMapData(zscores, "zscores", metadata$sig_id)
        attr(zscores, "zscoresFilename") <- as.character(zscores)
    }
    if (is.character(compoundInfo)) {
        compoundInfo <- loadCMapData(compoundInfo, "compoundInfo")
    }

    if (!is.null(zscores)) {
        attr(zscores, "genes") <- geneInfo$pr_gene_symbol[
            match(attr(zscores, "genes"), geneInfo$pr_gene_id)]
        attr(zscores, "metadata") <- metadata
        attr(zscores, "geneInfo") <- geneInfo
        attr(zscores, "compoundInfo") <- compoundInfo
        class(zscores) <- c("perturbationChanges", class(zscores))

        # Item information
        attr(zscores, "source") <- "CMap"
        attr(zscores, "type")   <- "perturbations"

        if (loadZscores) zscores <- loadCMapZscores(zscores, inheritAttrs=TRUE)
    }

    # Display summary message of loaded perturbations
    filters <- attr(metadata, "filter")
    summaryMsg <- sprintf(
        "\nSummary: %s CMap perturbations and %s genes",
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

# Compare against CMap perturbations -------------------------------------------

#' Calculate cell line mean
#'
#' @param data Data table: comparison against CMap data
#' @param cellLine Character: perturbation identifiers as names and respective
#' cell lines as values
#' @param metadata Data table: \code{data} metadata
#' @inheritParams rankSimilarPerturbations
#'
#' @importFrom dplyr bind_rows
#'
#' @return A list with two items:
#' \describe{
#' \item{\code{data}}{input \code{data} with extra rows containing cell line
#'   average scores (if calculated)}
#' \item{\code{rankingInfo}}{data table with ranking information}
#' \item{\code{metadata}}{metadata associated with output \code{data}, including
#'   for identifiers regarding mean cell line scores}
#' }
#' @keywords internal
calculateCellLineMean <- function(data, cellLine, metadata, rankPerCellLine) {
    scoreCol <- 2
    # Remove cell line information from the identifier
    allIDs <- parseCMapID(data[[1]], cellLine=FALSE)
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
        toRank <- rankPerCellLine | !isSummarised
        avgCellLines <- sapply(res, "[[", "cellLines")

        rankingInfo <- data.table(
            c(names(cellLine), names(avgCellLines)),
            c(toRank, rep(TRUE, length(avgCellLines))))

        # Append metadata associated with mean perturbation scores
        avgCellLinesMetadata <- metadata[
            match(names(avgCellLines), parseCMapID(metadata$sig_id)), ]
        avgCellLinesMetadata$sig_id  <- names(avgCellLines)
        avgCellLinesMetadata$cell_id <- avgCellLines
        avgCellLinesMetadata$distil_id <- NA
        metadataJoint <- rbind(avgCellLinesMetadata, metadata)
    } else {
        rankingInfo  <- data.table(names(cellLine), TRUE)
        dataJoint     <- data
        metadataJoint <- metadata
    }
    res <- list("reference"=dataJoint, "rankingInfo"=rankingInfo,
                "metadata"=metadataJoint)
    return(res)
}

#' Rank differential expression profile against CMap perturbations by similarity
#'
#' Compare differential expression results against CMap perturbations.
#'
#' @inherit rankAgainstReference
#' @param perturbations \code{perturbationChanges} object: CMap perturbations
#'   (check \code{\link{prepareCMapPerturbations}()})
#'
#' @aliases compareAgainstCMap
#' @family functions related with the ranking of CMap perturbations
#' @export
#'
#' @examples
#' # Example of a differential expression profile
#' data("diffExprStat")
#'
#' \dontrun{
#' # Download and load CMap perturbations to compare with
#' cellLine <- c("HepG2", "HUH7")
#' cmapMetadataCompounds <- filterCMapMetadata(
#'     "cmapMetadata.txt", cellLine=cellLine, timepoint="24 h",
#'     dosage="5 \u00B5M", perturbationType="Compound")
#'
#' cmapPerturbationsCompounds <- prepareCMapPerturbations(
#'     cmapMetadataCompounds, "cmapZscores.gctx", "cmapGeneInfo.txt",
#'     "cmapCompoundInfo_drugs.txt", loadZscores=TRUE)
#' }
#' perturbations <- cmapPerturbationsCompounds
#'
#' # Rank similar CMap perturbations (by default, Spearman's and Pearson's
#' # correlation are used, as well as GSEA with the top and bottom 150 genes of
#' # the differential expression profile used as reference)
#' rankSimilarPerturbations(diffExprStat, perturbations)
#'
#' # Rank similar CMap perturbations using only Spearman's correlation
#' rankSimilarPerturbations(diffExprStat, perturbations, method="spearman")
rankSimilarPerturbations <- function(input, perturbations,
                                     method=c("spearman", "pearson", "gsea"),
                                     geneSize=150, cellLineMean="auto",
                                     rankPerCellLine=FALSE, threads=1,
                                     chunkGiB=1, verbose=FALSE) {
    metadata  <- attr(perturbations, "metadata")
    cellLines <- length(unique(metadata$cell_id))
    rankedPerts <- rankAgainstReference(
        input, perturbations, method=method, geneSize=geneSize,
        cellLines=cellLines, cellLineMean=cellLineMean, rankByAscending=TRUE,
        rankPerCellLine=rankPerCellLine, threads=threads, chunkGiB=chunkGiB,
        verbose=verbose)

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
        colnames(rankedPerts)[[1]] <- id
    }
    class(rankedPerts) <- c("similarPerturbations", class(rankedPerts))
    return(rankedPerts)
}

# perturbationChanges object ---------------------------------------------------

#' Operations on a \code{perturbationChanges} object
#'
#' @param x \code{perturbationChanges} object
#' @param ... Extra arguments
#' @param perturbation Character (perturbation identifier) or a
#'   \code{similarPerturbations} table (from which the respective perturbation
#'   identifiers are retrieved)
#' @inheritParams compareWithAllMethods
#' @inheritParams plot.referenceComparison
#' @param title Character: plot title (if \code{NULL}, the default title depends
#'   on the context; ignored when plotting multiple perturbations)
#'
#' @importFrom methods is
#' @importFrom stats setNames
#'
#' @family functions related with the ranking of CMap perturbations
#' @return Subset, plot or return dimensions or names of a
#'   \code{perturbationChanges} object
#' @export
#'
#' @examples
#' data("diffExprStat")
#' data("cmapPerturbationsKD")
#'
#' compareKD <- rankSimilarPerturbations(diffExprStat, cmapPerturbationsKD)
#' EIF4G1knockdown <- grep("EIF4G1", compareKD[[1]], value=TRUE)
#' plot(cmapPerturbationsKD, EIF4G1knockdown, diffExprStat, method="spearman")
#' plot(cmapPerturbationsKD, EIF4G1knockdown, diffExprStat, method="pearson")
#' plot(cmapPerturbationsKD, EIF4G1knockdown, diffExprStat, method="gsea")
#'
#' data("cmapPerturbationsCompounds")
#' pert <- "CVD001_HEPG2_24H:BRD-A14014306-001-01-1:4.1"
#' plot(cmapPerturbationsCompounds, pert, diffExprStat, method="spearman")
#' plot(cmapPerturbationsCompounds, pert, diffExprStat, method="pearson")
#' plot(cmapPerturbationsCompounds, pert, diffExprStat, method="gsea")
#'
#' # Multiple cell line perturbations
#' pert <- "CVD001_24H:BRD-A14014306-001-01-1:4.1"
#' plot(cmapPerturbationsCompounds, pert, diffExprStat, method="spearman")
#' plot(cmapPerturbationsCompounds, pert, diffExprStat, method="pearson")
#' plot(cmapPerturbationsCompounds, pert, diffExprStat, method="gsea")
plot.perturbationChanges <- function(x, perturbation, input,
                                     method=c("spearman", "pearson", "gsea"),
                                     geneSize=150,
                                     genes=c("both", "top", "bottom"), ...,
                                     title=NULL) {
    plotPerturbationChanges(x=x, perturbation=perturbation, input=input,
                            method=method, geneset=NULL, geneSize=geneSize,
                            genes=genes, ..., title=title)
}

plotPerturbationChanges <- function(x, perturbation, input,
                                    method=c("spearman", "pearson", "gsea"),
                                    geneset=NULL, geneSize=150,
                                    genes=c("both", "top", "bottom"), ...,
                                    title=NULL) {
    method <- match.arg(method)
    if (is(perturbation, "similarPerturbations")) {
        perturbation <- perturbation[[1]]
    }

    cellLinePerts <- colnames(x)[
        parseCMapID(colnames(x), cellLine=FALSE) %in% perturbation]
    isSummaryPert <- length(cellLinePerts) > 0

    if (length(perturbation) == 0) {
        stop("a perturbation ID must be provided")
    } else if (length(perturbation) > 1) {
        stop("only one perturbation ID is currently supported")
    } else if (!perturbation %in% colnames(x) && !isSummaryPert) {
        stop("perturbation not found in the columns of 'x'")
    }

    if (!isSummaryPert) cellLinePerts <- perturbation
    names(cellLinePerts) <- cellLinePerts
    if (is.character(x)) {
        zscores <- loadCMapZscores(x[ , cellLinePerts], verbose=FALSE)
    } else {
        zscores <- unclass(x)
    }

    data <- lapply(cellLinePerts, function(pert, zscores) {
        sub <- zscores[ , pert, drop=FALSE]
        setNames(as.numeric(sub), rownames(sub))
    }, zscores)

    if (method != "gsea") {
        plot <- plotSingleCorr(data, perturbation, input, title=title)
    } else {
        if (is.null(geneset)) geneset <- prepareGSEAgenesets(input, geneSize)
        plot <- NULL
        areMultiplePerturbations <- length(seq(data)) > 1
        for (i in seq(data)) {
            dataset <- unclass(data[[i]])
            if (is.null(title) || areMultiplePerturbations) {
                title <- names(data)[[i]]
            }
            p <- plotGSEA(dataset, geneset, genes, title=title, ...,
                          compact=areMultiplePerturbations)
            plot <- c(plot, setNames(list(p), title))
        }
        names(plot) <- names(cellLinePerts)
        if (areMultiplePerturbations) {
            plot <- plot_grid(plotlist=plot, ncol=1)
        } else {
            plot <- plot[[1]]
        }
    }
    return(plot)
}

#' @rdname plot.perturbationChanges
#' @param i,j Character or numeric indexes specifying elements to extract
#' @param drop Boolean: coerce result to the lowest possible dimension?
#' @export
`[.perturbationChanges` <- function(x, i, j, drop=FALSE, ...) {
    if (is.character(x)) {
        out <- subsetData(x, i, j, "genes", "perturbations", nargs(), ...)
    } else {
        out <- NextMethod("[", drop=drop)
    }

    # Trim metadata to only contain subset information
    attrs <- attributes(x)
    if (!is.null(ncol(out)) && ncol(x) != ncol(out) &&
        !is.null(attrs$metadata)) {
        samples <- attrs$metadata$sig_id %in% colnames(out)
        attr(out, "metadata") <- NULL
        attrs$metadata <- attrs$metadata[samples, , drop=FALSE]
    }
    # Inherit input's attributes
    attrs <- attrs[!names(attrs) %in% names(attributes(out))]
    attributes(out) <- c(attributes(out), attrs)
    return(out)
}

#' @rdname plot.perturbationChanges
#' @export
dim.perturbationChanges <- function(x) {
    if (is.character(x)) {
        res <- vapply(dimnames(x), length, numeric(1))
    } else {
        res <- NextMethod("dim")
    }
    return(res)
}

#' @rdname plot.perturbationChanges
#' @export
dimnames.perturbationChanges <- function(x) {
    if (is.character(x)) {
        res <- list(attr(x, "genes"), attr(x, "perturbations"))
    } else {
        res <- NextMethod("dimnames")
    }
    return(res)
}

# similarPerturbations object --------------------------------------------------

#' Print a \code{similarPerturbations} object
#'
#' @param x \code{similarPerturbations} object
#' @param perturbation Character (perturbation identifier) or numeric
#'   (perturbation index)
#' @param ... Extra parameters passed to \code{print}
#'
#' @family functions related with the ranking of CMap perturbations
#' @return Information on \code{perturbationChanges} object or on specific
#'   perturbations
#' @export
print.similarPerturbations <- function(x, perturbation=NULL, ...) {
    if (is.null(perturbation)) {
        NextMethod("print")
    } else {
        if (is.numeric(perturbation)) perturbation <- x[[1]][perturbation]

        metadata <- attr(x, "metadata")
        if (!is.null(metadata)) {
            selectMetadata <- metadata[metadata$sig_id %in% perturbation]
            if (nrow(selectMetadata) == 0) {
                # Check to see if using identifiers referring to summary stats
                summaryID <- parseCMapID(metadata$sig_id, cellLine=FALSE)
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
