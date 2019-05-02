#' Load CMap data
#'
#' Load CMap data. If \code{file} does not exist, it will first be downloaded.
#'
#' @note If \code{type = compoundInfo}, two files from
#' \strong{The Drug Repurposing Hub} will be downloaded containing information
#' about drugs and perturbations. The files will be named \code{file} with
#' \code{_drugs} and \code{_samples} before their extension, respectively.
#'
#' @param file Character: path to file
#' @param type Character: type of data to load (\code{metadata},
#' \code{geneInfo}, \code{zscores} or \code{compoundInfo})
#' @param zscoresId Character: identifiers to partially load z-scores file
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
#' if (interactive()) {
#'   loadCmapData("cmapZscores.gctx.gz", "zscores",
#'                cmapMetadataKnockdown$sig_id)
#' }
loadCMapData <- function(file, type=c("metadata", "geneInfo", "zscores",
                                      "compoundInfo"),
                         zscoresId=NULL) {
    type <- match.arg(type)
    nas  <- c("NA", "na", "-666", "-666.0", "-666 -666", "-666 -666|-666 -666",
              "-666.000000", "-666.0|-666.000000")
    if (type == "metadata") {
        link <- paste0(
            "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&",
            "format=file&", "file=GSE92742_Broad_LINCS_sig_info.txt.gz")
        downloadIfNeeded(file, link)
        data <- fread(file, sep="\t", na.strings=nas)

        data$pert_dose[data$pert_dose == "300.0|300.000000"] <- 300
        data$pert_dose <- as.numeric(data$pert_dose)
        data$pert_idose[data$pert_idose == "300 ng|300 ng"] <- "300 ng"
    } else if (type == "geneInfo") {
        link <- paste0(
            "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&",
            "format=file&", "file=GSE92742_Broad_LINCS_gene_info.txt.gz")
        downloadIfNeeded(file, link)
        data <- fread(file, sep="\t", na.strings=nas)
    } else if (type == "zscores") {
        link <- paste0(
            "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&",
            "format=file&file=GSE92742_Broad_LINCS_Level5_COMPZ.",
            "MODZ_n473647x12328.gctx.gz")
        downloadIfNeeded(file, link)
        data <- new("GCT", src=file, rid=NULL, cid=zscoresId,
                    set_annot_rownames=FALSE, matrix_only=FALSE)@mat
    } else if (type == "compoundInfo") {
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
        downloadIfNeeded(file[["drugs"]], link)

        # Replace separation symbols for targets
        drugData <- readAfterComments(file[["drugs"]])
        drugData$target <- gsub("|", ", ", drugData$target, fixed=TRUE)

        # Process perturbation data
        link <- paste0(
            "https://s3.amazonaws.com/data.clue.io/repurposing/downloads/",
            "repurposing_samples_20180907.txt")
        downloadIfNeeded(file[["samples"]], link)
        pertData <- readAfterComments(file[["samples"]])
        pertData <- pertData[ , c("pert_iname", "expected_mass", "smiles",
                                  "InChIKey", "pubchem_cid")]
        pertData <- collapseDuplicatedRows(pertData, "pert_iname")
        data <- merge(drugData, pertData, all=TRUE)
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
        metadata$pert_idose, c("%", "nM", "µM", "µL", "ng", "ng/µL", "ng/mL"))
    timepoint <- sortNumericUnitChar(metadata$pert_itime)

    list("perturbationType"=pertTypes,
         "cellLine"=sort(unique(metadata$cell_id)),
         "dosage"=dose,
         "timepoint"=timepoint)
}

#' Correlate differential expression scores per cell line
#'
#' @inheritParams compareAgainstCMap
#' @param method Character: correlation method
#' @param filtered Character: perturbations to filter
#'
#' @importFrom pbapply pblapply
#' @importFrom data.table data.table
#' @importFrom stats cor.test p.adjust
#'
#' @return Data frame with correlations statistics, p-value and q-value
#' @keywords internal
correlatePerCellLine <- function(cellLine, diffExprGenes, perturbations,
                                 filtered, method, pAdjustMethod="BH") {
    # Suppress warnings to avoid "Cannot compute exact p-value with ties"
    corPert <- function(pert, filtered, perturbations, diffExprGenes, method) {
        thisPert <- perturbations[ , filtered[pert]]
        cor.test(thisPert, diffExprGenes, method=method)
    }
    cors <- suppressWarnings(pblapply(
        seq(filtered), corPert, filtered, perturbations, diffExprGenes, method))

    cor  <- sapply(cors, "[[", "estimate")
    pval <- sapply(cors, "[[", "p.value")
    qval <- p.adjust(pval, pAdjustMethod)
    names(cor) <- names(pval) <- names(qval) <- filtered

    res <- data.table(names(cor), cor, pval, qval)
    names(res) <- c("identifier", sprintf("%s_%s", method,
                                          c("coef", "pvalue", "qvalue")))
    return(res)
}

#' Perform gene set enrichment (GSA) per cell line
#'
#' @inheritParams compareAgainstCMap
#' @inheritParams fgsea::fgsea
#' @param filtered Character: perturbations to filter
#'
#' @importFrom fgsea fgsea
#' @importFrom data.table data.table
#' @importFrom pbapply pblapply
#' @importFrom dplyr bind_rows
#'
#' @return Data frame containing gene set enrichment analysis (GSEA) results per
#' cell line
#' @keywords internal
performGSAperCellLine <- function(cellLine, perturbations, filtered, pathways) {
    perturbations <- unclass(perturbations[ , filtered])

    performGSAwithPerturbationSignature <- function(k, perturbation, pathways) {
        signature <- perturbation[ , k]
        names(signature) <- rownames(perturbation)
        suppressWarnings(fgsea(pathways=pathways, stats=sort(signature),
                               minSize=15, maxSize=500, nperm=1))
    }
    gsa <- pblapply(
        seq(ncol(perturbations)), performGSAwithPerturbationSignature,
        perturbations, pathways)
    # gsa <- plyr::compact(gsa) # in case of NULL elements
    names(gsa) <- colnames(perturbations)[seq(ncol(perturbations))]

    gsaRes <- cbind(identifier=rep(names(gsa), each=2), bind_rows(gsa))

    # Based on CMap paper (page e8) - weighted connectivity score (WTCS)
    ES_list <- c()
    for (geneID in unique(gsaRes$identifier)) {
        geneID_GSA  <- gsaRes[gsaRes$identifier %in% geneID, ]
        topGenes    <- geneID_GSA$ES[grepl("top", geneID_GSA$pathway)]
        bottomGenes <- geneID_GSA$ES[grepl("bottom", geneID_GSA$pathway)]

        if (sign(topGenes) == sign(bottomGenes))
            ES_final <- 0
        else if (sign(topGenes) != sign(bottomGenes))
            ES_final <- (topGenes - bottomGenes) / 2
        ES_list <- c(ES_list, ES_final)
    }

    names(ES_list) <- unique(gsaRes[["identifier"]])
    results <- data.table(names(ES_list), ES_list)
    names(results)[1:2] <- c("identifier", "WTCS")
    return(results)
}

#' Compare with cell line (print progress)
#'
#' @param metadata Data frame: perturbation metadata
#'
#' @keywords internal
compareWithCellLineProgress <- function(cellLine, cellLines, FUN, method,
                                        perturbations, metadata, ...) {
    current <- match(cellLine, cellLines)
    total   <- length(cellLines)
    msg <- paste("Performing %s using %s's %s perturbations",
                 "(%s out of %s cell lines)...")
    methodStr <- switch(method,
                        "spearman"="Spearman's correlation",
                        "pearson" ="Pearson's correlation",
                        "gsea"    ="GSEA")

    # Filter perturbation based on currently selected cell line
    cellLinePerts <- metadata$sig_id[
        tolower(metadata$cell_id) == tolower(cellLine)]
    filtered <- colnames(perturbations)[
        colnames(perturbations) %in% cellLinePerts]
    pertNumber <- length(filtered)

    cat(sprintf(msg, methodStr, cellLine, pertNumber, current, total),
        fill=TRUE)
    FUN(cellLine, perturbations=perturbations, filtered=filtered, ...)
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
    names(pathways) <- paste0(c("top", "bottom"), geneSize)
    return(pathways)
}

#' Compare single method
#'
#' @inheritParams prepareGSEApathways
#' @param method Character: comparison method (\code{spearman}, \code{pearson}
#'   or \code{gsea}; multiple methods may be selected at once)
#' @param perturbations \code{cmapPerturbations} object: CMap perturbations
#'   (check \code{\link{loadCMapPerturbations}})
#' @param cellLineMean Boolean: add a column with the mean score across cell
#'   lines? If \code{cellLineMean = "auto"} (default) the mean score will be
#'   added if more than one cell line is available
#'
#' @importFrom utils head tail
#' @importFrom tidyr gather
#' @importFrom dplyr bind_rows
#'
#' @keywords internal
compareSingleMethod <- function(method, diffExprGenes=diffExprGenes,
                                perturbations=perturbations, geneSize=geneSize,
                                cellLineMean=cellLineMean) {
    startTime <- Sys.time()
    metadata  <- attr(perturbations, "metadata")
    cellLine  <- unique(metadata$cell_id)

    pathways <- NULL
    if (method %in% c("spearman", "pearson")) {
        cat("Subsetting perturbations based on intersecting genes...",
            fill=TRUE)
        genes <- intersect(names(diffExprGenes), rownames(perturbations))
        diffExprGenes        <- diffExprGenes[genes]
        class(perturbations) <- tail(class(perturbations), 1)
        perturbations        <- perturbations[genes, ]

        # Correlate per cell line
        cellLineRes <- lapply(
            cellLine, compareWithCellLineProgress, cellLines=cellLine,
            correlatePerCellLine, method, diffExprGenes=diffExprGenes,
            perturbations=perturbations, metadata=metadata, method)
        colnameSuffix <- sprintf("_%s_coef", method)
    } else if (method == "gsea") {
        pathways <- prepareGSEApathways(diffExprGenes, geneSize)
        cellLineRes <- lapply(
            cellLine, compareWithCellLineProgress, cellLines=cellLine,
            performGSAperCellLine, method, perturbations=perturbations,
            metadata=metadata, pathways=pathways)
        colnameSuffix <- "_WTCS"
    }
    names(cellLineRes) <- cellLine
    data <- bind_rows(cellLineRes)

    # Set whether to calculate the mean value across cell lines
    if (cellLineMean == "auto") cellLineMean <- length(cellLine) > 1

    if (cellLineMean) {
        scoreCol <- 2
        # Remove cell line information from the identifier
        allIDs <- gsub("\\_[A-Z].*\\_", "\\_", data$identifier)
        idsFromMultipleCellLines <- names(table(allIDs)[table(allIDs) > 1])
        names(idsFromMultipleCellLines) <- idsFromMultipleCellLines

        cellLine <- gsub(".*\\_([A-Z].*)\\_.*", "\\1", data$identifier)
        names(cellLine) <- data$identifier

        res <- pblapply(idsFromMultipleCellLines, function(id, allIDs, score,
                                                         cellLine) {
            list(cellLines=paste(cellLine[id == allIDs], collapse=", "),
                 mean=mean(score[id == allIDs]))
        }, allIDs=allIDs, score=data[[scoreCol]], cellLine=cellLine)

        avg <- sapply(res, "[[", 2)
        avgDF <- data.frame(names(avg), avg, stringsAsFactors=FALSE)
        colnames(avgDF) <- colnames(data)[c(1, scoreCol)]
        data <- bind_rows(list(data, avgDF))
        cellLines <- rbind(
            data.frame(cellLines=cellLine,
                       summarised=allIDs %in% idsFromMultipleCellLines),
            data.frame(cellLines=sapply(res, "[[", 1), summarised=TRUE))
        attr(data, "cellLines") <- cellLines
    }

    # Relabel the "identifier" column name to be more descriptive
    pertType <- unique(metadata$pert_type)
    if (length(pertType) == 1) {
        pertTypes <- getCMapPerturbationTypes()
        pertType  <- names(pertTypes[pertTypes == pertType])

        if (pertType == "Compound")
            id <- "compound_perturbation"
        else if (grepl("biological agents", pertType))
            id <- "biological_agent_perturbation"
        else
            id <- "gene_perturbation"
        colnames(data)[colnames(data) == "identifier"] <- id
    }
    rownames(data) <- data$genes

    # Add pathway information based on GSEA run if available
    if (method == "gsea" && !is.null(pathways)) {
        attr(data, "pathways") <- pathways
    }

    # Report run settings and time
    diffTime <- format(round(Sys.time() - startTime, 2))
    msg <- paste0("Comparison against %s perturbations using '%s' method %s",
                  "performed in %s")
    extra <- ifelse(method == "gsea",
                    sprintf("(gene size of %s) ", geneSize), "")
    message(sprintf(msg, ncol(perturbations), method, extra, diffTime))
    return(data)
}

#' Compare differential expression results against CMap perturbations
#'
#' Weighted connectivity scores (WTCS) are calculated when \code{method} is
#' \code{gsea}. For more information on WTCS, read
#' \url{https://clue.io/connectopedia/cmap_algorithms}.
#'
#' @details Order results according to the mean correlation coefficient (if
#' \code{method} is \code{spearman} or \code{pearson}) or the weighted
#' connectivity score (WTCS) score (if \code{method} is \code{gsea}) across cell
#' lines (if \code{cellLineMean} is \code{TRUE}; otherwise results are ordered
#' based on the first cell line alone).
#'
#' @inheritParams compareSingleMethod
#' @param rankPerturbationsByCellLine Boolean: when ranking, also rank
#'   perturbations regarding individual cell lines
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
#' compareAgainstCMap(diffExprStat, perturbations)
#'
#' # Compare using only Spearman's correlation
#' compareAgainstCMap(diffExprStat, perturbations, method="spearman")
compareAgainstCMap <- function(diffExprGenes, perturbations,
                               method=c("spearman", "pearson", "gsea"),
                               geneSize=150, cellLineMean="auto",
                               rankPerturbationsByCellLine=FALSE) {
    supported <- c("spearman", "pearson", "gsea")
    method <- unique(method)
    method <- method[method %in% supported]

    if (length(method) == 0) {
        stop(paste(
            "Method must contain one of the following supported comparison",
            "methods:", paste(supported, collapse=", ")))
    }

    names(method) <- method
    res <- lapply(method, compareSingleMethod,
                  diffExprGenes=diffExprGenes, perturbations=perturbations,
                  geneSize=geneSize, cellLineMean=cellLineMean)
    cellLineAttr <- attr(res[[1]], "cellLine")

    pathways <- NULL
    if (!is.null(res$gsea)) pathways <- attr(res$gsea, "pathways")
    merged <- Reduce(merge, res)

    # Rank perturbations
    rankPerturbations <- function(data) {
        dataDf <- data[ , -c(1), drop=FALSE]
        ranked <- apply(-dataDf, 2, rank, na.last="keep")
        colnames(ranked) <- paste0(colnames(dataDf), "_rank")
        mode(ranked) <- "integer"
        return(cbind(data, ranked))
    }
    ranked <- rankPerturbations(merged)

    # Inherit metadata from perturbations and other useful information
    attr(ranked, "metadata")      <- attr(perturbations, "metadata")
    attr(ranked, "geneInfo")      <- attr(perturbations, "geneInfo")
    attr(ranked, "compoundInfo")  <- attr(perturbations, "compoundInfo")
    attr(ranked, "diffExprGenes") <- diffExprGenes
    attr(ranked, "cellLine")      <- cellLineAttr
    attr(ranked, "pathways")      <- pathways

    class(ranked) <- c("cmapComparison", class(ranked))
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

    if (!is.null(cellLine))
        metadata <- metadata[tolower(metadata$cell_id) %in% tolower(cellLine), ]

    if (!is.null(timepoint))
        metadata <- metadata[metadata$pert_itime %in% timepoint, ]

    if (!is.null(dosage))
        metadata <- metadata[metadata$pert_idose %in% dosage, ]

    if (!is.null(perturbationType)) {
        tmp <- getCMapPerturbationTypes()[perturbationType]
        if (!is.na(tmp)) perturbationType <- tmp
        metadata <- metadata[metadata$pert_type %in% perturbationType, ]
    }

    return(metadata)
}

#' Load CMap perturbation data
#'
#' @param metadata Data frame (CMap metadata) or character (respective filepath)
#' @param zscores Data frame (GCTX z-scores) or character (respective filepath)
#' @param geneInfo Data frame (CMap gene info) or character (respective
#'   filepath)
#' @param compoundInfo Data frame (CMap compound info) or character (respective
#'   filepath)
#'
#' @importFrom R.utils gunzip
#' @importFrom methods new
#'
#' @return Perturbation data from CMap as a data table
#' @export
#' @examples
#' if (interactive()) {
#'   metadata <- loadCMapData("cmapMetadata.txt", "metadata")
#'   metadata <- filterCMapMetadata(metadata, cellLine="HepG2")
#'   zscores  <- loadCMapData("cmapZscores.gctx", "zscores", metadata$sig_id)
#'   geneInfo <- loadCMapData("cmapGeneInfo.txt", "geneInfo")
#'   loadCMapPerturbations(metadata, zscores, geneInfo)
#' }
loadCMapPerturbations <- function(metadata, zscores, geneInfo,
                                  compoundInfo=NULL) {
    if (is.character(metadata)) metadata <- loadCMapData(metadata, "metadata")
    if (is.character(geneInfo)) geneInfo <- loadCMapData(geneInfo, "geneInfo")
    if (is.character(zscores)) {
        zscores <- loadCMapData(zscores, "zscores", metadata$sig_id)
    }
    if (is.character(compoundInfo)) {
        compoundInfo <- loadCMapData(compoundInfo, "compoundInfo")
    }

    rownames(zscores) <- geneInfo$pr_gene_symbol[
        match(rownames(zscores), geneInfo$pr_gene_id)]
    attr(zscores, "metadata") <- metadata
    attr(zscores, "geneInfo") <- geneInfo
    attr(zscores, "compoundInfo") <- compoundInfo
    class(zscores) <- c("cmapPerturbations", class(zscores))
    return(zscores)
}

#' Get perturbation types
#'
#' @param controls Boolean: return perturbation types used as control?
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
        controlPerts <- c(
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
