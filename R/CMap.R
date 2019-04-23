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
    nas  <- c("NA", "na", "-666", "-666.0", "-666 -666", "-666 -666|-666 -666")
    if (type == "metadata") {
        link <- paste0(
            "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&",
            "format=file&", "file=GSE92742_Broad_LINCS_sig_info.txt.gz")
        downloadIfNeeded(file, link)
        data <- fread(file, sep="\t", na.strings=nas)
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
getCMapConditions <- function(metadata, cellLine=NULL, timepoint=NULL,
                              dosage=NULL, perturbationType=NULL,
                              control=FALSE) {
    metadata <- filterCMapMetadata(metadata, cellLine=cellLine,
                                   timepoint=timepoint, dosage=dosage,
                                   perturbationType=perturbationType)
    pertTypes <- getCMapPerturbationTypes()
    pertTypes <- names(pertTypes)[pertTypes %in% unique(metadata$pert_type)]
    if (!control) {
        pertTypes <- grep("Control", pertTypes, value=TRUE, invert=TRUE,
                          fixed=TRUE)
    }

    list("Perturbation type"=pertTypes,
         "Cell line"=unique(metadata$cell_id),
         "Dosage"=unique(metadata$pert_idose),
         "Time points"=unique(metadata$pert_itime))
}

#' Correlate differential expression scores per cell line
#'
#' @inheritParams compareAgainstCMap
#' @param metadata Data frame: perturbation metadata
#' @param method Character: correlation method
#'
#' @importFrom pbapply pblapply
#' @importFrom data.table data.table
#' @importFrom stats cor.test p.adjust
#'
#' @return Data frame with correlations statistics, p-value and q-value
#' @keywords internal
correlatePerCellLine <- function(cellLine, diffExprGenes, perturbations,
                                 metadata, method, pAdjustMethod="BH") {
    # Filter perturbation based on current cell line
    filtered <- colnames(perturbations)[
        tolower(metadata$cell_id) == tolower(cellLine)]

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
    names(res) <- c("identifier", sprintf("%s_%s_%s", cellLine, method,
                                          c("coef", "pvalue", "qvalue")))
    return(res)
}

#' Perform gene set enrichment (GSA) per cell line
#'
#' @inheritParams compareAgainstCMap
#' @inheritParams fgsea::fgsea
#'
#' @importFrom fgsea fgsea
#' @importFrom data.table data.table
#' @importFrom pbapply pblapply
#' @importFrom plyr rbind.fill
#'
#' @return Data frame containing gene set enrichment analysis (GSEA) results per
#' cell line
#' @keywords internal
performGSAperCellLine <- function(cellLine, perturbations, pathways) {
    perturbation <- perturbations[
        , tolower(attr(perturbations, "metadata")$cell_id) == tolower(cellLine)]
    perturbation <- unclass(perturbation)

    performGSAwithPerturbationSignature <- function(k, perturbation, pathways) {
        signature <- perturbation[ , k]
        names(signature) <- rownames(perturbation)
        suppressWarnings(fgsea(pathways=pathways$gsc, stats=sort(signature),
                               minSize=15, maxSize=500, nperm=1))
    }
    gsa <- pblapply(seq(ncol(perturbation)),
                    performGSAwithPerturbationSignature, perturbation, pathways)
    # gsa <- plyr::compact(gsa) # in case of NULL elements
    names(gsa) <- colnames(perturbation)[seq(ncol(perturbation))]

    gsaRes <- cbind(identifier=rep(names(gsa), each=2), rbind.fill(gsa))

    # based on CMap paper (page e8) - weighted connectivity score (WTCS)
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
    names(results)[1:2] <- c("identifier", paste0(cellLine, "_WTCS"))
    return(results)
}

#' Compare against CMap datasets
#'
#' Weighted connectivity scores (WTCS) are calculated when \code{method} is
#' \code{gsea}. For more information on WTCS, read
#' \url{https://clue.io/connectopedia/cmap_algorithms}.
#'
#' @details Order results according to the correlation coefficient (if
#' \code{method} is \code{spearman} or \code{pearson}) or the weighted
#' connectivity score (WTCS) score (if \code{method} is \code{gsea}) of the mean
#' across cell lines (if \code{cellLineMean} is \code{TRUE}; otherwise results
#' are ordered based on the first cell line alone).
#'
#' @param diffExprGenes Numeric: named vector of differentially expressed genes
#'   where the name of the vector are gene names and the values are a statistic
#'   that represents significance and magnitude of differentially expressed
#'   genes (e.g. t-statistics)
#' @param perturbations \code{cmapPerturbations} object: file with CMap loaded
#'   perturbations (check \code{\link{loadCMapPerturbations}})
#' @param cellLine Character: cell line(s)
#' @param method Character: comparison method (\code{spearman}, \code{pearson}
#'   or \code{gsea})
#' @param geneSize Number: top and bottom differentially expressed genes to use
#'   for gene set enrichment (GSE) (only used if \code{method} is \code{gsea})
#' @param pAdjustMethod Character: method for p-value adjustment (for more
#'   details, see \code{\link{p.adjust.methods}}; only used if \code{method} is
#'   \code{spearman} or \code{pearson})
#' @param cellLineMean Boolean: add a column with the mean across cell lines? If
#' \code{"auto"} (default) the mean will be added if more than one cell line is
#' available
#'
#' @importFrom data.table setkeyv
#' @importFrom piano loadGSC
#' @importFrom utils head tail
#'
#' @return Data table with correlation or GSEA results comparing differential
#' gene expression values with those associated with CMap perturbations
#' @export
#'
#' @examples
#' cellLine <- "HepG2"
#' data("cmapPerturbationsSmallMolecules")
#' perturbations <- cmapPerturbationsSmallMolecules
#' data("diffExprStat")
#'
#' # Compare against CMap using Spearman correlation
#' compareAgainstCMap(diffExprStat, perturbations, cellLine,
#'                    method="spearman")
#'
#' # Compare against CMap using Pearson correlation
#' compareAgainstCMap(diffExprStat, perturbations, cellLine,
#'                    method="pearson")
#'
#' # Compare against CMap using gene set enrichment analysis (GSEA)
#' compareAgainstCMap(diffExprStat, perturbations, cellLine, method="gsea")
compareAgainstCMap <- function(diffExprGenes, perturbations,
                               method=c("spearman", "pearson", "gsea"),
                               geneSize=150, pAdjustMethod="BH",
                               cellLineMean="auto") {
    startTime <- Sys.time()
    method   <- match.arg(method)
    metadata <- attr(perturbations, "metadata")
    cellLine <- unique(metadata$cell_id)

    compareWithCellLineProgress <- function(cellLine, cellLines, FUN, method,
                                            ...) {
        current <- match(cellLine, cellLines)
        total   <- length(cellLines)
        msg <- paste("Performing %s using %s's perturbations",
                     "(%s out of %s cell lines)...")
        methodStr <- switch(method,
                            "spearman"="Spearman's correlation",
                            "pearson" ="Pearson's correlation",
                            "gsea"    ="GSEA")
        cat(sprintf(msg, methodStr, cellLine, current, total), fill=TRUE)
        FUN(cellLine, ...)
    }

    if (method %in% c("spearman", "pearson")) {
        cat(paste("Subsetting perturbations based on intersecting genes for",
                  "comparison..."), fill=TRUE)
        genes <- intersect(names(diffExprGenes), rownames(perturbations))
        diffExprGenes        <- diffExprGenes[genes]
        class(perturbations) <- tail(class(perturbations), 1)
        perturbations        <- perturbations[genes, ]

        # Correlate per cell line
        cellLineRes <- lapply(
            cellLine, compareWithCellLineProgress, cellLine,
            correlatePerCellLine, method, diffExprGenes, perturbations,
            metadata, method, pAdjustMethod)
        colnameSuffix <- sprintf("_%s_coef", method)
    } else if (method == "gsea") {
        ordered     <- order(diffExprGenes, decreasing=TRUE)
        topGenes    <- names(diffExprGenes)[head(ordered, geneSize)]
        bottomGenes <- names(diffExprGenes)[tail(ordered, geneSize)]
        gsc <- loadGSC(matrix(c(
            c(topGenes, bottomGenes),
            c(rep(paste0("top", geneSize), geneSize),
              rep(paste0("bottom", geneSize), geneSize))), ncol=2))
        cellLineRes <- lapply(cellLine, compareWithCellLineProgress, cellLine,
                              performGSAperCellLine, method, perturbations, gsc)
        colnameSuffix <- "_WTCS"
        pathways <- gsc$gsc
    }
    names(cellLineRes) <- cellLine

    # Merge results per cell line
    merged <- data.table("identifier"=unique(names(diffExprGenes)))
    for (i in seq(cellLine)) {
        # Remove cell line information from the identifier
        cellLineRes[[i]]$identifier <- gsub("\\_[A-Z].*\\_", "\\_",
                                            cellLineRes[[i]]$identifier)
        merged <- merge(merged, cellLineRes[[i]], all=TRUE, on="identifier")
    }

    data <- merged[rowSums(is.na(merged)) != ncol(merged) - 1, ]

    # Set whether to calculate the mean value across cell lines
    if (cellLineMean == "auto") cellLineMean <- length(cellLine) > 1

    if (cellLineMean) {
        # Calculate mean across cell lines and use it when ordering data
        orderCol <- paste0("Average", colnameSuffix)
        data[ , orderCol] <- rowMeans(
            data[ , grep(colnameSuffix, names(data)), with=FALSE], na.rm=TRUE)
    } else {
        # Use the results for the first cell line when ordering data
        orderCol <- grep(colnameSuffix, names(data))[[1]]
    }
    data <- data[order(data[[orderCol]], decreasing=TRUE)]

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
    attr(data, "method") <- method
    attr(data, "diffExprGenes") <- diffExprGenes
    attr(data, "perturbations") <- perturbations
    if (method == "gsea") attr(data, "pathways") <- pathways
    class(data) <- c("cmapComparison", class(data))

    # Report run settings and time
    diffTime <- format(round(Sys.time() - startTime, 2))
    msg <- paste0("Comparison against %s perturbations using '%s' method %s",
                  "performed in %s")
    extra <- ifelse(method == "gsea",
                    sprintf("(gene size of %s) ", geneSize), "")
    message(sprintf(msg, ncol(perturbations), method, extra, diffTime))
    return(data)
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
loadCMapPerturbations <- function(metadata, zscores, geneInfo, compoundInfo) {
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
#' @return Perturbation types and respective codes as used by CMap datasets
#' @export
#'
#' @examples
#' getCMapPerturbationTypes()
getCMapPerturbationTypes <- function () {
    c("Compound"="trt_cp",
      "Peptides and other biological agents (e.g. cytokine)"="trt_lig",
      "shRNA for loss of function (LoF) of gene"="trt_sh",
      "Consensus signature from shRNAs targeting the same gene"="trt_sh.cgs",
      "cDNA for overexpression of wild-type gene"="trt_oe",
      "cDNA for overexpression of mutated gene"="trt_oe.mut",
      "CRISPR for LLoF"="trt_xpr",
      "Controls - vehicle for compound treatment (e.g DMSO)"="ctl_vehicle",
      "Controls - vector for genetic perturbation (e.g empty vector, GFP)"="ctl_vector",
      "Controls - consensus signature from shRNAs that share a common seed sequence"="trt_sh.css",
      "Controls - consensus signature of vehicles"="ctl_vehicle.cns",
      "Controls - consensus signature of vectors"="ctl_vector.cns",
      "Controls - consensus signature of many untreated wells"="ctl_untrt.cns",
      "Controls - Untreated cells"="ctl_untrt")
}
