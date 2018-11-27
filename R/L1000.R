#' Download L1000 data
#'
#' The data will be downloaded if not available
#'
#' @param file Character: filepath
#' @param type Character: type of data to load
#' @param zscoresId Character: identifiers to partially load z-scores file
#' (for performance reasons)
#'
#' @importFrom data.table fread
#'
#' @return Metadata as a data table
#' @export
#'
#' @examples
#' # Download L1000 metadata
#' l1000metadata <- downloadL1000data("l1000metadata.txt", "metadata")
#'
#' # Download L1000 gene info
#' downloadL1000data("l1000geneInfo.txt", "geneInfo")
#'
#' # Download L1000 zscores based on filtered metadata
#' l1000metadataKnockdown <- filterL1000metadata(
#'   l1000metadata, cellLine="HepG2",
#'   perturbationType="Consensus signature from shRNAs targeting the same gene")
#'
#' if (interactive()) {
#'     downloadL1000data("l1000zscores.gctx.gz", "zscores",
#'                       l1000metadataKnockdown$sig_id)
#' }
downloadL1000data <- function(file, type=c("metadata", "geneInfo", "zscores"),
                              zscoresId=NULL) {
    type <- match.arg(type)
    nas  <- c("NA", "na", "-666", "-666.0", "-666 -666", "-666 -666|-666 -666")
    if (type == "metadata") {
        downloadIfNeeded(file, paste0(
            "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&",
            "format=file&",
            "file=GSE92742_Broad_LINCS_sig_info.txt.gz"))

        message("Loading L1000 metadata...")
        data <- fread(file, sep="\t", na.strings=nas)
    } else if (type == "geneInfo") {
        downloadIfNeeded(
            file,
            paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&",
                   "format=file&",
                   "file=GSE92742_Broad_LINCS_gene_info.txt.gz"))
        data <- fread(file, sep="\t", na.strings=nas)
    } else if (type == "zscores") {
        downloadIfNeeded(
            file, paste0(
                "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&",
                "format=file&file=GSE92742_Broad_LINCS_Level5_COMPZ.",
                "MODZ_n473647x12328.gctx.gz"))
        data <- new("GCT", src=file, rid=NULL, cid=zscoresId,
                    set_annot_rownames=FALSE, matrix_only=FALSE)@mat
    }
    return(data)
}

#' List available conditions in L1000 datasets
#'
#' Downloads metadata if not available
#'
#' @param metadata Data table: L1000 metadata
#' @param control Boolean: show controls for perturbation types?
#'
#' @return List of conditions in L1000 datasets
#' @export
#'
#' @examples
#' data("l1000metadata")
#' # l1000metadata <- downloadL1000data("l1000metadata.txt", "metadata")
#' getL1000conditions(l1000metadata)
getL1000conditions <- function(metadata, control=FALSE) {
    pertTypes <- getL1000perturbationTypes()
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
#' @inheritParams compareAgainstL1000
#' @param method Character: correlation method
#'
#' @importFrom pbapply pblapply
#' @importFrom data.table data.table
#' @importFrom stats cor.test p.adjust
#'
#' @return Data frame with correlations statistics, p-value and q-value
#' @keywords internal
correlatePerCellLine <- function(cellLine, diffExprGenes, perturbations,
                                 method, pAdjustMethod="BH") {
    cat(paste("Comparing with cell line", cellLine), fill=TRUE)
    perturbation <- perturbations[
        , tolower(attr(perturbations, "cellLines")) == tolower(cellLine)]

    # Select intersecting genes
    genes <- intersect(names(diffExprGenes), rownames(perturbation))
    diffExprGenes <- diffExprGenes[genes]

    # setkeyv(perturbation, "Gene")
    # genesKey <- match(rownames(perturbation), genes)
    # genesKey <- genesKey[!is.na(genesKey)]
    ref <- perturbation[genes, ]

    # Suppress warnings to avoid "Cannot compute exact p-value with ties"
    cors <- suppressWarnings(
        pblapply(seq(ncol(ref)),
                 function(k) cor.test(ref[ , k], diffExprGenes, method=method)))

    cor <- sapply(cors, "[[", "estimate")
    pval <- sapply(cors, "[[", "p.value")
    qval <- p.adjust(pval, pAdjustMethod)
    names(cor) <- names(pval) <- names(qval) <- colnames(perturbation)

    res <- data.table(names(cor), cor, pval, qval)
    names(res) <- c("genes", sprintf("%s_%s_%s", cellLine, method,
                                     c("coef", "pvalue", "qvalue")))
    attr(res, "perturbation") <- ref
    return(res)
}

#' Perform gene set enrichment (GSA) per cell line
#'
#' @inheritParams compareAgainstL1000
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
        , tolower(attr(perturbations, "cellLines")) == tolower(cellLine)]

    # Run GSA with top 150 genes
    performGSAwithPerturbationSignature <- function(k, perturbation, pathways) {
        signature <- perturbation[ , k]
        names(signature) <- rownames(perturbation)
        suppressWarnings(fgsea(pathways=pathways$gsc, stats=sort(signature),
                               minSize=15, maxSize=500, nperm=1))
    }

    cat("Performing GSA using perturbation signatures...", fill=TRUE)
    gsa <- pblapply(seq(ncol(perturbation)),
                    performGSAwithPerturbationSignature, perturbation, pathways)
    # gsa <- plyr::compact(gsa) # in case of NULL elements
    names(gsa) <- colnames(perturbation)[seq(ncol(perturbation))]

    gsaRes <- cbind(genes=rep(names(gsa), each=2), rbind.fill(gsa))

    # based on L1000 paper (page e8) - weighted connectivity score (WTCS)
    ES_list <- c()
    for (geneID in unique(gsaRes$genes)) {
        geneID_GSA  <- gsaRes[gsaRes$genes %in% geneID, ]
        topGenes    <- geneID_GSA$ES[grepl("top", geneID_GSA$pathway)]
        bottomGenes <- geneID_GSA$ES[grepl("bottom", geneID_GSA$pathway)]

        if (sign(topGenes) == sign(bottomGenes))
            ES_final <- 0
        else if (sign(topGenes) != sign(bottomGenes))
            ES_final <- (topGenes - bottomGenes) / 2
        ES_list <- c(ES_list, ES_final)
    }

    names(ES_list) <- unique(gsaRes[["genes"]])
    results <- data.table(names(ES_list), ES_list)
    names(results)[1:2] <- c("genes", paste0(cellLine, "_WTCS"))
    return(results)
}

#' Compare against L1000 datasets
#'
#' @param diffExprGenes Numeric: named vector of differentially expressed genes
#'   where the name of the vector are gene names and the values are a statistic
#'   that represents significance and magnitude of differentially expressed
#'   genes (e.g. t-statistics)
#' @param perturbations \code{l1000perturbations} object: file with L1000 loaded
#'   perturbations (check \code{\link{loadL1000perturbations}})
#' @param cellLine Character: cell line(s)
#' @param method Character: comparison method (\code{spearman}, \code{pearson}
#'   or \code{gsea})
#' @param geneSize Number: top and bottom differentially expressed genes to use
#'   for gene set enrichment (GSE) (only used if \code{method} is \code{gsea})
#' @param pAdjustMethod Character: method for p-value adjustment (for more
#'   details, see \code{\link{p.adjust.methods}}; only used if \code{method} is
#'   \code{spearman} or \code{pearson})
#'
#' @importFrom data.table setkeyv
#' @importFrom piano loadGSC
#' @importFrom utils head tail
#'
#' @return Data table with correlation or GSEA results comparing differential
#' gene expression values with those associated with L1000 perturbations
#' @export
#'
#' @examples
#' cellLine <- "HepG2"
#' data("l1000perturbationsSmallMolecules")
#' perturbations <- l1000perturbationsSmallMolecules
#' data("diffExprStat")
#'
#' # Compare against L1000 using Spearman correlation
#' compareAgainstL1000(diffExprStat, perturbations, cellLine,
#'                     method="spearman")
#'
#' # Compare against L1000 using Pearson correlation
#' compareAgainstL1000(diffExprStat, perturbations, cellLine,
#'                     method="pearson")
#'
#' # Compare against L1000 using gene set enrichment analysis (GSEA)
#' compareAgainstL1000(diffExprStat, perturbations, cellLine, method="gsea")
compareAgainstL1000 <- function(diffExprGenes, perturbations, cellLine,
                                method=c("spearman", "pearson", "gsea"),
                                geneSize=150, pAdjustMethod="BH") {
    method <- match.arg(method)
    if (method %in% c("spearman", "pearson")) {
        cellLineRes <- lapply(cellLine, correlatePerCellLine, diffExprGenes,
                              perturbations, method, pAdjustMethod)
        colnameSuffix <- sprintf("_%s_coef", method)
    } else if (method == "gsea") {
        ordered     <- order(diffExprGenes, decreasing=TRUE)
        topGenes    <- names(diffExprGenes)[head(ordered, geneSize)]
        bottomGenes <- names(diffExprGenes)[tail(ordered, geneSize)]
        gsc <- loadGSC(matrix(c(
            c(topGenes, bottomGenes),
            c(rep(paste0("top", geneSize), geneSize),
              rep(paste0("bottom", geneSize), geneSize))), ncol=2))
        cellLineRes <- lapply(cellLine, performGSAperCellLine, perturbations,
                              gsc)
        colnameSuffix <- "_WTCS"
        pathways <- gsc$gsc
    }
    names(cellLineRes) <- cellLine

    # Merge results per cell line
    merged <- data.table(genes=unique(names(diffExprGenes)))
    for (i in seq(cellLine))
        merged <- merge(merged, cellLineRes[[i]], all=TRUE, on="genes")

    data <- merged[rowSums(is.na(merged)) != ncol(merged) - 1, ]

    # Add mean calculated across cell lines
    data[ , paste0("Average", colnameSuffix)] <- rowMeans(
        data[ , grep(colnameSuffix, names(data)), with=FALSE], na.rm=TRUE)

    rownames(data) <- data$genes
    attr(data, "method") <- method
    attr(data, "diffExprGenes") <- diffExprGenes
    attr(data, "perturbations") <- perturbations
    if (method == "gsea") attr(data, "pathways") <- pathways
    class(data) <- c("l1000comparison", class(data))
    return(data)
}

#' Filter L1000 metadata
#'
#' @param metadata Data frame: metadata
#' @param cellLine Character: cell line (if \code{NULL}, all values are loaded)
#' @param timepoint Character: timepoint (if \code{NULL}, all values are loaded)
#' @param dosage Character: dosage (if \code{NULL}, all values are loaded)
#' @param perturbationType Character: type of perturbation (if \code{NULL}, all
#' perturbation types are loaded)
#'
#' @return Filtered L1000 metadata
#'
#' @export
#' @examples
#' data("l1000metadata")
#' # l1000metadata <- downloadL1000data("l1000metadata.txt", "metadata")
#' filterL1000metadata(l1000metadata, cellLine="HEPG2", timepoint="2 h",
#'                     dosage="25 ng/mL")
filterL1000metadata <- function(metadata, cellLine=NULL, timepoint=NULL,
                                dosage=NULL, perturbationType=NULL) {
    if (!is.null(cellLine))
        metadata <- metadata[tolower(metadata$cell_id) %in% tolower(cellLine), ]

    if (!is.null(timepoint))
        metadata <- metadata[metadata$pert_itime %in% timepoint, ]

    if (!is.null(dosage))
        metadata <- metadata[metadata$pert_idose %in% dosage, ]

    if (!is.null(perturbationType)) {
        tmp <- getL1000perturbationTypes()[perturbationType]
        if (!is.na(tmp)) perturbationType <- tmp
        metadata <- metadata[metadata$pert_type %in% perturbationType, ]
    }

    return(metadata)
}

#' Load L1000 perturbation data
#'
#' @param metadata Data frame: L1000 Metadata
#' @param zscores Data frame: GCTX z-scores
#' @param geneInfo Data frame: L1000 gene info
#' @param sanitizeCompoundNames Boolean: replace identifiers with compound names
#'
#' @importFrom R.utils gunzip
#' @importFrom methods new
#'
#' @return Perturbation data from L1000 as a data table
#' @export
#' @examples
#' if (interactive()) {
#'   metadata <- downloadL1000data("l1000metadata.txt", "metadata")
#'   metadata <- filterL1000metadata(metadata, cellLine="HepG2")
#'   zscores  <- downloadL1000data("l1000zscores.gctx", "zscores",
#'       metadata$sig_id)
#'   geneInfo <- downloadL1000data("l1000geneInfo.txt", "geneInfo")
#'   loadL1000perturbations(metadata, zscores, geneInfo)
#' }
loadL1000perturbations <- function(metadata, zscores, geneInfo,
                                   sanitizeCompoundNames=FALSE) {
    if (sanitizeCompoundNames) {
        # Change colnames per drug
        colnames(zscores) <- metadata$pert_iname[
            match(colnames(zscores), metadata$sig_id)]
    }

    # Change rownames
    rownames(zscores) <- geneInfo$pr_gene_symbol[
        match(rownames(zscores), geneInfo$pr_gene_id)]

    # Add attribute containing cell lines
    attr(zscores, "cellLines") <- metadata$cell_id
    class(zscores) <- c("l1000perturbations", class(zscores))
    return(zscores)
}

#' Get perturbation types
#'
#' @return Perturbation types and respective codes as used by L1000 datasets
#' @export
#'
#' @examples
#' getL1000perturbationTypes()
getL1000perturbationTypes <- function () {
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
