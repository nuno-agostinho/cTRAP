#' Load L1000 metadata
#'
#' The data will be downloaded if not available
#'
#' @param l1000metadataFile Character: filepath to L1000 metadata
#'
#' @importFrom data.table fread
#'
#' @return Metadata as a data table
#' @examples
#' loadL1000metadata("L1000metadata.txt")
loadL1000metadata <- function(l1000metadataFile) {
    # Download metadata (if not available)
    downloadIfNeeded(l1000metadataFile, paste0(
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&",
        "format=file&file=GSE92742%5FBroad%5FLINCS%5Fsig%5Finfo%2Etxt%2Egz"))

    cat("Loading L1000 metadata...", fill=TRUE)
    metadata <- fread(l1000metadataFile, sep="\t",
                      na.strings=c("NA", "na", "-666"))
    return(metadata)
}

#' List available conditions in L1000 datasets
#'
#' @inheritParams loadL1000metadata
#'
#' @export
#' @return List of conditions in L1000 datasets
getL1000Conditions <- function(l1000metadataFile) {
    info <- loadL1000metadata(l1000metadataFile)

    pertTypes <- getL1000PerturbationTypes()
    pertTypes <- names(pertTypes)[pertTypes %in% unique(info$pert_type)]

    list("Perturbation type"=pertTypes,
         "Cell line"=unique(info$cell_id),
         "Dosage"=unique(info$pert_idose),
         "Time points"=unique(info$pert_itime))
}

#' Calculate statistics
#'
#' @param perturbations
#' @param cellLine
#' @param gene
#' @param coef
#'
#' @return List of statistics
calculateStatistics <- function(perturbations, cellLine, gene, coef) {
    coef <- paste0(cellLine, coef)

    percentile_vector <- ecdf(unlist(perturbations[, coef, with=FALSE]))
    perc <- percentile_vector(perturbations[
        perturbations$genes %in% gene, coef, with=FALSE])

    res <- list(tStat_mean=mean(perturbations[[coef]], na.rm=TRUE),
                tStat_gene=as.numeric(
                    perturbations[perturbations$genes %in% gene, coef,
                                  with=FALSE]),
                tStat_percentile=perc[[cellLine]])
    return(res)
}

#' Correlate differential expression scores per cell line
#'
#' @param cellLine Character: cell line
#' @param diffExprGenes Numeric vector: statistic for differentially expressed
#' genes
#' @param perturbations
#' @param method Character: correlation method
#'
#' @importFrom pbapply pblapply
#' @importFrom data.table data.table
correlatePerCellLine <- function(cellLine, diffExprGenes, perturbations,
                                 method) {
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
    qval <- p.adjust(pval)
    names(cor) <- names(pval) <- names(qval) <- colnames(perturbation)

    res <- data.table(names(cor), cor, pval, qval)
    names(res) <- c("genes", sprintf("%s_t_%s_%s", cellLine, method,
                                     c("coef", "pvalue", "qvalue")))
    return(res)
}

#' Perform gene set enrichment (GSA) per cell line
#'
#' @importFrom fgsea fgsea
#' @importFrom data.table data.table
#' @importFrom pbapply pblapply
#' @importFrom plyr rbind.fill
#'
#' @return Data frame containing gene set enrichment analysis (GSEA) results per
#' cell line
performGSAperCellLine <- function(cellLine, countTable, perturbations, gsc) {
    perturbation <- perturbations[
        , tolower(attr(perturbations, "cellLines")) == tolower(cellLine)]

    # Run GSA with top 150 genes
    performGSAwithPerturbationSignature <- function(k, perturbation, gsc) {
        signature <- perturbation[ , k]
        names(signature) <- rownames(perturbation)
        fgsea(pathways=gsc$gsc, stats=sort(signature),
              minSize=15, maxSize=500, nperm=1)
    }

    cat("Performing GSA using perturbation signatures...", fill=TRUE)
    gsa <- pblapply(seq(ncol(perturbation)),
                    performGSAwithPerturbationSignature, perturbation, gsc)
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
#' where the name of the vector are gene names and the values are a statistic
#' that represents significance and magnitude of differentially expressed genes
#' (e.g. t-statistics)
#' @param geneSize Number: top and bottom differentially expressed genes to use
#' for gene set enrichment (GSE); if \code{method} is not \code{gsea}, this
#' argument does nothing
#'
#' @importFrom data.table setkeyv
#' @importFrom piano loadGSC
#'
#' @export
#' @examples
#' cellLines <- c("PC3", "VCAP", "A375")
#' timepoint <- c("96 h")
#' loadL1000perturbations("L1000", cellLines, timepoint)
#'
#' # Compare against L1000 using Spearman correlation
#' compareAgainstL1000(diffExprGenes, perturbations, cellLines,
#'                     method="spearman")
#'
#' # Compare against L1000 using Pearson correlation
#' compareAgainstL1000(diffExprGenes, perturbations, cellLines,
#'                     method="pearson")
#'
#' # Compare against L1000 using gene set enrichment analysis (GSEA)
#' compareAgainstL1000(diffExprGenes, perturbations, cellLines, method="gsea")
compareAgainstL1000 <- function(diffExprGenes, perturbations, cellLines,
                                method=c("spearman", "pearson", "gsea"),
                                geneSize=150) {
    method <- match.arg(method)
    if (method %in% c("spearman", "pearson")) {
        cellLineRes <- lapply(cellLines, correlatePerCellLine,
                              diffExprGenes, perturbations, method)
        colnameSuffix <- sprintf("_t_%s_coef", method)
    } else if (method == "gsea") {
        ordered     <- order(diffExprGenes)
        topGenes    <- names(diffExprGenes)[head(ordered, geneSize)]
        bottomGenes <- names(diffExprGenes)[tail(ordered, geneSize)]
        gsc <- loadGSC(matrix(c(
            c(topGenes, bottomGenes),
            c(rep(paste0("top", geneSize), geneSize),
              rep(paste0("bottom", geneSize), geneSize))), ncol=2))
        cellLineRes <- lapply(cellLines, performGSAperCellLine, countTable,
                              perturbations, gsc)
        colnameSuffix <- "_WTCS"
    }
    names(cellLineRes) <- cellLines

    # Merge results per cell line
    merged <- data.table(genes=unique(geneAnnot$`Associated Gene Name`))
    for (i in seq(cellLines))
        merged <- merge(merged, cellLineRes[[i]], all=TRUE, on="genes")

    data <- merged[rowSums(is.na(merged)) != ncol(merged) - 1, ]

    # Add mean calculated across cell lines
    data[ , paste0("Average", colnameSuffix)] <- rowMeans(
        data[ , grep(colnameSuffix, names(data)), with=FALSE], na.rm=TRUE)

    rownames(data) <- data$genes
    attr(data, "method") <- method
    class(data) <- c("L1000comparison", class(data))
    return(data)
}

calculate <- function(perturbations, cellLines, gene) {
    if (method %in% c("spearman", "pearson")) {
        colnameSuffix <- sprintf("_t_%s_coef", method)
    } else if (method == "gsea") {
        colnameSuffix <- "_WTCS"
    }

    cols <- c(cellLines, "Average")
    outputPerCellLine <- lapply(cols, calculateStatistics,
                                perturbations=perturbations, gene=gene,
                                coef=colnameSuffix)
    names(outputPerCellLine) <- cols
    return(cbind(Gene=gene, as.data.frame(outputPerCellLine)))
}

#' Load L1000 perturbation data
#'
#' @inheritParams loadL1000metadata
#' @param l1000zscoresFile Character: path to GCTX z-scores file
#' @param cellLine Character: cell line (if \code{NULL}, all values are loaded)
#' @param timepoint Character: timepoint (if \code{NULL}, all values are loaded)
#' @param dosage Character: dosage (if \code{NULL}, all values are loaded)
#' @param perturbationType Character: type of perturbation (if \code{NULL}, all
#' perturbation types are loaded)
#' @param sanitizeCompoundNames Boolean: replace identifiers with compound names
#'
#' @importFrom R.utils gunzip
#' @importFrom methods new
#'
#' @export
#' @return Perturbation data from L1000 as a data table
loadL1000perturbations <- function(l1000metadataFile, l1000zscoresFile,
                                   l1000geneFile, cellLine=NULL, timepoint=NULL,
                                   dosage=NULL, perturbationType=NULL,
                                   sanitizeCompoundNames=FALSE) {
    metadata <- loadL1000metadata(l1000metadataFile)

    # Filter elements of interest
    if (!is.null(cellLine))
        metadata <- metadata[tolower(metadata$cell_id) %in% tolower(cellLine), ]

    if (!is.null(timepoint))
        metadata <- metadata[metadata$pert_itime %in% timepoint, ]

    if (!is.null(dosage))
        metadata <- metadata[metadata$pert_idose %in% dosage, ]

    if (!is.null(perturbationType)) {
        tmp <- getL1000PerturbationTypes()[perturbationType]
        if (!is.na(tmp)) perturbationType <- tmp
        metadata <- metadata[metadata$pert_type %in% perturbationType, ]
    }

    # Get data from GCTX files based on the corresponding sig_ids
    sig_ids <- metadata$sig_id

    downloadIfNeeded(
        l1000zscoresFile,
        paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&",
               "format=file&file=GSE92742%5FBroad%5FLINCS%5FLevel5%5FCOMPZ%2",
               "EMODZ%5Fn473647x12328%2Egctx%2Egz"))
    selected_ds <- cmapR::parse.gctx(l1000zscoresFile, cid=sig_ids)@mat

    if (sanitizeCompoundNames) {
        # Change colnames per drug
        colnames(selected_ds) <- metadata$pert_iname[
            match(colnames(selected_ds), metadata$sig_id)]
    }

    # Change rownames
    downloadIfNeeded(
        l1000geneFile,
        paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92742&",
               "format=file&",
               "file=GSE92742%5FBroad%5FLINCS%5Fgene%5Finfo%2Etxt%2Egz"))
    geneInfo <- fread(l1000geneFile, sep="\t", na.strings=c("NA", "na", "-666"))

    rownames(selected_ds) <- geneInfo$pr_gene_symbol[
        match(rownames(selected_ds), geneInfo$pr_gene_id)]

    # Add attribute containing cell lines
    attr(selected_ds, "cellLines") <- metadata$cell_id
    return(selected_ds)
}

#' Get perturbation types
#'
#' @return Perturbation types and respective codes as used by L1000 datasets
getL1000PerturbationTypes <- function () {
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
