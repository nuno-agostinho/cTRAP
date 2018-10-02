#' Load ENCODE metadata
#'
#' @importFrom data.table fread
loadENCODEmetadata <- function(file, cellLine=NULL, lab=NULL, outputType=NULL) {
    data <- fread(file)

    control <- data$`Experiment target` == "Non-specific target control-human"
    data$`Experiment target` <- gsub("\\-.*", "", data$`Experiment target`)
    data$`Experiment target`[control] <- "Non-specific target control-human"

    if (!is.null(cellLine))
        data <- data[data$`Biosample term name` %in% cellLine]

    if (!is.null(lab))
        data <- data[data$Lab %in% lab]

    if (!is.null(outputType))
        data <- data[data$`Output type` %in% outputType]

    return(data)
}

#' Load data from L1000 files
#'
#' @importFrom data.table fread setnames
loadL1000file <- function(cell_line, timepoint, folder) {
    file <- file.path(folder, sprintf(
        "L1000_condition_%sh_%s_zscores_Diff_expr.txt", timepoint, cell_line))
    header <- read.table(file, header = TRUE, nrow = 1)
    perturbationList <- fread(file, skip=1, header=FALSE)
    setnames(perturbationList, c("Gene", colnames(header)))
    return(perturbationList)
}

#' @importFrom data.table fread
loadENCODEsample <- function (data, gene, replicate, control=FALSE) {
    rep    <- data$`Biological replicate(s)` == replicate
    target <- data$`Experiment target` %in% gene

    data <- data[rep & target]
    if (control) {
        sample <- sapply(strsplit(data$`Control Experiments`, ", "),
                         `[`, replicate)
    } else {
        sample <- data$`File accession`
    }
    sample <- paste0(sample)

    outfile <- paste0(sample, ".tsv")
    if (!file.exists(paste0(sample, ".tsv"))) {
        file <- sprintf(
            "https://www.encodeproject.org/files/%s/@@download/%s.tsv",
            sample, sample)
        download.file(file, outfile)
    }

    fread(outfile)
}

#' Load an ENCODE experiment for a given gene
loadENCODEexperiment <- function(data, gene) {
    table_rep1     <- loadENCODEsample(data, gene, replicate=1)
    table_rep2     <- loadENCODEsample(data, gene, replicate=2)
    table_control1 <- loadENCODEsample(data, gene, replicate=1, control=TRUE)
    table_control2 <- loadENCODEsample(data, gene, replicate=2, control=TRUE)

    # Check if transcripts are identical across samples
    sameTranscriptsAcrossSamples <- all(sapply(lapply(
        list(table_rep1, table_rep2, table_control1, table_control2),
        "[[", "transcript_id(s)"), identical, table_rep1$`transcript_id(s)`))
    if (!all(sameTranscriptsAcrossSamples))
        stop("Not all samples share the same transcript identifiers")

    # Merge gene counts from the different samples to a single table
    countTable <- cbind(table_rep1[ , c(1:2, 5)], table_rep2[ , 5],
                        table_control1[ , 5], table_control2[ , 5])
    names(countTable)[3:6] <- c("shRNA1", "shRNA2", "control1", "control2")
    rownames(countTable) <- countTable$gene_id
    return(countTable)
}

#' Perform differential gene expression
#'
#' @importFrom stats model.matrix aggregate
#' @importFrom limma voom lmFit eBayes topTable
performDifferentialExpression <- function(counts, geneAnnot) {
    counts <- data.frame(counts)
    rownames(counts) <- counts$gene_id

    # Design matrix
    Sample_info <- data.frame(
        sample = c("shRNA1", "shRNA2", "control1", "control2"),
        condition = c("shRNA", "shRNA", "control", "control"))
    design <- model.matrix(~ condition, Sample_info)
    rownames(design) <- Sample_info$sample

    # Check: identical(names(counts[ , 3:6]), rownames(design_matrix))
    voom <- voom(counts[ , 3:6], design=design, plot=FALSE,
                 normalize.method="quantile")

    # Fit linear model
    fit     <- lmFit(voom[ , colnames(voom$E) %in% rownames(design)],
                     design=design)
    ebayes  <- eBayes(fit)
    results <- topTable(ebayes, coef = 2, number = nrow(ebayes),
                        sort.by = "logFC", resort.by = "p")

    # Convert to gene symbol
    results$Gene_symbol <- geneAnnot$`Associated Gene Name`[
        match(sapply(strsplit(rownames(results), "\\."), `[`, 1),
              geneAnnot$`Ensembl Gene ID`)]

    # Mean-aggregation per gene symbol to compare unique gene knockdowns
    results2 <- aggregate(results[ , 1:6], data=results, FUN=mean,
                          by=list(Gene_symbol=results$Gene_symbol))
    return(results2)
}

#' Calculate statistics
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
#' @importFrom pbapply pblapply
#' @importFrom data.table data.table
correlatePerCellLine <- function(cellLine, diffExprGenes, perturbations,
                                 method) {
    cat(paste("Comparing with cell line", cellLine), fill=TRUE)
    perturbation <- perturbations[[cellLine]]

    # Select intersecting genes
    genes <- intersect(names(diffExprGenes), perturbation$Gene)
    diffExprGenes <- diffExprGenes[genes]

    setkeyv(perturbation, "Gene")
    ref  <- perturbation[genes, ]

    # Suppress warnings to avoid "Cannot compute exact p-value with ties"
    cors <- suppressWarnings(
        pblapply(seq(ncol(ref))[-1],
                 function(k) cor.test(ref[[k]], diffExprGenes, method=method)))

    cor <- sapply(cors, "[[", "estimate")
    names(cor) <- colnames(perturbation)[-c(1)]

    pval <- sapply(cors, "[[", "p.value")
    names(pval) <- colnames(perturbation)[-c(1)]

    res <- data.table(names(cor), cor, pval)
    names(res) <- c("genes",
                    sprintf("%s_t_%s_coef",   cellLine, method),
                    sprintf("%s_t_%s_pvalue", cellLine, method))
    return(res)
}

#' Perform gene set enrichment (GSA) per cell line
#'
#' @importFrom fgsea fgsea
#' @importFrom data.table data.table
#' @importFrom pbapply pblapply
#' @importFrom plyr rbind.fill
performGSAperCellLine <- function(cellLine, countTable, perturbations, gsc) {
    perturbation  <- perturbations[[cellLine]]

    # Run GSA with top 150 genes
    performGSAwithPerturbationSignature <- function(k, perturbation, gsc) {
        signature <- perturbation[[k]]
        names(signature) <- perturbation$Gene
        fgsea(pathways = gsc$gsc, stats = sort(signature),
              minSize = 15, maxSize = 500, nperm = 1)
    }

    cat("Performing GSA using perturbation signatures...", fill=TRUE)
    gsaRes <- pblapply(2:ncol(perturbation),
                       performGSAwithPerturbationSignature, perturbation, gsc)
    # gsaRes <- plyr::compact(gsaRes) # in case of NULL elements
    names(gsaRes) <- names(perturbation)[2:ncol(perturbation)]

    gsaRes_result <- cbind(genes=rep(names(gsaRes), each=2), rbind.fill(gsaRes))

    # based on L1000 paper (page e8) - weighted connectivity score (WTCS)
    ES_list <- c()
    for (geneID in unique(gsaRes_result$genes)) {
        geneID_GSA  <- gsaRes_result[gsaRes_result$genes %in% geneID, ]
        topGenes    <- geneID_GSA$ES[grepl("top", geneID_GSA$pathway)]
        bottomGenes <- geneID_GSA$ES[grepl("bottom", geneID_GSA$pathway)]

        if (sign(topGenes) == sign(bottomGenes))
            ES_final=0
        else if (sign(topGenes) != sign(bottomGenes))
            ES_final = (topGenes - bottomGenes) / 2
        ES_list <- c(ES_list, ES_final)
    }

    names(ES_list) <- unique(gsaRes_result[["genes"]])
    results <- data.table(names(ES_list), ES_list)
    names(results)[1:2] <- c("genes", paste0(cellLine, "_WTCS"))
    return(results)
}

#' Compare against L1000 datasets
#'
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

    browser()
    data <- merged[rowSums(is.na(merged)) != ncol(merged) - 1, ]

    # Add mean calculated across cell lines
    data[ , paste0("Average", colnameSuffix)] <- rowMeans(
        data[ , grep(colnameSuffix, names(data)), with=FALSE], na.rm=TRUE)
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
