#' Get experiments files for a given control
#'
#' @param control Character: control identifier
#' @param table Data frame
#'
#' @return Character vector with respective experiment identifiers
#' @keywords internal
getENCODEcontrols <- function(control, table) {
    sub <- table[table$`Experiment accession` == control, ]
    exp <- sub$`File accession`
    # exp <- paste(exp, collapse=", ")
    return(exp)
}

#' Download metadata for ENCODE knockdown experiments
#'
#' @param cellLine Character: cell line
#' @param gene Character: target gene
#' @param file Character: RDS file with metadata (if file doesn't exist, it will
#' be created)
#'
#' @importFrom httr content GET
#' @importFrom data.table fread
#'
#' @family functions related with using ENCODE expression data
#' @return Data frame containing ENCODE knockdown experiment metadata
#' @export
#' @examples
#' downloadENCODEknockdownMetadata("HepG2", "EIF4G1")
downloadENCODEknockdownMetadata <- function(cellLine=NULL, gene=NULL,
                                            file="ENCODEmetadata.RDS") {
    if (!file.exists(file)) {
        # Retrieve metadata for ENCODE KD experiments (JSON format) ------------
        message("Downloading metadata for ENCODE knockdown experiments...")
        url <- paste(
            sep="&", "https://www.encodeproject.org/search/?type=Experiment",
            "searchTerm=knockdown", "limit=all", "frame=object", "format=json")
        jsonMetadata <- content(GET(url))

        # Parse an experiment's control ----------------------------------------
        # E.g. for a given experiment, we will obtain the possible controls with
        # strings such as "/experiments/ENCSR942UNX/"; strings are then parsed
        # to retrieve but the control identifier (in this case, "ENCSR942UNX")
        parseControl <- function(experiment) {
            control <- experiment$"possible_controls"
            control <- gsub("/experiments/(.*)/", "\\1", control)
            return(control)
        }
        control <- sapply(jsonMetadata$"@graph", parseControl)
        names(control) <- sapply(jsonMetadata$"@graph", "[[", "accession")

        # Retrieve metadata for ENCODE KD experiments (table format) -----------
        url <- paste(
            sep="&", "https://www.encodeproject.org/metadata/type=Experiment",
            "limit=all", "searchTerm=knockdown/metadata.tsv")
        table <- suppressWarnings(suppressMessages(fread(url)))
        table <- table[table$`File assembly` == "hg19" &
                           table$`Output type` == "gene quantifications" &
                           table$Lab == "ENCODE Processing Pipeline", ]

        # Retrieve gene quantification experiment files per control ------------
        controlUnique     <- unique(unlist(control))
        controlExperiment <- sapply(controlUnique, getENCODEcontrols, table)
        names(controlExperiment) <- controlUnique
        controlExperiment2  <- sapply(controlExperiment, paste, collapse=", ")

        # Add controls in a table column ---------------------------------------
        controlCollapsed <- sapply(control, paste, collapse=", ")
        table$Control    <- controlCollapsed[table$`Experiment accession`]
        controlAll <- control[table$`Experiment accession`]
        controlAll[sapply(controlAll, length) == 0] <- NA

        # Parse control-specific experiment files and add in a table column ----
        index <- rep(seq(controlAll), sapply(controlAll, length))
        exp   <- controlExperiment2[unlist(controlAll)]
        exp   <- split(exp, index)
        table$`Control Experiments` <- sapply(exp, paste, collapse=", ")

        # Sanitize experiment targets ------------------------------------------
        term <- "Non-specific target control-human"
        control <- table$`Experiment target` == term
        table$`Experiment target` <- gsub(
            "\\-.*", "", table$`Experiment target`)
        table$`Experiment target`[control] <- term
        saveRDS(table, file)
    } else {
        table <- readRDS(file)
    }

    if (!is.null(gene)) table <- table[table$`Experiment target` == gene, ]
    if (!is.null(cellLine)) table <- table[
        tolower(table$`Biosample term name`) == tolower(cellLine), ]

    return(table)
}

#' Load ENCODE sample
#'
#' @param metadata Data frame: ENCODE metadata
#' @param replicate Number: replicate
#' @param control Boolean: load control experiment?
#'
#' @importFrom data.table fread
#' @return Data table with ENCODE sample data
#' @keywords internal
loadENCODEsample <- function (metadata, replicate, control=FALSE) {
    metadata <- metadata[metadata$`Biological replicate(s)` == replicate, ]

    if (control) {
        sample <- sapply(strsplit(metadata$`Control Experiments`, ", "),
                         `[`, replicate)
    } else {
        sample <- metadata$`File accession`
    }
    sample <- paste0(sample)

    outfile <- paste0(sample, ".tsv")
    link <- sprintf("https://www.encodeproject.org/files/%s/@@download/%s.tsv",
                    sample, sample)
    downloadIfNotFound(link, outfile)
    fread(outfile)
}

#' Load ENCODE samples
#'
#' Samples are automatically downloaded if they are not found in the current
#' working directory.
#'
#' @param metadata Character: ENCODE metadata
#'
#' @importFrom pbapply pblapply
#'
#' @family functions related with using ENCODE expression data
#' @return List of loaded ENCODE samples
#' @export
#'
#' @examples
#' if (interactive()) {
#'   # Load ENCODE metadata for a specific cell line and gene
#'   cellLine <- "HepG2"
#'   gene <- c("EIF4G1", "U2AF2")
#'   ENCODEmetadata <- downloadENCODEknockdownMetadata(cellLine, gene)
#'
#'   # Load samples based on filtered ENCODE metadata
#'   loadENCODEsamples(ENCODEmetadata)
#' }
loadENCODEsamples <- function(metadata) {
    loadENCODEsamplePerGene <- function(metadata) {
        gene <- list()
        reps <- as.numeric(metadata$`Biological replicate(s)`)
        for (rep in reps) {
            sample  <- loadENCODEsample(metadata, replicate=rep)
            control <- loadENCODEsample(metadata, replicate=rep, control=TRUE)
            gene <- c(gene, rep=list(sample), control=list(control))
        }
        names(gene) <- paste0(names(gene), rep(reps, each=max(reps)))
        return(gene)
    }

    metadataPerGene <- split(metadata, sprintf("%s_%s_%s",
                                               metadata$`Biosample term name`,
                                               metadata$`Experiment target`,
                                               metadata$`Experiment accession`))
    res <- pblapply(metadataPerGene, loadENCODEsamplePerGene)
    return(res)
}

#' Load ENCODE gene expression data
#'
#' @param samples List of loaded ENCODE samples
#'
#' @seealso \code{\link{convertENSEMBLtoGeneSymbols}()}
#'
#' @family functions related with using ENCODE expression data
#' @return Data frame containing gene read counts
#' @export
#'
#' @examples
#' if (interactive()) {
#'   # Load ENCODE metadata for a specific cell line and gene
#'   cellLine <- "HepG2"
#'   gene <- "EIF4G1"
#'   ENCODEmetadata <- downloadENCODEknockdownMetadata(cellLine, gene)
#'
#'   # Load samples based on filtered ENCODE metadata
#'   ENCODEsamples <- loadENCODEsamples(ENCODEmetadata)[[1]]
#'
#'   prepareENCODEgeneExpression(ENCODEsamples)
#' }
prepareENCODEgeneExpression <- function(samples) {
    # Check if transcripts are identical across samples
    sameTranscriptsAcrossSamples <- all(sapply(lapply(
        samples, "[[", "transcript_id(s)"),
        identical, samples$rep1$`transcript_id(s)`))
    if (!all(sameTranscriptsAcrossSamples))
        stop("not all samples share the same transcript identifiers")

    # Merge gene counts from the different samples to a single table
    countTable <- cbind(samples$rep1[ , c(1:2, 5)], samples$rep2[ , 5],
                        samples$control1[ , 5], samples$control2[ , 5])
    names(countTable)[3:6] <- c("shRNA1", "shRNA2", "control1", "control2")
    rownames(countTable) <- countTable$gene_id
    return(countTable)
}

#' Perform differential gene expression based on ENCODE data
#'
#' @param counts Data frame: gene expression
#'
#' @importFrom stats model.matrix aggregate
#' @importFrom limma voom lmFit eBayes topTable
#'
#' @family functions related with using ENCODE expression data
#' @return Data frame with differential gene expression results between
#' knockdown and control
#' @export
#'
#' @examples
#' if (interactive()) {
#'   # Download ENCODE metadata for a specific cell line and gene
#'   cellLine <- "HepG2"
#'   gene <- "EIF4G1"
#'   ENCODEmetadata <- downloadENCODEknockdownMetadata(cellLine, gene)
#'
#'   # Download samples based on filtered ENCODE metadata
#'   ENCODEsamples <- loadENCODEsamples(ENCODEmetadata)[[1]]
#'
#'   counts <- prepareENCODEgeneExpression(ENCODEsamples)
#'
#'   # Remove low coverage (at least 10 counts shared across two samples)
#'   minReads   <- 10
#'   minSamples <- 2
#'   filter <- rowSums(counts[ , -c(1, 2)] >= minReads) >= minSamples
#'   counts <- counts[filter, ]
#'
#'   # Convert ENSEMBL identifier to gene symbol
#'   counts$gene_id <- convertENSEMBLtoGeneSymbols(counts$gene_id)
#'
#'   # Perform differential gene expression analysis
#'   diffExpr <- performDifferentialExpression(counts)
#' }
performDifferentialExpression <- function(counts) {
    counts <- data.frame(counts)

    # Design matrix
    Sample_info <- data.frame(
        sample = c("shRNA1", "shRNA2", "control1", "control2"),
        condition = c("shRNA", "shRNA", "control", "control"))
    design <- model.matrix(~ condition, Sample_info)
    rownames(design) <- Sample_info$sample

    # Check: identical(names(counts[ , 3:6]), rownames(design_matrix))
    voom <- voom(counts[ , -c(1, 2)], design=design, plot=FALSE,
                 normalize.method="quantile")

    # Fit linear model
    fit <- lmFit(voom[ , colnames(voom$E) %in% rownames(design)], design=design)
    ebayes <- eBayes(fit)
    results <- topTable(ebayes, coef=2, number=nrow(ebayes),
                        genelist=counts$gene_id)

    # Mean-aggregation per gene symbol to compare unique gene knockdowns
    meanAggr <- aggregate(results[ , -1], data=results, FUN=mean,
                          by=list(Gene_symbol=results$ID))

    # Remove non-matching genes (if any)
    meanAggr           <- meanAggr[meanAggr$Gene_symbol != "", ]
    rownames(meanAggr) <- meanAggr$Gene_symbol
    return(meanAggr)
}
