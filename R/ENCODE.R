#' Load ENCODE metadata
#'
#' @param file Character: file path to ENCODE metadata file
#' @param cellLine Character: cell line
#' @param gene Character: target gene
#' @param lab Character: laboratory
#' @param outputType Character: output type
#'
#' @importFrom data.table fread
#'
#' @export
#' @return ENCODE metadata loaded as a data table
loadENCODEmetadata <- function(file, cellLine=NULL, gene=NULL, lab=NULL,
                               outputType=NULL) {
    data <- fread(file)

    control <- data$`Experiment target` == "Non-specific target control-human"
    data$`Experiment target` <- gsub("\\-.*", "", data$`Experiment target`)
    data$`Experiment target`[control] <- "Non-specific target control-human"

    if (!is.null(cellLine))
        data <- data[tolower(data$`Biosample term name`) %in% tolower(cellLine)]
    if (!is.null(gene)) data <- data[data$`Experiment target` %in% gene]
    if (!is.null(lab)) data <- data[data$Lab %in% lab]
    if (!is.null(outputType)) data <- data[data$`Output type` %in% outputType]
    return(data)
}

#' Load ENCODE sample
#'
#' @param metadata Data frame: ENCODE metadata
#' @param replicate Number: replicate
#' @param control Boolean: load control experiment?
#'
#' @importFrom data.table fread
loadENCODEsample <- function (metadata, replicate, control=FALSE) {
    metadata <- metadata[metadata$`Biological replicate(s)` == replicate]

    if (control) {
        sample <- sapply(strsplit(metadata$`Control Experiments`, ", "),
                         `[`, replicate)
    } else {
        sample <- metadata$`File accession`
    }
    sample <- paste0(sample)

    outfile <- paste0(sample, ".tsv")
    downloadIfNeeded(outfile, sprintf(
        "https://www.encodeproject.org/files/%s/@@download/%s.tsv",
        sample, sample), gz=FALSE)
    fread(outfile)
}

#' Load an ENCODE gene expression data for a given gene
#'
#' @param metadata Character: ENCODE metadata
#' @param gene Character: gene
#'
#' @export
#' @return Data frame containing gene read counts for which ENCODE knockdown
#' experiment of a given gene
loadENCODEgeneExpression <- function(metadata) {
    table_rep1     <- loadENCODEsample(metadata, replicate=1)
    table_rep2     <- loadENCODEsample(metadata, replicate=2)
    table_control1 <- loadENCODEsample(metadata, replicate=1, control=TRUE)
    table_control2 <- loadENCODEsample(metadata, replicate=2, control=TRUE)

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
