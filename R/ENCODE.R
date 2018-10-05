#' Get experiments files for a given control
#'
#' @param control Character: control identifier
#' @param table Data frame
#'
#' @return Character vector with respective experiment identifiers
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
#'
#' @importFrom httr content GET
#' @importFrom readr read_tsv
#'
#' @return Data frame containing ENCODE knockdown experiment metadata
#' @export
#' @examples
#' downloadENCODEknockdownMetadata("HepG2", "EIF4G1")
downloadENCODEknockdownMetadata <- function(cellLine=NULL, gene=NULL) {
    # Retrieve metadata for knockdown experiments from ENCODE (JSON format) ----
    cat("Downloading metadata for ENCODE knockdown experiments...", fill=TRUE)
    url <- paste(
        sep="&", "https://www.encodeproject.org/search/?type=Experiment",
        "searchTerm=knock", "limit=all", "frame=object", "format=json")
    jsonMetadata <- content(GET(url))

    # Parse an experiment's control --------------------------------------------
    # E.g. for a given experiment, we will obtain the possible controls with
    # strings such as "/experiments/ENCSR942UNX/"; strings are then parsed to
    # retrieve but the control identifier (in this case, "ENCSR942UNX")
    parseControl <- function(experiment) {
        control <- experiment$"possible_controls"
        control <- gsub("/experiments/(.*)/", "\\1", control)
        return(control)
    }
    control <- sapply(jsonMetadata$"@graph", parseControl)
    names(control) <- sapply(jsonMetadata$"@graph", "[[", "accession")

    # Retrieve metadata for ENCODE knockdown experiments (table format) --------
    url <- paste(
        sep="&", "https://www.encodeproject.org/metadata/type=Experiment",
        "limit=all", "searchTerm=knock/metadata.tsv")
    table <- suppressMessages(read_tsv(url))
    table <- table[table$`Assembly` == "hg19" &
                       table$`Output type` == "gene quantifications" &
                       table$Lab == "ENCODE Processing Pipeline", ]

    # Retrieve gene quantification experiment files per control ----------------
    controlUnique     <- unique(unlist(control))
    controlExperiment <- sapply(controlUnique, getENCODEcontrols, table)
    names(controlExperiment) <- controlUnique
    controlExperiment2  <- sapply(controlExperiment, paste, collapse=", ")

    # Add controls in a table column -------------------------------------------
    controlCollapsed <- sapply(control, paste, collapse=", ")
    table$Control    <- controlCollapsed[table$`Experiment accession`]
    controlAll <- control[table$`Experiment accession`]
    controlAll[sapply(controlAll, length) == 0] <- NA

    # Parse experiment files of each control and add them in a table column ----
    index <- rep(seq(controlAll), sapply(controlAll, length))
    exp   <- controlExperiment2[unlist(controlAll)]
    exp   <- split(exp, index)
    table$`Control Experiments` <- sapply(exp, paste, collapse=", ")

    # Sanitize experiment targets ----------------------------------------------
    control <- table$`Experiment target` == "Non-specific target control-human"
    table$`Experiment target` <- gsub("\\-.*", "", table$`Experiment target`)
    table$`Experiment target`[control] <- "Non-specific target control-human"

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
    downloadIfNeeded(outfile, sprintf(
        "https://www.encodeproject.org/files/%s/@@download/%s.tsv",
        sample, sample), gz=FALSE)
    fread(outfile)
}

#' Load an ENCODE gene expression data for a given gene
#'
#' @param metadata Character: ENCODE metadata
#'
#' @return Data frame containing gene read counts
#' @export
#'
#' @examples
#' cellLine <- "HepG2"
#' gene <- "EIF4G1"
#' ENCODEmetadata <- downloadENCODEknockdownMetadata(cellLine, gene)
#' loadENCODEgeneExpression(ENCODEmetadata)
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
