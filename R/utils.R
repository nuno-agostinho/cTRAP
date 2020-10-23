#' Strip non-alpha-numeric characters from a string
#'
#' @param str Character
#'
#' @return Character without non-alphanumeric values
#' @keywords internal
stripStr <- function(str) {
    str <- as.character(str)
    str <- gsub("[^[:alnum:]]", "", str)
    return(str)
}

#' Download data if given file is not found
#'
#' @param file Character: filepath
#' @param link Character: link to download file
#' @param ask Boolean: ask to download file?
#' @param toExtract Character: files to extract (if \code{NULL}, extract all)
#'
#' @importFrom utils download.file askYesNo unzip
#' @importFrom tools file_path_sans_ext file_ext
#' @importFrom R.utils isGzipped gunzip
#'
#' @return Download file if file is not found
#' @keywords internal
downloadIfNotFound <- function(link, file, ask=FALSE, toExtract=NULL) {
    extracted <- file_path_sans_ext(file)
    if (file.exists(extracted)) file <- extracted

    folder <- dirname(file)
    if (!file.exists(file)) {
        if (!dir.exists(folder)) {
            # Create folder based on file path
            message(sprintf("Creating folder %s...", folder))
            dir.create(folder)
        }

        # Warn or ask user about data download
        if (ask) {
            download <- askYesNo(
                paste(file, "not found: download file?"), FALSE)
            if (!download) return(stop(paste(file, "not found")))
        } else {
            message(paste(file, "not found: downloading data..."))
        }

        isBinary <- function(file) {
            formats <- c("gz", "bz2", "xz", "tgz", "zip", "rda", "rds", "RData")
            return(any(file_ext(file) %in% formats))
        }

        # Clean link if data is stored in Dropbox
        processed <- link
        processed <- gsub("\\?raw=1$", "", link)

        if (isGzipped(processed)) {
            if (!isGzipped(file)) file <- paste0(file, ".gz")
            mode <- "wb"
        } else if (isBinary(processed)) {
            mode <- "wb"
        } else {
            mode <- "w"
        }
        download.file(link, file, mode=mode)
    }

    # Extract data if GZ or ZIP
    extractionMsg <- sprintf("Extracting %s...", basename(file))
    if (isGzipped(file)) {
        message(extractionMsg)
        file <- gunzip(file, overwrite=TRUE)
    } else if (grepl("\\.zip$", file)) {
        message(extractionMsg)
        zipped <- file
        file <- unzip(zipped, exdir=folder, junkpaths=TRUE, files=toExtract)
        unlink(zipped)
    }
    return(file)
}

#' Convert ENSEMBL gene identifiers to gene symbols
#'
#' @param genes Character: ENSEMBL gene identifiers
#' @param dataset Character: \code{biomaRt} dataset name
#' @param mart Character: \code{biomaRt} database name
#'
#' @importFrom biomaRt useDataset useMart getBM
#'
#' @return Named character vector where names are the input ENSEMBL gene
#'   identifiers and the values are the matching gene symbols
#' @export
#' @examples
#' convertENSEMBLtoGeneSymbols(c("ENSG00000112742", "ENSG00000130234"))
convertENSEMBLtoGeneSymbols <- function(genes, dataset="hsapiens_gene_ensembl",
                                        mart="ensembl") {
    mart      <- useDataset(dataset, useMart(mart))
    processed <- sapply(strsplit(genes, "\\."), `[`, 1)
    geneConversion <- getBM(
        filters="ensembl_gene_id", values=processed, mart=mart,
        attributes=c("ensembl_gene_id", "hgnc_symbol"))
    converted <- geneConversion$hgnc_symbol[
        match(processed, geneConversion$ensembl_gene_id)]
    converted <- setNames(ifelse(converted != "", converted, genes), genes)
    return(converted)
}
