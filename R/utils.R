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

#' Subset rows or columns based on a given index
#' @return Subset rows/columns
#' @keywords internal
subsetDim <- function(k, dims, nargs, areCols=TRUE) {
    hasK <- !missing(k)
    # Allow to search based on characters
    names(dims) <- dims
    if (hasK && nargs == 2) {
        dims <- dims[k]
    } else if (hasK && areCols && nargs == 1) {
        dims <- dims[k]
    }
    if (anyNA(dims)) {
        stop(ifelse(areCols, "columns", "rows"), " out of bounds")
    }
    return(unname(dims))
}

#' Subset data by rows and/or columns
#'
#' @return Subset data
#' @keywords internal
subsetData <- function(x, i, j, rowAttr, colAttr, nargs, ...) {
    nargs <- nargs - length(list(...)) - 1
    rows  <- attr(x, rowAttr)
    rows  <- subsetDim(i, rows, nargs, areCols=FALSE)
    attr(x, rowAttr) <- rows

    # If no j is provided explicitly, replace j with i
    if (missing(j) && nargs == 1) j <- i
    cols <- attr(x, colAttr)
    cols <- subsetDim(j, cols, nargs, areCols=TRUE)
    attr(x, colAttr) <- cols
    return(x)
}

#' Faster version of \code{shiny::HTML}
#'
#' @param text Character: text
#'
#' @return HTML element
#' @keywords internal
HTMLfast <- function(text) {
    attr(text, "html") <- TRUE
    class(text) <- c("html", "character")
    return(text)
}

#' Create word break opportunities (for HTML) using given characters
#'
#' @param str Character: text
#' @param pattern Character: pattern(s) of interest to be used as word break
#' opportunities
#' @param html Boolean: convert to HTML?
#'
#' @importFrom shiny HTML
#'
#' @return String containing HTML elements
#' @keywords internal
prepareWordBreak <- function(str, pattern=c(".", "-", "\\", "/", "_", ",",
                                            " ", "+", "="),
                             html=TRUE) {
    res <- str
    # wbr: word break opportunity
    for (p in pattern) res <- gsub(p, paste0(p, "<wbr>"), res, fixed=TRUE)

    if (html) {
        if (length(res) == 1) {
            res <- HTML(res)
        } else {
            res <- lapply(res, HTMLfast)
        }
    }
    return(res)
}
