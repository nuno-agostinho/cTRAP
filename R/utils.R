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
#' @return Named character vector where names are the input ENSEMBL gene
#'   identifiers and the values are the matching gene symbols
#' @export
convertENSEMBLtoGeneSymbols <- function(genes, dataset="hsapiens_gene_ensembl",
                                        mart="ensembl") {
    .Deprecated("convertGeneIdentifiers")
    
    if (!require("biomaRt")) {
        stop("This function requires the 'biomaRt' package installed")
    }
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

#' Convert gene identifiers
#'
#' @param annotation \code{OrgDb} with genome wide annotation for an organism or
#'   \code{character} with species name to query \code{OrgDb}, e.g.
#'   \code{"Homo sapiens"}
#' @param genes Character: genes to be converted
#' @param key Character: type of identifier used, e.g. \code{ENSEMBL}; read
#' \code{?AnnotationDbi::columns}
#' @param target Character: type of identifier to convert to; read
#' \code{?AnnotationDbi::columns}
#' @param ignoreDuplicatedTargets Boolean: if \code{TRUE}, identifiers that
#' share targets with other identifiers will not be converted
#'
#' @importFrom AnnotationDbi select
#' @importFrom data.table data.table
#' @importFrom AnnotationHub AnnotationHub query
#'
#' @family functions for gene expression pre-processing
#' @return Character vector of the respective targets of gene identifiers. The
#' previous identifiers remain other identifiers have the same target (in case
#' \code{ignoreDuplicatedTargets = TRUE}) or if no target was found.
#' @export
#'
#' @examples
#' genes <- c("ENSG00000012048", "ENSG00000083093", "ENSG00000141510",
#'            "ENSG00000051180")
#' convertGeneIdentifiers(genes)
#' convertGeneIdentifiers(genes, key="ENSEMBL", target="UNIPROT")
#' 
#' # Explicit species name to automatically look for its OrgDb database
#' sp <- "Homo sapiens"
#' genes <- c("ENSG00000012048", "ENSG00000083093", "ENSG00000141510",
#'            "ENSG00000051180")
#' convertGeneIdentifiers(genes, sp)
#'
#' # Alternatively, set the annotation database directly
#' ah <- AnnotationHub::AnnotationHub()
#' sp <- AnnotationHub::query(ah, c("OrgDb", "Homo sapiens"))[[1]]
#' columns(sp) # these attributes can be used to change the attributes
#'
#' convertGeneIdentifiers(genes, sp)
convertGeneIdentifiers <- function(genes, annotation="Homo sapiens",
                                   key="ENSEMBL", target="SYMBOL",
                                   ignoreDuplicatedTargets=TRUE) {
    if (is.character(annotation)) {
        ah <- AnnotationHub()
        annotation <- query(ah, c("OrgDb", annotation))[[1]]
        if (length(annotation) == 0) {
            stop(sprintf("No query found for species '%s'", annotation))
        }
    } else if (!is(annotation, "OrgDb")) {
        stop("Annotation needs to be a 'character' or 'OrgDb' object")
    }

    if (key == "ENSEMBL") {
        # Remove ENSEMBL identifiers
        genesClean <- gsub("\\..*", "", genes)
        # Keep version for gene identifier containing the string "PAR_Y"
        par_y <- grep("PAR", genes)
        genesClean[par_y] <- genes[par_y]
    } else {
        genesClean <- genes
    }

    match <- tryCatch(
        suppressMessages(select(annotation, genesClean, target, key)),
        error=return)

    if (is(match, "error")) return(setNames(genes, genes))
    match <- data.table(match, key=key)

    # Ignore missing values
    match <- match[!is.na(match[[target]]), ]

    # Collapse genes with more than one matching target
    colnames(match)[2] <- "target"
    collapsed <- match[
        , list(target=paste(unique(target), collapse="/")), by=key]

    if (ignoreDuplicatedTargets) {
        # Ignore genes sharing the same target
        geneTargets <- collapsed[["target"]]
        collapsed   <- collapsed[
            !geneTargets %in% unique(geneTargets[duplicated(geneTargets)]), ]
    }

    # Replace identifiers by their matching targets (if possible)
    converted <- collapsed[["target"]][match(genesClean, collapsed[[key]])]
    genes[!is.na(converted)] <- converted[!is.na(converted)]
    names(genes) <- genesClean
    return(genes)
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
