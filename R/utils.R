#' Perform differential gene expression based on ENCODE data
#'
#' @param counts Data frame: gene expression
#'
#' @importFrom stats model.matrix aggregate
#' @importFrom limma voom lmFit eBayes topTable
#'
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

#' Download data if given file is not found
#'
#' @param file Character: filepath
#' @param link Character: link to download file
#' @param ask Boolean: ask to download file?
#' @param toExtract Character: files to extract (if \code{NULL}, extract all)
#' @param fixExtension Boolean: fix extension used based on link
#'
#' @importFrom utils download.file askYesNo
#' @importFrom tools file_path_sans_ext
#' @importFrom R.utils isGzipped gunzip
#'
#' @return Download file if file is not found
#' @keywords internal
downloadIfNotFound <- function(link, file, ask=FALSE, toExtract=NULL,
                               forceCorrectExtension=TRUE) {
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

        # Download data
        if (isGzipped(link)) {
            if (!isGzipped(file)) file <- paste0(file, ".gz")
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

#' Parse CMap identifier
#'
#' @param id Character: CMap identifier
#' @param cellLine Boolean: if \code{TRUE}, return cell line information from
#'   CMap identifier; else, return the CMap identifier without the cell line
#'
#' @return Character vector with information from CMap identifiers
#' @export
#'
#' @examples
#' id <- c("CVD001_HEPG2_24H:BRD-K94818765-001-01-0:4.8",
#'         "CVD001_HEPG2_24H:BRD-K96188950-001-04-5:4.3967",
#'         "CVD001_HUH7_24H:BRD-A14014306-001-01-1:4.1")
#' parseCMapID(id, cellLine=TRUE)
#' parseCMapID(id, cellLine=FALSE)
parseCMapID <- function(id, cellLine=FALSE) {
    if (cellLine) {
        # Retrieve cell line
        res <- gsub(".*\\_([A-Z].*)\\_.*", "\\1", id)
        # Assign missing values to identifiers of summarised perturbation scores
        res <- ifelse(grepl(":", res), NA, res)
    } else {
        # Remove cell line identifier
        res <- gsub("\\_[A-Z].*\\_", "\\_", id)
    }
    names(res) <- id
    return(res)
}

# perturbationChanges object ---------------------------------------------------

#' Subset a \code{perturbationChanges} object
#'
#' @param x \code{perturbationChanges} object
#' @param i,j Character or numeric indexes specifying elements to extract
#' @param drop Boolean: coerce result to the lowest possible dimension?
#' @param ... Extra parameters passed to \code{`[`}
#'
#' @return \code{perturbationChanges} object with subset data
#' @export
`[.perturbationChanges` <- function(x, i, j, drop=FALSE, ...) {
    if (is.character(x)) {
        out <- x
        nargs <- nargs() - length(list(...)) - 1

        hasI <- !missing(i)
        hasJ <- !missing(j)
        genes <- attr(out, "genes")
        perts <- attr(out, "perturbations")
        # Allow to search based on characters
        names(genes) <- genes
        names(perts) <- perts

        if (nargs == 2) {
            if (hasI) genes <- genes[i]
            if (hasJ) perts <- perts[j]
        } else if (hasI && nargs == 1) {
            perts <- perts[i]
        }
        if (anyNA(perts) || anyNA(genes)) stop("subscript out of bounds")
        attr(out, "genes") <- unname(genes)
        attr(out, "perturbations") <- unname(perts)
    } else {
        out <- NextMethod("[", drop=drop)
    }

    # Trim metadata to only contain subset information
    attrs <- attributes(x)
    if (!is.null(ncol(out)) && ncol(x) != ncol(out) &&
        !is.null(attrs$metadata)) {
        samples <- attrs$metadata$sig_id %in% colnames(out)
        attr(out, "metadata") <- NULL
        attrs$metadata <- attrs$metadata[samples, , drop=FALSE]
    }
    # Inherit input's attributes
    attrs <- attrs[!names(attrs) %in% names(attributes(out))]
    attributes(out) <- c(attributes(out), attrs)
    return(out)
}

#' Dimensions of a \code{perturbationChanges} object
#'
#' @param x \code{perturbationChanges} object
#'
#' @return Dimensions of a \code{perturbationChanges} object
#' @export
dim.perturbationChanges <- function(x) {
    if (is.character(x)) {
        res <- vapply(dimnames(x), length, numeric(1))
    } else {
        res <- NextMethod("dim")
    }
    return(res)
}

#' Dimnames of a \code{perturbationChanges} object
#'
#' @param x \code{perturbationChanges} object
#'
#' @return Retrieve dimnames of a \code{perturbationChanges} object
#' @export
dimnames.perturbationChanges <- function(x) {
    if (is.character(x)) {
        res <- list(attr(x, "genes"), attr(x, "perturbations"))
    } else {
        res <- NextMethod("dimnames")
    }
    return(res)
}

# similarPerturbations object --------------------------------------------------

#' Print a \code{similarPerturbations} object
#'
#' @param x \code{similarPerturbations} object
#' @param perturbation Character (perturbation identifier) or numeric
#'   (perturbation index)
#' @param ... Extra parameters passed to \code{print}
#'
#' @return Information on \code{perturbationChanges} object or on specific
#'   perturbations (if \code{perturbation} is set)
#' @export
print.similarPerturbations <- function(x, perturbation=NULL, ...) {
    if (is.null(perturbation)) {
        NextMethod("print")
    } else {
        if (is.numeric(perturbation)) perturbation <- x[[1]][perturbation]

        metadata <- attr(x, "metadata")
        if (!is.null(metadata)) {
            selectMetadata <- metadata[metadata$sig_id %in% perturbation]
            if (nrow(selectMetadata) == 0) {
                # Check to see if using identifiers referring to summary stats
                summaryID <- parseCMapID(metadata$sig_id, cellLine=FALSE)
                selectMetadata <- metadata[summaryID %in% perturbation]
            }
        }

        compoundInfo <- attr(x, "compoundInfo")
        if (!is.null(compoundInfo)) {
            compound <- selectMetadata$pert_iname
            selectCompounds <- compoundInfo[
                compoundInfo$pert_iname %in% compound]
            res <- list(metadata=selectMetadata, compoundInfo=selectCompounds)
        } else {
            selectCompounds <- NULL
            res <- list(metadata=selectMetadata)
        }
        return(res)
    }
}

#' Cross Tabulation and Table Creation
#'
#' @param x \code{similarPerturbations} object
#' @param ... Extra parameters passed to \code{table}
#' @param clean Boolean: only show certain columns (to avoid redundancy)?
#'
#' @return Complete table with metadata based on a \code{similarPerturbations}
#'   object
as.table.similarPerturbations <- function(x, ..., clean=TRUE) {
    metadata <- attr(x, "metadata")
    if (!is.null(metadata)) {
        nonCellID <- "non_cell_id"

        summaryID <- parseCMapID(metadata$sig_id, cellLine=FALSE)
        metadata[[nonCellID]] <- summaryID
        metadataSubset <- metadata[unique(match(summaryID, summaryID)), ]
        metadataSubset[ , c("cell_id", "sig_id", "distil_id")] <- NULL

        x[[nonCellID]] <- parseCMapID(x[[1]], cellLine=FALSE)
        res <- merge(x, metadataSubset, all.x=TRUE, by=nonCellID)
        res[[nonCellID]] <- NULL

        compoundInfo <- attr(x, "compoundInfo")
        if (!is.null(compoundInfo)) {
            res <- merge(res, compoundInfo, by="pert_iname", all.x=TRUE)
            # Place "pert_iname" column after "pert_id" one
            m <- match("pert_id", colnames(res))
            res <- res[ , c(2:m, 1, (m+1):(ncol(res))), with=FALSE]
        }
    } else {
        res <- x
    }

    if (clean) {
        hideCols <- c(colnames(res)[endsWith(colnames(res), "value") |
                                        endsWith(colnames(res), "value_rank")],
                      "pert_dose", "pert_dose_unit",
                      "pert_time", "pert_time_unit")
        hideCols <- hideCols[hideCols %in% colnames(res)]
        if (length(hideCols) > 0) res <- res[ , -hideCols, with=FALSE]
    }

    return(res)
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
