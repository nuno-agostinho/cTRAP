# Data loading from CTRP -------------------------------------------------------

#' Load CTRP data
#'
#' If given paths direct to non-existing files, those files will be downloaded
#'
#' @param geneExpressionFile Character: path to file with gene expression
#' @param geneInfoFile Character: path to file with gene information
#' @param cellLineInfoFile Character: path to file with cell line information
#'
#' @importFrom data.table fread
#' @importFrom reshape2 dcast
#' @importFrom stats setNames
#' @importFrom R.utils renameFile
#'
#' @return Data frame
#' @keywords internal
loadCTRPgeneExpression <- function(
    geneExpressionFile="CTRP 2.1/geneExpr.txt",
    geneInfoFile="CTRP 2.1/geneInfo.txt",
    cellLineInfoFile="CTRP 2.1/cellLineInfo.txt") {

    if (!any(file.exists(geneExpressionFile, geneInfoFile))) {
        url <- file.path("ftp://caftpd.nci.nih.gov/pub/OCG-DCC/CTD2/Broad",
                         "CTRPv2.1_2016_pub_NatChemBiol_12_109",
                         "CTRPv2.1_2016_pub_NatChemBiol_12_109.zip")
        folder <- dirname(geneExpressionFile)
        toExtract <- c("v21.data.gex_avg_log2.txt",
                       "v21.meta.gex_features.txt",
                       "v21.meta.per_cell_line.txt")
        downloadIfNotFound(url, file.path(folder, basename(url)),
                           toExtract=toExtract)

        renameFile(file.path(folder, toExtract[[1]]), geneExpressionFile)
        renameFile(file.path(folder, toExtract[[2]]), geneInfoFile)
        renameFile(file.path(folder, toExtract[[3]]), cellLineInfoFile)
    }
    # Prepare gene expression (in wide format)
    message(sprintf("Loading CTRP gene expression from %s...",
                    geneExpressionFile))
    geneExpr <- fread(geneExpressionFile)
    wide     <- dcast(geneExpr, master_ccl_id ~ idx_gene_feature,
                      fun.aggregate=mean, value.var="mrna_expression_avg_log2")
    colnames(wide)[1] <- "cellLine"

    # Correctly set gene names
    message(sprintf("Loading CTRP gene metadata from %s...", geneInfoFile))
    geneInfo  <- fread(geneInfoFile)
    geneNames <- setNames(geneInfo$gene_primary_name, geneInfo$idx_gene_feature)
    colnames(wide)[-1]     <- geneNames[colnames(wide)[-1]]
    attr(wide, "geneInfo") <- geneInfo[geneInfo$gene_primary_name %in%
                                           colnames(wide)[-1], ]

    # Cell line metadata
    cellLineInfo <- fread(cellLineInfoFile)
    colnames(cellLineInfo)[1]  <- "cellLine"
    attr(wide, "cellLineInfo") <- cellLineInfo
    return(wide)
}

#' @rdname loadCTRPgeneExpression
#'
#' @param drugSensitivityFile Character: path to file with drug sensitivity
#' @param experimentFile Character: path to file with experiment information
#' @param compoundFile Character: path to file with compound information
#'
#' @importFrom data.table fread
#' @importFrom reshape2 dcast
#' @importFrom stats setNames
#' @importFrom R.utils renameFile
loadCTRPdrugSensitivity <- function(
    drugSensitivityFile="CTRP 2.1/drugSensitivity.txt",
    experimentFile="CTRP 2.1/experimentInfo.txt",
    compoundFile="CTRP 2.1/compoundInfo.txt") {

    if (!any(file.exists(drugSensitivityFile, experimentFile, compoundFile))) {
        url <- paste0("ftp://caftpd.nci.nih.gov/pub/OCG-DCC/CTD2/Broad/",
                      "CTRPv2.0_2015_ctd2_ExpandedDataset/",
                      "CTRPv2.0_2015_ctd2_ExpandedDataset.zip")
        folder <- dirname(drugSensitivityFile)
        toExtract <- c("v20.data.curves_post_qc.txt",
                       "v20.meta.per_experiment.txt",
                       "v20.meta.per_compound.txt")
        downloadIfNotFound(url, file.path(folder, basename(url)),
                           toExtract=toExtract)

        renameFile(file.path(folder, toExtract[[1]]), drugSensitivityFile)
        renameFile(file.path(folder, toExtract[[2]]), experimentFile)
        renameFile(file.path(folder, toExtract[[3]]), compoundFile)
    }
    message(sprintf("Loading CTRP drug sensitivity from %s...",
                    drugSensitivityFile))
    drugSensitivity <- fread(drugSensitivityFile)

    # Get cell lines from experiment file
    message(sprintf("Loading CTRP experiment metadata from %s...",
                    experimentFile))
    experimentInfo        <- fread(experimentFile)
    experimentInfo$run_id <- experimentInfo$experiment_date <- NULL

    merged <- merge(drugSensitivity, unique(experimentInfo))
    wide   <- dcast(merged, master_ccl_id ~ master_cpd_id, fun.aggregate=mean,
                    value.var="area_under_curve")
    attr(wide, "drugActivityMetric") <- "z-score(-log10(GI50))"
    attr(wide, "isDrugActivityDirectlyProportionalToSensitivity") <- TRUE
    colnames(wide)[1] <- "cellLine"

    # Correctly set compound names
    message(sprintf("Loading CTRP compound metadata from %s...", compoundFile))
    compoundInfo <- loadCTRPcompoundInfo(compoundFile)
    attr(wide, "compoundInfo") <- compoundInfo
    return(wide)
}

#' @rdname loadCTRPgeneExpression
#' @keywords internal
loadCTRPcompoundInfo <- function(compoundFile="CTRP 2.1/compoundInfo.txt") {
    compoundInfo <- fread(compoundFile)
    cols <- c("master_cpd_id"="id",
              "cpd_name"="name",
              "broad_cpd_id"="broad id",
              "cpd_smiles"="SMILES",
              "cpd_status"="FDA status",
              "target_or_activity_of_compound"="MOA",
              "gene_symbol_of_protein_target"="target")
    colnames(compoundInfo)[match(names(cols), colnames(compoundInfo))] <- cols
    colnames(compoundInfo) <- gsub("_", " ", colnames(compoundInfo))
    return(compoundInfo)
}

# Data loading from NCI60 ------------------------------------------------------

#' @rdname loadCTRPgeneExpression
#' @importFrom readxl read_excel
#' @importFrom data.table data.table transpose
#' @importFrom R.utils renameFile
loadNCI60geneExpression <- function(file="NCI60/geneExpr.xls",
                                    cellLineInfoFile="cellLineInfo.xls") {
    if (!file.exists(file)) {
        url <- paste0(
            "https://discover.nci.nih.gov/cellminerdata/normalizedArchives/",
            "nci60_RNA__RNA_seq_composite_expression.zip")
        folder <- dirname(file)
        downloadIfNotFound(url, file.path(folder, basename(url)))
        renameFile(file.path(folder, "RNA__RNA_seq_composite_expression.xls"),
                   file)
        renameFile(file.path(folder, "NCI60_CELL_LINE_METADATA.xls"),
                   cellLineInfoFile)
    }
    message(sprintf("Loading NCI60 gene expression from %s...", file))
    geneExpr <- read_excel(file, skip=10, na=c("na", "-"))

    # Convert data to numeric
    mat <- data.matrix(geneExpr[7:66])
    df  <- data.table(mat)

    # Transpose data
    trans <- cbind(colnames(mat), transpose(df))
    colnames(trans) <- c("cellLine", geneExpr[[1]])

    # Cell line metadata
    message(sprintf("Loading NCI60 cell line metadata from %s...", file))
    cellLineInfo <- read_excel(cellLineInfoFile, skip=7, n_max=60)
    # Remove footnote letters
    colnames(cellLineInfo) <- gsub(" .(,.)?$", "", colnames(cellLineInfo))
    colnames(cellLineInfo)[1] <- "cellLine"
    cellLineInfo$cellLine <- gsub(" .$", "", cellLineInfo$cellLine)

    # Fix cell names to be the exact same between datasets
    cellLineInfo$cellLine <- gsub("_", "-", cellLineInfo$cellLine)
    fix <- c("BR:HS578T"="BR:HS 578T",
             "BR:T47D"="BR:T-47D",
             "CO:COLO205"="CO:COLO 205",
             "LE:HL-60"="LE:HL-60(TB)",
             "ME:LOXIMVI"="ME:LOX IMVI",
             "LC:A549"="LC:A549/ATCC",
             "OV:NCI-ADR-RES"="OV:NCI/ADR-RES",
             "RE:RXF-393"="RE:RXF 393")
    cellLines <- match(names(fix), cellLineInfo$cellLine)
    cellLineInfo$cellLine[cellLines] <- fix
    colnames(cellLineInfo)[1] <- "cellLine"
    attr(trans, "cellLineInfo") <- cellLineInfo
    return(trans)
}

#' @inherit loadNCI60geneExpression
#'
#' @importFrom readxl read_excel
#' @importFrom data.table data.table transpose
#' @importFrom R.utils renameFile
#'
#' @keywords internal
loadNCI60drugSensitivity <- function(file="NCI60/drugSensitivity.xls") {
    if (!file.exists(file)) {
        url <- paste0("https://discover.nci.nih.gov/cellminerdata/",
                      "normalizedArchives/DTP_NCI60_ZSCORE.zip")
        toExtract <- "DTP_NCI60_ZSCORE.xls"
        folder <- dirname(file)
        downloadIfNotFound(url, file.path(folder, basename(url)))
        renameFile(file.path(folder, toExtract), file)
    }
    message(sprintf("Loading NCI60 drug sensitivity from %s...", file))
    drugSensitivity <- read_excel(file, skip=8, na=c("na", "-"))

    # Convert data to numeric
    mat <- data.matrix(drugSensitivity[7:66])
    mat <- -mat # Additive inverse of AUC values
    df  <- data.table(mat)

    # Transpose data
    trans <- cbind(colnames(mat), transpose(df))
    colnames(trans) <- c("cellLine", drugSensitivity[[1]])
    attr(trans, "drugActivityMetric") <- "-AUC (area under dose-response curve)"
    attr(trans, "isDrugActivityDirectlyProportionalToSensitivity") <- TRUE

    # Prepare compound information
    loadNCI60compoundInfo <- function(drugSensitivity) {
        # NCI60_CIDs: internal table with matches between PubChem CIDs and SIDs
        compoundInfo <- drugSensitivity[1:6]
        compoundInfo$`PubChem CID` <- NCI60_CIDs$CID[match(
            compoundInfo$`PubChem SID`, NCI60_CIDs$SID)]
        colnames(compoundInfo) <- gsub(" [a-z]$", "", colnames(compoundInfo))
        compoundInfo <- compoundInfo[c(1:5, 7, 6)]

        cols <- c("NSC #"="id", "Drug name"="name", "Mechanism of action"="MOA")
        colnames(compoundInfo)[match(names(cols),
                                     colnames(compoundInfo))] <- cols
        return(compoundInfo)
    }
    compoundInfo <- loadNCI60compoundInfo(drugSensitivity)
    attr(trans, "compoundInfo") <- compoundInfo
    return(trans)
}

# Data loading from GDSC -------------------------------------------------------

#' @rdname loadCTRPgeneExpression
#' @importFrom readxl read_excel
#'
#' @param file Character: file path
loadGDSC7file <- function(file, filename, type, ...) {
    if (!file.exists(file)) {
        url <- paste0(
            "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/",
            "release-7.0/", filename)
        downloadIfNotFound(url, file)
    }
    message(sprintf("Loading GDSC %s from %s...", type, file))
    data <- read_excel(file, ...)
    return(data)
}

#' @rdname loadCTRPgeneExpression
loadGDSC7cellLineInfo <- function(file="GDSC_7/cellLineInfo.xlsx") {
    cellLineInfo <- loadGDSC7file(file, "Cell_Lines_Details.xlsx",
                                  "cell line information", sheet=2)
    colnames(cellLineInfo)[2] <- "cellLine"
    return(cellLineInfo)
}

#' @rdname loadCTRPgeneExpression
loadGDSC7compoundInfo <- function(file="GDSC_7/compoundInfo.xlsx") {
    compoundInfo <- loadGDSC7file(file, "Screened_Compounds.xlsx",
                                 "compound information")
    cols <- c("DRUG_ID"="id", "DRUG_NAME"="name")
    colnames(compoundInfo)[match(names(cols), colnames(compoundInfo))] <- cols
    colnames(compoundInfo) <- gsub("_", " ", tolower(colnames(compoundInfo)))
    return(compoundInfo)
}

#' @rdname loadCTRPgeneExpression
#' @importFrom data.table fread transpose
loadGDSC7geneExpression <- function(file="GDSC_7/geneExpr.txt") {
    if (!file.exists(file)) {
        url <- paste0(
            "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/",
            "release-7.0/sanger1018_brainarray_ensemblgene_rma.txt.gz")
        downloadIfNotFound(url, file)
    }
    message(sprintf("Loading GDSC gene expression from %s...", file))
    geneExpr <- fread(file, header=TRUE)

    gext  <- transpose(geneExpr[ , -1])
    gext2 <- cbind(colnames(geneExpr)[-1], gext)
    colnames(gext2) <- c("cellLine", geneExpr[[1]])

    # Convert from ENSEMBL to gene symbols
    colnames(gext2) <- convertENSEMBLtoGeneSymbols(colnames(gext2))
    return(gext2)
}

#' @rdname loadCTRPgeneExpression
#' @importFrom reshape2 dcast
loadGDSC7drugSensitivity <- function(file="GDSC_7/drugs.xlsx") {
    drugSensitivity <- loadGDSC7file(file, "v17.3_fitted_dose_response.xlsx",
                                     "drug sensitivity")
    # Convert to wide format: compound (rows) vs cell lines (columns)
    drugSensitivity$minus_LN_IC50 <- -drugSensitivity$LN_IC50
    wide <- dcast(drugSensitivity, COSMIC_ID ~ DRUG_ID, fun.aggregate=mean,
                  value.var="minus_LN_IC50")
    colnames(wide)[1] <- "cellLine"
    attr(wide, "drugActivityMetric") <- "-log(IC50)"
    attr(wide, "isDrugActivityDirectlyProportionalToSensitivity") <- TRUE
    return(wide)
}

# Correlate gene expression and drug sensitivity -------------------------------

#' @importFrom stats cor
correlateGEandDrugSensitivity <- function(geneExpr, drugSensitivity,
                                          method="spearman") {
    message("Correlating gene expression and drug sensitivity...")
    # Match cell lines between datasets
    cellLines       <- intersect(geneExpr$cellLine, drugSensitivity$cellLine)
    geneExpr        <- geneExpr[match(cellLines, geneExpr$cellLine), ]
    drugSensitivity <- drugSensitivity[match(cellLines,
                                             drugSensitivity$cellLine), ]
    time <- Sys.time()
    res <- cor(geneExpr[ , -1], drugSensitivity[ , -1], method=method,
               use="pairwise.complete.obs")
    diffTime <- format(round(Sys.time() - time, 2))
    msg <- "Gene expression and drug sensitivity correlation performed in %s\n"
    message(sprintf(msg, diffTime))

    attr(res, "cellLines") <- cellLines
    attr(res, "method")    <- method
    attr(res, "type")      <- "compounds"
    attr(res, "date")      <- Sys.time()
    attr(res, "runtime")   <- diffTime
    return(res)
}

#' Prepare gene expression and drug sensitivity correlation matrix
#'
#' @param dataset Character: dataset to use (\code{CTRP}, \code{GDSC} or
#'   \code{NCI60})
#' @param method Character: correlation method to use between gene expression
#'   and drug sensitivity
#'
#' @details If path directs to non-existing files, data will be downloaded.
#'
#' @return Correlation matrix between gene expression and drug sensitivity
#' @keywords internal
prepareExpressionDrugSensitivityAssociation <- function(
    dataset=c("GDSC 7", "CTRP 2.1", "NCI60"), method="spearman") {

    dataset <- match.arg(dataset)

    source <- dataset
    if (dataset == "GDSC 7") {
        geneExpr        <- loadGDSC7geneExpression()
        geneInfo        <- NULL
        cellLineInfo    <- loadGDSC7cellLineInfo()
        drugSensitivity <- loadGDSC7drugSensitivity()
        compoundInfo    <- loadGDSC7compoundInfo()
    } else if (dataset == "CTRP 2.1") {
        geneExpr        <- loadCTRPgeneExpression()
        geneInfo        <- attr(geneExpr, "geneInfo")
        cellLineInfo    <- attr(geneExpr, "cellLineInfo")
        drugSensitivity <- loadCTRPdrugSensitivity()
        compoundInfo    <- attr(drugSensitivity, "compoundInfo")
    } else if (dataset == "NCI60") {
        geneExpr        <- loadNCI60geneExpression()
        geneInfo        <- NULL
        cellLineInfo    <- attr(geneExpr, "cellLineInfo")
        drugSensitivity <- loadNCI60drugSensitivity()
        compoundInfo    <- attr(drugSensitivity, "compoundInfo")
    }
    cor <- correlateGEandDrugSensitivity(geneExpr, drugSensitivity, method)

    # Include metadata
    cellLines    <- intersect(geneExpr$cellLine, drugSensitivity$cellLine)
    cellLineInfo <- cellLineInfo[cellLineInfo$cellLine %in% cellLines, ]
    attr(cor, "cellLineInfo") <- cellLineInfo

    attr(cor, "geneInfo")     <- geneInfo
    attr(cor, "compoundInfo") <- compoundInfo
    attr(cor, "source")       <- source
    attr(cor, "drugActivityMetric") <- attr(drugSensitivity,
                                            "drugActivityMetric")
    attr(cor, "isDrugActivityDirectlyProportionalToSensitivity") <- attr(
        drugSensitivity, "isDrugActivityDirectlyProportionalToSensitivity")
    return(cor)
}

#' @keywords internal
#' @importFrom utils askYesNo
#' @importFrom rhdf5 h5createFile h5write h5createDataset
writeExpressionDrugSensitivityCorHDF5 <- function(
    data, filename="expressionDrugSensitivityCorNCI60.h5") {
    # Create file with given name
    if (file.exists(filename)) {
        continue <- askYesNo(paste(filename, "already exists. Overwrite file?"),
                             default=FALSE)
        if (continue) {
            message("Overwriting file...")
        } else {
            return(invisible(FALSE))
        }
    }
    time <- Sys.time()
    unlink(filename)
    h5createFile(filename)

    attrs <- attributes(data)
    classes <- sapply(sapply(attributes(data), class), "[[", 1)
    attr(data, "classes") <- unname(cbind(names(classes), classes))

    # Save data.frame as datasets (data.frame attributes are unsupported)
    attrsDF <- names(attrs)[sapply(attrs, is.data.frame)]
    attributes(data)[attrsDF] <- NULL
    for (k in attrsDF) h5write(attrs[[k]], filename, k)

    # Save dim names as datasets (large character attributes are unsupported)
    h5write(colnames(data), filename, "colnames")
    h5write(rownames(data), filename, "rownames")
    attr(data, "dimnames") <- NULL

    # Convert logical values to 1 or 0 (logical attributes are unsupported)
    attributes(data)[sapply(attributes(data), isTRUE)]  <- 1
    attributes(data)[sapply(attributes(data), isFALSE)] <- 0

    # Convert date to string (POSIXt attributes are unsupported)
    attr(data, "date") <- as.character(attr(data, "date"))

    # Split data in appropriate chunks (faster reading of smaller subsets)
    chunk <- c( nrow(data), ceiling(ncol(data) / 10) )
    h5createDataset(filename, "data", dims=dim(data), chunk=chunk)

    # h5write(data, filename, "data", write.attributes=TRUE)
    h5write(data, filename, "data", write.attributes=TRUE)
    message(paste("HDF5 file written in", format(round(Sys.time() - time, 2))))
    return(invisible(TRUE))
}

# Load gene expression/drug sensitivity associations ---------------------------

#' List available gene expression and drug sensitivity correlation matrices
#'
#' @param url Boolean: return download link?
#'
#' @family functions related with the prediction of targeting drugs
#' @return Character vector of available gene expression and drug sensitivity
#'   correlation matrices
#' @export
#'
#' @examples
#' listExpressionDrugSensitivityAssociation()
listExpressionDrugSensitivityAssociation <- function(url=FALSE) {
    options <- c(
        "GDSC 7"="5q0dazbtnpojw2m/expressionDrugSensitivityCorGDSC7.rds",
        "CTRP 2.1"="zj53pxwiwdwo133/expressionDrugSensitivityCorCTRP2.1.rds",
        "NCI60"="20ko9lyyyoilfz6/expressionDrugSensitivityCorNCI60.h5")
    link <- sprintf("https://www.dropbox.com/s/%s?raw=1", options)
    names(link) <- names(options)

    res <- link
    if (!url) res <- names(res)
    return(res)
}

#' @keywords internal
#' @importFrom tibble tibble
#' @importFrom rhdf5 h5ls h5read h5readAttributes
readExpressionDrugSensitivityCorHDF5 <- function(
    filename="expressionDrugSensitivityCorNCI60.h5", rows=NULL, cols=NULL,
    loadValues=FALSE) {
    # Check which datasets are attributes or not
    info   <- h5ls(filename)
    isAttr <- info[["name"]] != "data"

    # Match character rows/cols with respective index
    readDimNames <- function(filename, type, index=NULL) {
        attr <- h5read(filename, type)
        if (!is.null(index)) attr <- attr[index]
        return(as.character(attr))
    }
    matchIndex <- function(dim, type) {
        if (is.character(dim)) {
            dim <- sort(unique(match(dim, readDimNames(filename, type))))
        }
        return(dim)
    }
    rows <- matchIndex(rows, "rownames")
    cols <- matchIndex(cols, "colnames")

    # Read data from file
    name <- "data"
    attrs <- h5readAttributes(filename, name)
    if (loadValues) {
        data <- h5read(filename, name, index=list(rows, cols),
                       read.attributes=TRUE)
    } else {
        data <- filename
        mostattributes(data) <- attrs
    }

    # Convert attributes to original classes
    classes <- attr(data, "classes")
    classes <- setNames(classes[ , 2], classes[ , 1])
    classes <- classes[names(classes) %in% names(attributes(data))]
    attr(data, "classes") <- NULL

    for (k in seq(classes)) {
        this <- classes[[k]]
        name <- names(classes[k])
        if ( this == "character" ) {
            attr(data, name) <- as.character(attr(data, name))
        } else if ( this == "POSIXct" ) {
            attr(data, name) <- as.POSIXct(attr(data, name))
        } else if ( this == "logical" ) {
            attr(data, name) <- as.logical(attr(data, name))
        }
    }

    # Convert COMPOUND (data.frames) to tibbles and add them as attributes
    attrsDF <- info[["name"]][info[["dclass"]] == "COMPOUND"]
    for (i in attrsDF) {
        attr(data, i) <- suppressWarnings(tibble(h5read(filename, i)))
    }

    # Add rownames and colnames as attributes
    attr(data, "compounds") <- readDimNames(filename, "colnames", cols)
    attr(data, "genes")     <- readDimNames(filename, "rownames", rows)
    if (loadValues) {
        colnames(data) <- attr(data, "compounds")
        rownames(data) <- attr(data, "genes")
    }
    return(data)
}

#' Operations on \code{expressionDrugSensitivityAssociation} objects
#' @param x An \code{expressionDrugSensitivityAssociation} object
#' @export
dimnames.expressionDrugSensitivityAssociation <- function(x) {
    if (is.character(x)) {
        res <- list(attr(x, "genes"), attr(x, "compounds"))
    } else {
        res <- NextMethod("dimnames")
    }
    return(res)
}

#' @rdname dimnames.expressionDrugSensitivityAssociation
#' @export
dim.expressionDrugSensitivityAssociation <- function(x) {
    if (is.character(x)) {
        res <- vapply(dimnames(x), length, numeric(1))
    } else {
        res <- NextMethod("dim")
    }
    return(res)
}

#' @rdname dimnames.expressionDrugSensitivityAssociation
#'
#' @param i,j Character or numeric indexes specifying elements to extract
#' @param drop Boolean: coerce result to the lowest possible dimension?
#' @param ... Extra arguments given to other methods
#'
#' @export
`[.expressionDrugSensitivityAssociation` <- function(x, i, j, drop=FALSE, ...) {
    if (is.character(x)) {
        out <- x
        nargs <- nargs() - length(list(...)) - 1
        genes <- subsetDim(i, attr(out, "genes"), nargs, areCols=FALSE)
        drugs <- subsetDim(j, attr(out, "compounds"), nargs, areCols=TRUE)
        attr(out, "genes")     <- genes
        attr(out, "compounds") <- drugs
    } else {
        out <- NextMethod("[", drop=drop)
    }
    return(out)
}

#' Load gene expression and drug sensitivity correlation matrix
#'
#' @param source Character: source of matrix to load; see
#'   \code{\link{listExpressionDrugSensitivityAssociation}}
#' @param file Character: filepath to gene expression and drug sensitivity
#'   association dataset (automatically downloaded if file does not exist)
#' @param rows Character or integer: rows
#' @param cols Character or integer: columns
#' @param loadValues Boolean: load data values? Only available for some sources
#'
#' @family functions related with the prediction of targeting drugs
#' @return Correlation matrix between gene expression (rows) and drug
#'   sensitivity (columns)
#' @export
#'
#' @importFrom tools file_ext
#'
#' @examples
#' gdsc <- listExpressionDrugSensitivityAssociation()[[1]]
#' loadExpressionDrugSensitivityAssociation(gdsc)
loadExpressionDrugSensitivityAssociation <- function(source, file=NULL,
                                                     rows=NULL, cols=NULL,
                                                     loadValues=FALSE) {
    available <- listExpressionDrugSensitivityAssociation(url=TRUE)
    source    <- match.arg(source, names(available))
    link      <- available[source]

    if (is.null(file)) file <- gsub("\\?.*", "", basename(link))
    downloadIfNotFound(link, file)
    message(sprintf("Loading data from %s...", file))
    if (file_ext(file) == "rds") {
        if (is.null(cols)) cols <- TRUE
        if (is.null(rows)) rows <- TRUE
        res <- readRDS(file)
        cor <- res[rows, cols, drop=FALSE]
        attrs <- attributes(res)
        attrs <- attrs[!names(attrs) %in% names(attributes(cor))]
        attributes(cor) <- c(attributes(cor), attrs)
    } else {
        cor <- readExpressionDrugSensitivityCorHDF5(file, rows=rows, cols=cols,
                                                    loadValues=loadValues)
    }
    attr(cor, "filename") <- normalizePath(file)
    class(cor) <- c("expressionDrugSensitivityAssociation", class(cor))
    return(cor)
}

# Predict targeting drugs ------------------------------------------------------

#' Predict targeting drugs
#'
#' Identify compounds that may target the phenotype associated with a
#' user-provided differential expression profile by comparing such against a
#' correlation matrix of gene expression and drug sensitivity.
#'
#' @inheritParams rankAgainstReference
#' @param expressionDrugSensitivityCor Matrix: correlation matrix of gene
#'   expression (rows) and drug sensitivity (columns) across cell lines.
#'   Pre-prepared gene expression and drug sensitivity associations are
#'   available to download using
#'   \code{\link{loadExpressionDrugSensitivityAssociation}()}.
#' @param isDrugActivityDirectlyProportionalToSensitivity Boolean: are the
#'   values used for drug activity directly proportional to drug sensitivity?
#'   See details.
#'
#' @importFrom pbapply pbapply
#'
#' @inheritSection rankSimilarPerturbations GSEA score
#'
#' @details
#'   If \code{isDrugActivityDirectlyProportionalToSensitivity = NULL}, the
#'   attribute \code{isDrugMetricDirectlyProportionalToSensitivity} of
#'   \code{expressionDrugSensitivityCor} is used (see
#'   \code{\link{loadExpressionDrugSensitivityAssociation}()}).
#'
#' @family functions related with the prediction of targeting drugs
#' @return Data table with correlation or GSEA results comparing differential
#'   expression values against gene expression and drug sensitivity associations
#' @export
#'
#' @examples
#' # Example of a differential expression profile
#' data("diffExprStat")
#'
#' # Load expression and drug sensitivity association derived from GDSC data
#' gdsc <- loadExpressionDrugSensitivityAssociation("GDSC 7")
#'
#' # Predict targeting drugs on a differential expression profile
#' predictTargetingDrugs(diffExprStat, gdsc)
predictTargetingDrugs <- function(
    input, expressionDrugSensitivityCor,
    method=c("spearman", "pearson", "gsea"), geneSize=150,
    isDrugActivityDirectlyProportionalToSensitivity=NULL, threads=1) {

    cellLines <- length(attr(expressionDrugSensitivityCor, "cellLines"))

    # Check if drug metric is proportional to sensitivity
    if (is.null(isDrugActivityDirectlyProportionalToSensitivity)) {
        isDrugActivityDirectlyProportionalToSensitivity <- attr(
            expressionDrugSensitivityCor,
            "isDrugActivityDirectlyProportionalToSensitivity")
    }

    if (is.null(isDrugActivityDirectlyProportionalToSensitivity)) {
        stop("attribute 'isDrugActivityProportionalToSensitivity' for argument",
             " 'expressionDrugSensitivityCor' cannot be NULL")
    }

    rankedDrugs <- rankAgainstReference(
        input, expressionDrugSensitivityCor, method=method, geneSize=geneSize,
        cellLines=cellLines, cellLineMean=FALSE, threads=threads,
        rankByAscending=isDrugActivityDirectlyProportionalToSensitivity)
    colnames(rankedDrugs)[[1]] <- "compound"

    # Inherit input settings and other relevant information
    attr(rankedDrugs, "expressionDrugSensitivityCor") <- attributes(
        expressionDrugSensitivityCor)
    class(rankedDrugs) <- c("targetingDrugs", class(rankedDrugs))
    return(rankedDrugs)
}

#' @importFrom ggplot2 ggplot geom_point geom_rug geom_density_2d xlab ylab
#'   theme_bw theme
#'
#' @keywords internal
plotTargetingDrug <- function(x, drug, method=c("spearman", "pearson", "gsea"),
                              genes=c("both", "top", "bottom"),
                              data=NULL, title=NULL) {
    method <- match.arg(method)
    drug <- as.character(drug)

    if (is.null(data)) {
        info <- attr(x, "expressionDrugSensitivityCor")
        data <- loadExpressionDrugSensitivityAssociation(info$source,
                                                         info$filename)
    }
    cor   <- data[ , drug]
    input <- attr(x, "input")

    if (method %in% c("spearman", "pearson")) {
        intersecting <- intersect(names(input), names(cor))
        df <- data.frame(input=input[intersecting], cor=cor[intersecting])

        plot <- ggplot(df, aes_string("input", "cor")) +
            geom_point(alpha=0.1) +
            geom_rug(alpha=0.1) +
            geom_density_2d() +
            xlab("Differentially expressed genes (input)") +
            ylab(drug) +
            theme_bw() +
            theme(legend.position="bottom")
        if (!is.null(title)) {
            plot <- plot +
                ggtitle(title) +
                theme(plot.title = element_text(hjust = 0.5))
        }
    } else if (method == "gsea") {
        geneset <- c(attr(x, "pathways"), # legacy
                     attr(x, "geneset"))
        if (is.null(title)) title <- drug
        plot <- plotGSEA(cor, geneset, genes, title=title)
    }
    return(plot)
}
