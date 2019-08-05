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
loadCTRPgeneExpression <- function(geneExpressionFile="CTRP/geneExpr.txt",
                                   geneInfoFile="CTRP/geneInfo.txt",
                                   cellLineInfoFile="CTRP/cellLineInfo.txt") {
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
    drugSensitivityFile="CTRP/drugSensitivity.txt",
    experimentFile="CTRP/experimentInfo.txt",
    compoundFile="CTRP/compoundInfo.txt") {

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
    colnames(wide)[1] <- "cellLine"

    # Correctly set compound names
    message(sprintf("Loading CTRP compound metadata from %s...", compoundFile))
    compoundInfo <- loadCTRPcompoundInfo(compoundFile)
    attr(wide, "compoundInfo") <- compoundInfo
    return(wide)
}

#' @rdname loadCTRPgeneExpression
#' @keywords internal
loadCTRPcompoundInfo <- function(compoundFile="CTRP/compoundInfo.txt") {
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
    df  <- data.table(mat)

    # Transpose data
    trans <- cbind(colnames(mat), transpose(df))
    colnames(trans) <- c("cellLine", drugSensitivity[[1]])
    attr(trans, "drugActivityMetric") <- "area under dose-response curve (AUC)"

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
loadGDSCfile <- function(file, filename, type, ...) {
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
loadGDSCcellLineInfo <- function(file="GDSC/cellLineInfo.xlsx") {
    cellLineInfo <- loadGDSCfile(file, "Cell_Lines_Details.xlsx",
                                 "cell line information", sheet=2)
    colnames(cellLineInfo)[2] <- "cellLine"
    return(cellLineInfo)
}

#' @rdname loadCTRPgeneExpression
loadGDSCcompoundInfo <- function(file="GDSC/compoundInfo.xlsx") {
    compoundInfo <- loadGDSCfile(file, "Screened_Compounds.xlsx",
                                 "compound information")
    cols <- c("DRUG_ID"="id", "DRUG_NAME"="name")
    colnames(compoundInfo)[match(names(cols), colnames(compoundInfo))] <- cols
    colnames(compoundInfo) <- gsub("_", " ", tolower(colnames(compoundInfo)))
    return(compoundInfo)
}

#' @rdname loadCTRPgeneExpression
#' @importFrom data.table fread transpose
loadGDSCgeneExpression <- function(file="GDSC/geneExpr.txt") {
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
loadGDSCdrugSensitivity <- function(file="GDSC/drugs.xlsx") {
    drugSensitivity <- loadGDSCfile(file, "v17.3_fitted_dose_response.xlsx",
                                    "drug sensitivity")
    # Convert to wide format: compound (rows) vs cell lines (columns)
    drugSensitivity$minus_LN_IC50 <- -drugSensitivity$LN_IC50
    wide <- dcast(drugSensitivity, COSMIC_ID ~ DRUG_ID, fun.aggregate=mean,
                  value.var="minus_LN_IC50")
    colnames(wide)[1] <- "cellLine"
    attr(wide, "drugActivityMetric") <- "-log(IC50)"
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
#' @details If path direct to non-existing files, respective data will be
#'   downloaded.
#'
#' @return Correlation matrix between gene expression and drug sensitivity
#' @keywords internal
prepareExpressionDrugSensitivityAssociation <- function(
    dataset=c("GDSC", "CTRP", "NCI60"), method="spearman") {

    if (dataset == "GDSC") {
        geneExpr        <- loadGDSCgeneExpression()
        geneInfo        <- NULL
        cellLineInfo    <- loadGDSCcellLineInfo()
        drugSensitivity <- loadGDSCdrugSensitivity()
        compoundInfo    <- loadGDSCcompoundInfo()
    } else if (dataset == "CTRP") {
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
    attr(cor, "source")       <- dataset
    attr(cor, "drugActivityMetric") <- attr(drugSensitivity,
                                            "drugActivityMetric")
    return(cor)
}

# Load gene expression/drug sensitivity associations ---------------------------

#' Load gene expression and drug sensitivity correlation matrix
#'
#' @param source Character: source (\code{CTRP}, \code{GDSC} or \code{NCI60})
#' @param file Character: filepath to gene expression and drug sensitivity
#'   association dataset (automatically downloaded if file does not exist)
#'
#' @family functions related with the prediction of targeting drugs
#' @return Correlation matrix between gene expression (rows) and drug
#'   sensitivity (columns)
#' @export
#'
#' @examples
#' loadExpressionDrugSensitivityAssociation("GDSC")
loadExpressionDrugSensitivityAssociation <- function(source, file=NULL) {
    source <- match.arg(source, c("GDSC", "CTRP", "NCI60"))

    if (source == "GDSC") {
        url <- "5q0dazbtnpojw2m/expressionDrugSensitivityCorGDSC.rds"
    } else if (source == "CTRP") {
        url <- "zj53pxwiwdwo133/expressionDrugSensitivityCorCTRP.rds"
    } else if (source == "NCI60") {
        url <- "p6596ym2f08zroh/expressionDrugSensitivityCorNCI60.rds"
    }
    link <- sprintf("https://www.dropbox.com/s/%s?raw=1", url)

    if (is.null(file))
        file <- sprintf("expressionDrugSensitivityCor%s.rds", source)

    downloadIfNotFound(link, file)
    message(sprintf("Loading data from %s...", file))
    cor <- readRDS(file)
    return(cor)
}

# Predict targeting drugs ------------------------------------------------------

#' Predict targeting drugs
#'
#' Identify compounds that may target the phenotype associated with a
#' user-provided differential expression profile by comparing such against a
#' correlation matrix of gene expression and drug sensitivity.
#'
#' @inheritParams compareAgainstReference
#' @param expressionDrugSensitivityCor Matrix: correlation matrix of gene
#'   expression (rows) and drug sensitivity (columns) across cell lines.
#'   Pre-prepared gene expression and drug sensitivity associations are
#'   available to download using
#'   \code{\link{loadExpressionDrugSensitivityAssociation}}.
#'
#' @importFrom pbapply pbapply
#'
#' @inheritSection rankSimilarPerturbations GSEA score
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
#' gdsc <- loadExpressionDrugSensitivityAssociation("GDSC")
#'
#' # Predict targeting drugs on a differential expression profile
#' predictTargetingDrugs(diffExprStat, gdsc)
predictTargetingDrugs <- function(diffExprGenes, expressionDrugSensitivityCor,
                                  method=c("spearman", "pearson", "gsea"),
                                  geneSize=150) {
    cellLines <- length(attr(expressionDrugSensitivityCor, "cellLines"))

    # Check if drug metric is proportional to sensitivity
    isDrugMetricProportionalToSensitivity <- attr(
        expressionDrugSensitivityCor, "isDrugActivityProportionalToSensitivity")
    drugActivityMetric <- attr(expressionDrugSensitivityCor,
                               "drugActivityMetric")

    if (is.null(isDrugMetricProportionalToSensitivity) &&
        !is.null(drugActivityMetric)) {
        if (grepl("-log(IC50)", drugActivityMetric, fixed=TRUE)) {
            isDrugMetricProportionalToSensitivity <- TRUE
        } else if (grepl("-log10(GI50)", drugActivityMetric, fixed=TRUE)) {
            isDrugMetricProportionalToSensitivity <- TRUE
        } else if (grepl("area under dose-response curve",
                         drugActivityMetric)) {
            isDrugMetricProportionalToSensitivity <- FALSE
        }
    }

    if (is.null(isDrugMetricProportionalToSensitivity)) {
        stop("Attribute 'isDrugActivityProportionalToSensitivity' for argument",
             " 'expressionDrugSensitivityCor' cannot be NULL.")
    }

    rankedDrugs <- compareAgainstReference(
        diffExprGenes, expressionDrugSensitivityCor, method=method,
        geneSize=geneSize, cellLines=cellLines, cellLineMean=FALSE,
        rankByAscending=isDrugMetricProportionalToSensitivity)
    colnames(rankedDrugs)[[1]] <- "compound"

    # Inherit input settings and other relevant information
    attr(rankedDrugs, "expressionDrugSensitivityCor") <- attributes(
        expressionDrugSensitivityCor)
    class(rankedDrugs) <- c("targetingDrugs", class(rankedDrugs))
    return(rankedDrugs)
}
