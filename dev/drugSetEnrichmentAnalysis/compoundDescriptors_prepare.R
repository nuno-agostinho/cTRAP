# Nuno Agostinho (NMorais Lab), Instituto de Medicina Molecular, 07 Aug 2019
# Convert molecular descriptor datasets to processed data in RDS format
library(cTRAP)
library(data.table)
descriptors <- list()

# NCI60 compound 2D descriptors ------------------------------------------------
prepareCompoundDescriptors_NCI60_2D <- function(file) {
    descriptors <- fread(file)
    # Remove last row and last column (only contain NAs)
    descriptors <- descriptors[seq(nrow(descriptors) - 1),
                               seq(ncol(descriptors) - 1), with=FALSE]

    # Remove columns not interesting to be used as molecular descriptors
    cols <- c("Drug_name", "PubChem_SID", "Structure of SMILES [idcode]",
              "SMILES")
    descriptors <- descriptors[ , -cols, with=FALSE]

    # Replace empty strings and values of 99999 with missing values
    descriptors[descriptors == ""] <- NA
    descriptors[descriptors == 99999] <- NA

    # Add compound information from NCI60
    nci60 <- loadExpressionDrugSensitivityAssociation("NCI60")
    compoundInfo <- attr(nci60, "compoundInfo")
    compoundInfo[[1]] <- as.numeric(compoundInfo[[1]])

    descriptors  <- cbind(compound=compoundInfo[[1]], descriptors)
    attr(descriptors, "compoundInfo") <- compoundInfo
    attr(descriptors, "source")       <- "NCI60"
    attr(descriptors, "type")         <- "2D"
    return(descriptors)
}

file_NCI60_2D <- "compound_descriptors_NCI60_2D.txt"
descriptors$NCI60_2D <- prepareCompoundDescriptors_NCI60_2D(file_NCI60_2D)
saveRDS(descriptors$NCI60_2D, paste0(file_path_sans_ext(file_NCI60_2D), ".rds"))

# NCI60 compound 3D descriptors ------------------------------------------------
prepareCompoundDescriptors_NCI60_3D <- function(file) {
    descriptors     <- fread(file)
    descriptors$id  <- NULL
    descriptors$idx <- NULL

    # Match compound information from NCI60 based on name and SMILES
    nci60               <- loadExpressionDrugSensitivityAssociation("NCI60")
    compoundInfo        <- attr(nci60, "compoundInfo")
    compoundInfo$id     <- as.numeric(compoundInfo$id)
    compoundInfo$SMILES <- gsub("\\\\", "", compoundInfo$SMILES)

    merged <- merge(compoundInfo, descriptors, all=TRUE,
                    by.x=c("name", "SMILES"), by.y=c("name", "smiles"))
    merged <- merged[order(merged$id), ]
    browser()

    # Remove rows with mostly missing values
    merged <- merged[rowSums(is.na(merged)) < round(ncol(merged) * 0.9), ]
    merged <- merged[!is.na(merged$id), ]

    # Split between molecular descriptors and compound information
    cols <- colnames(compoundInfo)
    compoundInfo <- merged[ , cols]
    rownames(compoundInfo) <- NULL

    descriptors <- data.table(merged)[ , -cols, with=FALSE]

    # Replace empty strings and values of 99999 with missing values
    # descriptors[descriptors == ""] <- NA
    descriptors[descriptors == 99999] <- NA

    descriptors <- cbind(compound=compoundInfo$id, descriptors)
    attr(descriptors, "compoundInfo") <- compoundInfo
    attr(descriptors, "source")       <- "NCI60"
    attr(descriptors, "type")         <- "3D"
    return(descriptors)
}

file_NCI60_3D <- "compound_descriptors_NCI60_3D_2.txt"
descriptors$NCI60_3D <- prepareCompoundDescriptors_NCI60_3D(file_NCI60_3D)
saveRDS(descriptors$NCI60_3D, paste0(file_path_sans_ext(file_NCI60_3D), ".rds"))

# CMap compound 2D descriptors -------------------------------------------------
prepareCMapCompoundInfo <- function(id) {
    metadata     <- loadCMapData("metadata.txt", "metadata")
    metadata     <- unique(metadata[ , c("pert_id", "pert_iname"), with=FALSE])
    metadata     <- metadata[metadata$pert_id %in% id, ]

    compoundInfo <- loadCMapData("compoundInfo.txt", "compoundInfo")
    compoundInfo <- merge(metadata, compoundInfo, by="pert_iname", all.x=TRUE)
    compoundInfo <- compoundInfo[match(id, compoundInfo$pert_id), ]
    return(compoundInfo)
}

prepareCompoundDescriptors_CMap_2D <- function(file) {
    descriptors <- fread(file)
    descriptors <- descriptors[ , seq(ncol(descriptors) - 1), with=FALSE]

    # Remove columns not interesting to be used as molecular descriptors
    descriptors <- descriptors[ , - c(
        "pert_iname", "pert_type", "inchi_key_prefix", "inchi_key",
        "Structure of canonical_smiles [idcode]", "canonical_smiles",
        "pubchem_cid", "JuanCarlos_Match")]

    # Replace empty strings and values of 99999 with missing values
    descriptors[descriptors == ""] <- NA
    descriptors[descriptors == 99999] <- NA

    # Add compound information
    attr(descriptors,"compoundInfo") <- prepareCMapCompoundInfo(
        descriptors$pert_id)
    attr(descriptors, "source")      <- "CMap"
    attr(descriptors, "type")        <- "2D"
    return(descriptors)
}

file_CMap_2D <- "compound_descriptors_CMap_2D.txt"
descriptors$CMap_2D <- prepareCompoundDescriptors_CMap_2D(file_CMap_2D)
saveRDS(descriptors$CMap_2D, paste0(file_path_sans_ext(file_CMap_2D), ".rds"))

# CMap compound 3D descriptors -------------------------------------------------
prepareCompoundDescriptors_CMap_3D <- function(file) {
    descriptors <- fread(file)
    id <- descriptors$id

    # Remove columns not interesting to be used as molecular descriptors
    descriptors <- descriptors[ , -c("V1", "name", "smiles", "id")]
    descriptors  <- cbind(compound=id, descriptors)

    # Replace empty strings and values of 99999 with missing values
    descriptors[descriptors == ""] <- NA
    descriptors[descriptors == 99999] <- NA

    # Add compound information
    attr(descriptors,"compoundInfo") <- prepareCMapCompoundInfo(id)
    attr(descriptors, "source")      <- "CMap"
    attr(descriptors, "type")        <- "3D"
    return(descriptors)
}

file_CMap_3D <- "compound_descriptors_CMap_3D.txt"
descriptors$CMap_3D <- prepareCompoundDescriptors_CMap_3D(file_CMap_3D)
saveRDS(descriptors$CMap_3D, paste0(file_path_sans_ext(file_CMap_3D), ".rds"))
