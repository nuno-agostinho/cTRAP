# Nuno Agostinho (NMorais Lab), Instituto de Medicina Molecular, 07 Aug 2019
# Convert molecular descriptor datasets to processed data in RDS format

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(cTRAP)
library(data.table)

# NCI60 compound 2D descriptors ------------------------------------------------
prepareCompoundDescriptors_NCI60_2D <- function(file) {
    descriptors <- fread(file)
    descriptors <- descriptors[seq(nrow(descriptors) - 1),
                               seq(ncol(descriptors) - 1), with=FALSE]

    # Remove columns not interesting to be used as molecular descriptors
    descriptors <- descriptors[ , - c("Drug_name", "PubChem_SID",
                                      "Structure of SMILES [idcode]", "SMILES")]

    # Replace empty strings with missing values
    descriptors[descriptors == ""] <- NA

    # Add compound information from NCI60
    nci60 <- loadExpressionDrugSensitivityAssociation("NCI60")
    compoundInfo <- attr(nci60, "compoundInfo")
    descriptors  <- cbind(compound=compoundInfo$id, descriptors)
    attr(descriptors, "compoundInfo") <- compoundInfo
    attr(descriptors, "source")       <- "NCI60"
    attr(descriptors, "type")         <- "2D"
    return(descriptors)
}

file <- "compound_descriptors_NCI60_2D.txt"
descriptors <- prepareCompoundDescriptors_NCI60_2D(file)
saveRDS(descriptors, paste0(file_path_sans_ext(file), ".rds"))

# NCI60 compound 3D descriptors ------------------------------------------------
prepareCompoundDescriptors_NCI60_3D <- function(file) {
    descriptors <- fread(file)
    name <- descriptors[["smiles"]]
    cols <- c("id", "name", "smiles")
    compoundInfo <- descriptors[ , cols, with=FALSE]

    # Remove columns not interesting to be used as molecular descriptors
    descriptors <- descriptors[ , -cols, with=FALSE]

    # Replace empty strings with missing values
    descriptors[descriptors == ""] <- NA

    # Add compound information from NCI60
    # nci60 <- loadExpressionDrugSensitivityAssociation("NCI60")
    # compoundInfo <- attr(nci60, "compoundInfo")
    descriptors  <- cbind(compound=compoundInfo$id, descriptors)
    attr(descriptors, "compoundInfo") <- compoundInfo
    attr(descriptors, "source")       <- "NCI60"
    attr(descriptors, "type")         <- "3D"
    return(descriptors)
}

file <- "compound_descriptors_NCI60_3D_2.txt"
descriptors <- prepareCompoundDescriptors_NCI60_3D(file)
saveRDS(descriptors, paste0(file_path_sans_ext(file), ".rds"))

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

    # Replace empty strings with missing values
    descriptors[descriptors == ""] <- NA

    # Add compound information
    attr(descriptors,"compoundInfo") <- prepareCMapCompoundInfo(
        descriptors$pert_id)
    attr(descriptors, "source")      <- "CMap"
    attr(descriptors, "type")        <- "2D"
    return(descriptors)
}

file <- "compound_descriptors_CMap_2D.txt"
descriptors <- prepareCompoundDescriptors_CMap_2D(file)
saveRDS(descriptors, paste0(file_path_sans_ext(file), ".rds"))

# CMap compound 3D descriptors -------------------------------------------------
prepareCompoundDescriptors_CMap_3D <- function(file) {
    descriptors <- fread(file)
    id <- descriptors$id

    # Remove columns not interesting to be used as molecular descriptors
    descriptors <- descriptors[ , -c("V1", "name", "smiles", "id")]
    descriptors  <- cbind(compound=id, descriptors)

    # Replace empty strings with missing values
    descriptors[descriptors == ""] <- NA

    # Add compound information
    attr(descriptors,"compoundInfo") <- prepareCMapCompoundInfo(id)
    attr(descriptors, "source")      <- "CMap"
    attr(descriptors, "type")        <- "3D"
    return(descriptors)
}

file <- "compound_descriptors_CMap_3D.txt"
descriptors <- prepareCompoundDescriptors_CMap_3D(file)
saveRDS(descriptors, paste0(file_path_sans_ext(file), ".rds"))
