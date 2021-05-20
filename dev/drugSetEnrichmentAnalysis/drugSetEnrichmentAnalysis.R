data("diffExprStat")

metadata <- loadCMapData("cmapMetadata.txt")
filtered <- filterCMapMetadata(metadata, perturbationType="Compound")
perts <- prepareCMapPerturbations(
    metadata, "cmapZscores.gctx", "cmapGeneInfo.txt",
    "cmapCompoundInfo_drugs.txt")
# devtools::document(); devtools::load_all()
rankedPerts <- rankSimilarPerturbations(diffExprStat, perts)

compoundInfo <- attr(perts, "compoundInfo")
cmap_pubchem_cid <- as.numeric(compoundInfo$pubchem_cid)
length(unique(cmap_pubchem_cid))

# Drug set enrichment analysis -------------------------------------------------

# Create drug sets from 2D and 3D descriptors
descriptors2D <- data.table::fread("compound_descriptors_NCI60_2D.txt",
                                   na.strings=c("", "NA"))
descriptors2D$V45 <- NULL
descriptors2D <- descriptors2D[seq(nrow(descriptors2D) - 1), ]

# Add NCI60 identifier as first column
nci60_compoundInfo <- attr(loadNCI60drugSensitivity(), "compoundInfo")
descriptors2D <- cbind(nci60_compoundInfo$id, descriptors2D)

# Convert PubChem SIDs to CIDs using https://pubchem.ncbi.nlm.nih.gov/idexchange
conversion <- cTRAP::NCI60_CIDs
descriptors2D$PubChem_CID <- conversion$CID[match(descriptors2D$PubChem_SID,
                                                  conversion$SID)]

predicted <- predictTargetingDrugs(diffExprStat, nci60cor)
drugSets <- prepareDrugSets(descriptors2D, "PubChem_CID")
test <- analyseDrugSetEnrichment(predicted, drugSets)

# CMap intersection ------------------------------------------------------------

# # 233 matches out of 5790 PubChem CIDs from CMap:
# table(na.omit(unique(cmap_pubchem_cid)) %in%
#           na.omit(descriptors2D$PubChem_CID))

# Intersection by compound name ------------------------------------------------
cmap_name <- unique(na.omit(compoundInfo$pert_iname))
descriptors2D_name <- unique(na.omit(descriptors2D$Drug_name))
table(cmap_name %in% descriptors2D_name)
# 61 matches out of 6127 CMap compound names

stripNonAlphaNumericChr <- function(str) gsub("[^[:alnum:] ]", "", str)
table(stripNonAlphaNumericChr(cmap_name) %in%
          stripNonAlphaNumericChr(descriptors2D_name))
# 68 matches out of 6059 CMap compound names

cmap_name_2 <- unique(na.omit(attr(perts, "metadata")$pert_iname))
table(cmap_name_2 %in% descriptors2D_name)
# 62 matches out of 28926 CMap compound names

table(stripNonAlphaNumericChr(cmap_name_2) %in%
          stripNonAlphaNumericChr(descriptors2D_name))
# 63 matches out of 28926 CMap compound names

drugSets <- prepareDrugSets(descriptors2D, "PubChem_CID")
analyseDrugSetEnrichment(rankedPerts, drugSets)

# 3D molecular descriptors -----------------------------------------------------
# descriptors3D <- data.table::fread("compound_descriptors_3D.csv",
#                                    na.strings=c("", "NA"))
# analyseDrugSetEnrichment(rankedPerts, drugSets)
