# genomicsDrugSensitivityCorGDSC <- prepareGenomicsDrugSensitivityAssociation("GDSC")
# genomicsDrugSensitivityCorCTRP <- prepareGenomicsDrugSensitivityAssociation("CTRP")
# genomicsDrugSensitivityCorNCI60 <- prepareGenomicsDrugSensitivityAssociation("NCI60")

ctrp  <- prepareExpressionDrugSensitivityAssociation("CTRP")
gdsc  <- prepareExpressionDrugSensitivityAssociation("GDSC")
nci60 <- prepareExpressionDrugSensitivityAssociation("NCI60")

corCTRP  <- predictTargetingDrugs(diffExprStat, ctrp)
corGDSC  <- predictTargetingDrugs(diffExprStat, gdsc)
corNCI60 <- predictTargetingDrugs(diffExprStat, nci60)

# Plot CMap perturbation similarity with drug sensitivity results --------------
getCMapResults <- function(diffExprStat, ..., subset=FALSE) {
    metadata <- loadCMapData("cmapMetadata.txt")
    metadataCompound <- filterCMapMetadata(metadata, ...)

    if (subset) {
        metadataSplit <- split(metadataCompound, metadataCompound$cell_id)
        if (length(metadataSplit) < 2) {
            sample <- rep(TRUE, nrow(metadataSplit[[1]]))[1:100]
        } else {
            sample <- parseCMapID(metadataCompound$sig_id) %in%
                intersect(parseCMapID(metadataSplit[[1]]$sig_id),
                          parseCMapID(metadataSplit[[2]]$sig_id))[1:100]
        }
        metadataCompound <- metadataCompound[sample, ]
    }

    perts <- prepareCMapPerturbations(metadataCompound, "cmapZscores.gctx",
                                      "cmapGeneInfo.txt",
                                      "cmapCompoundInfo_drugs.txt")
    cmp <- rankSimilarPerturbations(diffExprStat, perts)
    return(cmp)
}
cmapResults <- getCMapResults(diffExprStat, perturbationType="Compound")

plot_targetingDrugsVScmapResults(corCTRP,  cmapResults2, "spearman_rank",
                                 labelColumn="target")
plot_targetingDrugsVScmapResults(corGDSC,  cmapResults2, "spearman_rank")
plot_targetingDrugsVScmapResults(corNCI60, cmapResults2, "spearman_rank")

# Drug set enrichment analysis -------------------------------------------------

# Create drug sets from 2D and 3D descriptors
descriptors2D <- data.table::fread("compound_descriptors_NCI60_2D.txt",
                                   na.strings=c("", "NA"))
descriptors2D$V45 <- NULL
descriptors2D <- descriptors2D[seq(nrow(descriptors2D) - 1)]
drugSets <- prepareDrugSets(descriptors2D, "PubChem_SID")

descriptors3D <- data.table::fread("compound_descriptors_3D.csv",
                                   na.strings=c("", "NA"))
# analyseDrugSetEnrichment(cmapResults, drugSets)
