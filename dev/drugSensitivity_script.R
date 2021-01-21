# Prepare expression and drug sensitivity association data ---------------------
# ctrp  <- prepareExpressionDrugSensitivityAssociation("CTRP")
# gdsc  <- prepareExpressionDrugSensitivityAssociation("GDSC")
# nci60 <- prepareExpressionDrugSensitivityAssociation("NCI60")

# Save (and manually upload) objects as RDS/h5 files ---------------------------
saveRDS(ctrp,  "expressionDrugSensitivityCorCTRP2.1.rds")
saveRDS(gdsc,  "expressionDrugSensitivityCorGDSC7.rds")
# saveRDS(nci60, "expressionDrugSensitivityCorNCI60.rds")

# Download and load expression and drug sensitivity association data -----------
ctrp  <- loadExpressionDrugSensitivityAssociation("CTRP")
gdsc  <- loadExpressionDrugSensitivityAssociation("GDSC")
nci60 <- loadExpressionDrugSensitivityAssociation("NCI60")

# Predict targeting drugs ------------------------------------------------------
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
