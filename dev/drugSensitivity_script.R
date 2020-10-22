# genomicsDrugSensitivityCorGDSC <- prepareGenomicsDrugSensitivityAssociation("GDSC")
# genomicsDrugSensitivityCorCTRP <- prepareGenomicsDrugSensitivityAssociation("CTRP")
# genomicsDrugSensitivityCorNCI60 <- prepareGenomicsDrugSensitivityAssociation("NCI60")

ctrp  <- loadGenomicsDrugSensitivityAssociation("CTRP")
gdsc  <- loadGenomicsDrugSensitivityAssociation("GDSC")
nci60 <- loadGenomicsDrugSensitivityAssociation("NCI60")

corCTRP  <- predictTargetingDrug(diffExprStat, ctrp)
corGDSC  <- predictTargetingDrug(diffExprStat, gdsc)
corNCI60 <- predictTargetingDrug(diffExprStat, nci60)

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

# analyseDrugSetEnrichment <- function(stats, sets) {
#     if (is(stats, "similarPerturbations")) {
#         # Use negative values to invert perturbation relevance
#         rank   <- -cmp[["rankProduct_rank"]]
#         values <- na.omit(setNames(rank, cmp[[1]])) # Discard missing values
#
#         metadata     <- attr(stats, "metadata")
#         compoundInfo <- attr(stats, "compoundInfo")
#
#         # Convert to perturbation name (discard if non-matching)
#         names(values) <- metadata[
#             metadata$sig_id == names(values), ][["pert_iname"]]
#         values <- values[!is.na(names(values))]
#
#         # Convert to PubChem identifier
#         names(values) <- compoundInfo[
#             match(names(values), compoundInfo$pert_iname), ][["pubchem_cid"]]
#
#         values <- values[!is.na(names(values))]
#         dt <- data.table(PubChem=names(values), rank=values)[
#             , .("rank"=max(rank)), by=PubChem]
#         stats <- setNames(dt$rank, dt$PubChem)
#     }
#     stats <- sort(stats)
#     browser()
#     fgsea(sets, stats, nperm=10000, maxSize=500)
# }
# analyseDrugSetEnrichment(cmapResults, drugSets)
