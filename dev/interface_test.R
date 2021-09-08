# Action verbs: launch, start, open, run
# - launchDiffExprLoader()
# - launchCMapDataLoader()
# - launchMetadataViewer()
# - launchResultPlotter()

diffExpr <- launchDiffExprLoader(cellLine="HepG2", gene="EIF4G1")

# Rank similar perturbations ---------------------------------------------------
cmapKD <- launchCMapDataLoader(
    cellLine="HepG2",
    perturbationType="Consensus signature from shRNAs targeting the same gene")
cmapCompounds <- launchCMapDataLoader(
    cellLine="HepG2",
    perturbationType="Compound")
cmapPerts <- launchCMapDataLoader(cellLine="HepG2")

launchMetadataViewer(cmapKD, cmapCompounds, cmapPerts)

compareKD        <- rankSimilarPerturbations(diffExpr, cmapKD)
compareCompounds <- rankSimilarPerturbations(diffExpr, cmapCompounds)
comparePerts     <- rankSimilarPerturbations(diffExpr, cmapPerts)

launchResultPlotter(compareCompounds, compareKD, comparePerts)

# Predict targeting drugs ------------------------------------------------------
listExpressionDrugSensitivityAssociation()
assocMatrix <- listExpressionDrugSensitivityAssociation()[[1]]
assoc       <- loadExpressionDrugSensitivityAssociation(assocMatrix)
predicted   <- predictTargetingDrugs(diffExpr, assoc)
launchResultPlotter(predicted)

# Plot targeting drugs vs similar perturbations
launchResultPlotter(predicted, compareCompounds)

# Analyse drug set enrichment --------------------------------------------------
descriptors <- loadDrugDescriptors("NCI60", "3D")
drugSets    <- prepareDrugSets(descriptors)
launchDrugSetEnrichmentAnalyser(drugSets, predicted)
