# Action verbs: launch, start, open, run
# - launchDiffExprLoader()
# - launchCMapDataLoader()
# - launchMetadataViewer()
# - launchResultPlotter()

diffExpr <- launchDiffExprLoader()

# Rank similar perturbations ---------------------------------------------------
cmapKD <- launchCMapDataLoader(
    cellLine="HepG2",
    perturbationType="Consensus signature from shRNAs targeting the same gene")
cmapCompounds <- launchCMapDataLoader(
    cellLine="HepG2",
    perturbationType="Compound")

launchMetadataViewer(cmapKD)

compareKD        <- rankSimilarPerturbations(diffExpr, cmapKD)
compareCompounds <- rankSimilarPerturbations(diffExpr, cmapCompounds)

launchResultPlotter(compareKD)
launchResultPlotter(compareCompounds)
launchResultPlotter(compareCompounds, compareKD)

launchMetadataViewer(compareKD)

# Predict targeting drugs ------------------------------------------------------
listExpressionDrugSensitivityAssociation()
assocMatrix <- listExpressionDrugSensitivityAssociation()[[1]]
assoc       <- loadExpressionDrugSensitivityAssociation(assocMatrix)
predicted   <- predictTargetingDrugs(diffExpr, assoc)
launchResultPlotter(predicted)

# Plot targeting drugs vs similar perturbations
launchResultPlotter(predicted, compareCompounds)

# Analyse drug set enrichment --------------------------------------------------
descriptors <- loadDrugDescriptors("NCI60", "2D")
drugSets    <- prepareDrugSets(descriptors)
launchDrugSetEnrichmentAnalyser(drugSets, predicted)
