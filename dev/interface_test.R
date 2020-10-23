# Action verbs: launch, start, open, run

# launchDiffExprLoader() -------------------------------------------------------
diffExpr <- launchDiffExprLoader()

# launchCMapDataLoader() -------------------------------------------------------
cmapKD <- launchCMapDataLoader(
    cellLine="HepG2",
    perturbationType="Consensus signature from shRNAs targeting the same gene")
cmapCompounds <- launchCMapDataLoader(
    cellLine="HepG2",
    perturbationType="Compound")

launchMetadataViewer(cmapPerturbationsKD)

# Rank similar perturbations ---------------------------------------------------
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
predicted   <- predictTargetingDrugs(diffExprStat, assoc)
launchResultPlotter(predicted)

# Plot targeting drugs vs similar perturbations
launchResultPlotter(predicted, compareCompounds)

# Drug set enrichment analysis -------------------------------------------------
descriptors <- loadDrugDescriptors("NCI60", "2D")
drugSets    <- prepareDrugSets(descriptors)
dsea        <- analyseDrugSetEnrichment(drugSets, predicted)

# launchDrugSetEnrichmentAnalysis(drugSets, predicted) -------------------------
plotDrugSetEnrichment(drugSets, predicted, selectedSets=head(dsea$pathway, 6))
