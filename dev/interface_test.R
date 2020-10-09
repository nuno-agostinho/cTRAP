# Action verbs: launch, start, open, run

# launchDiffExprLoader() -------------------------------------------------------
diffExpr <- launchDiffExprLoader()

# launchCMapDataLoader() -------------------------------------------------------
cmapPerturbationsKD <- launchCMapDataLoader(
    cellLine="HepG2",
    perturbationType="Consensus signature from shRNAs targeting the same gene")
cmapPerturbationsCompounds <- launchCMapDataLoader(
    cellLine="HepG2", perturbationType="Compound")

launchMetadataViewer(cmapPerturbationsKD)

# Rank similar perturbations ---------------------------------------------------
compareKD        <- rankSimilarPerturbations(diffExprStat, cmapPerturbationsKD)
compareCompounds <- rankSimilarPerturbations(diffExprStat,
                                             cmapPerturbationsCompounds)

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
descriptors <- loadDrugDescriptors("CMap", "2D")
drugSets    <- prepareDrugSets(descriptors)
dsea        <- analyseDrugSetEnrichment(drugSets, predicted)

# launchDrugSetEnrichmentAnalysis(drugSets, predicted) -------------------------
plotDrugSetEnrichment(drugSets, predicted, selectedSets=head(dsea$pathway, 6))
