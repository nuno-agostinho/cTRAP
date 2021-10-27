#!/usr/bin/env Rscript

# Nuno Agostinho, 2 December 2020

lapTime <- function(msg) message(Sys.time(), " - ", msg)
lapTime("script started")

args <- commandArgs(trailingOnly=TRUE)
file <- args[1]

# Predict targeting drugs
diffExpr <- readRDS(file)

library(cTRAP)
nci60 <- listExpressionDrugSensitivityAssociation()[[3]]

lapTime("loading NCI60 expressiong and drug sensitivity association")
assoc <- loadExpressionDrugSensitivityAssociation(
    nci60, "input/expressionDrugSensitivityCorNCI60.h5")

lapTime("predicting targeting drugs")
predicted <- predictTargetingDrugs(diffExpr, assoc)

out <- sprintf("output/predictedTargetingDrugs_%s.rds",
               gsub(" ", "_", Sys.time()))
saveRDS(predicted, out)

# Molecular descriptor enrichment analysis
lapTime("loading CMap 3D drug descriptors")
descriptors <- loadDrugDescriptors(
    "NCI60", "3D", "input/molecular_descriptors_NCI60_3D.rds")

lapTime("preparing drug sets")
drugSets <- prepareDrugSets(descriptors)

lapTime("analysing drug set enrichment")
dsea <- analyseDrugSetEnrichment(drugSets, predicted)

lapTime("script ended")
