#!/usr/bin/env Rscript

lapTime <- function(msg) message(Sys.time(), " - ", msg)
lapTime("script started")

# Process arguments ------------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
file                <- args[1]
perturbationTypeStr <- args[2]
loadZscores         <- args[3]
threads             <- args[4]

basename <- tools::file_path_sans_ext(file)
perturbationType <- switch(
    perturbationTypeStr,
    "all"=NULL,
    "compound"="Compound",
    "knockdown"="Consensus signature from shRNAs targeting the same gene",
    "overexpression"="cDNA for overexpression of wild-type gene")
loadZscores <- loadZscores == "TRUE"
threads     <- as.numeric(threads)

headerMessage <- function(text) {
    width <- nchar(text) + 1
    message("\n", text, " ", rep("=", max(80 - width, 4)))
}

headerMessage("Script arguments")
message("- File: ", file)
message("- Perturbation type: ", perturbationType)
message("- Load Z scores: ", loadZscores)
message("- Threads: ", threads)
lapTime("arguments processed")

# Load processed data ----------------------------------------------------------
headerMessage("Load processed data")
message(sprintf("Loading input data from file %s...", file))
diffExpr <- readRDS(file)

message(sprintf("\nData for first genes out of %s in total:", length(diffExpr)))
print(head(diffExpr))

message("\nSummary statistics:")
print(summary(diffExpr))

lapTime("user-provided data loaded")

# Load CMap compound perturbtions ----------------------------------------------
library(cTRAP)
headerMessage("Rank similar CMap perturbations")

message("Filtering CMap metadata...")
out <- sprintf("output/cTRAP_%s_out_diffExpr_vs_%s_loadZscores-%s.rds",
               packageVersion("cTRAP"), perturbationTypeStr, loadZscores)
metadata <- filterCMapMetadata("input/cmapMetadata.txt",
                               perturbationType=perturbationType)
lapTime("CMap metadata filtered")

message("Preparing CMap perturbations...")
compoundPerts <- prepareCMapPerturbations(
    metadata,
    "input/cmapZscores.gctx",
    "input/cmapGeneInfo.txt",
    "input/cmapCompoundInfo.txt",
    loadZscores=loadZscores)
lapTime("CMap perturbations prepared")

message("Ranking CMap perturbations...")
cmp <- rankSimilarPerturbations(diffExpr, compoundPerts, threads=threads)
lapTime("CMap perturbations ranked")

message(sprintf("Saving results to %s...", out))
saveRDS(cmp, out)

message(sprintf("Script finished successfully in %s",
                format(Sys.time() - time)))
