#!/usr/bin/env Rscript

options(keep.source=TRUE)
message(sprintf("Running CMap perturbation ranking at %s...",
                format(Sys.time())))
time <- Sys.time()

# Process arguments ------------------------------------------------------------
args <- commandArgs()
file <- args[length(args) - 2]
basename <- tools::file_path_sans_ext(file)

perturbationTypeStr <- args[length(args) - 1]
perturbationType <- switch(
    perturbationTypeStr,
    "all"=NULL,
    "compound"="Compound",
    "knockdown"="Consensus signature from shRNAs targeting the same gene",
    "overexpression"="cDNA for overexpression of wild-type gene")

loadZscores <- args[length(args)]
loadZscores <- loadZscores == "TRUE"

headerMessage <- function(text) {
    width <- nchar(text) + 1
    message("\n", text, " ", rep("=", max(80 - width, 4)))
}

library(profvis)
prof <- profvis({
# Load processed data ----------------------------------------------------------
headerMessage("Load processed data")
message(sprintf("Loading input data from file %s...", file))
diffExpr <- readRDS(file)

message(sprintf("\nData for first genes out of %s in total:", length(diffExpr)))
print(head(diffExpr))

message("\nSummary statistics:")
print(summary(diffExpr))

# Load CMap compound perturbtions ----------------------------------------------
library(cTRAP)
out <- sprintf("output/out_diffExpr_vs_%s_loadZscores-%s.rds",
               perturbationTypeStr, loadZscores)
headerMessage("Rank similar CMap perturbations")
metadata <- filterCMapMetadata("input/cmapMetadata.txt",
                               perturbationType=perturbationType)
compoundPerts <- prepareCMapPerturbations(
    metadata,
    "input/cmapZscores.gctx",
    "input/cmapGeneInfo.txt",
    "input/cmapCompoundInfo.txt",
    loadZscores=loadZscores)

cmp <- rankSimilarPerturbations(diffExpr, compoundPerts)
message(sprintf("Saving results to %s...", out))
saveRDS(cmp, out)
})

htmlwidgets::saveWidget(prof, sprintf('profvis_%s.html', out))
message(sprintf("Script finished successfully in %s",
                format(Sys.time() - time)))
