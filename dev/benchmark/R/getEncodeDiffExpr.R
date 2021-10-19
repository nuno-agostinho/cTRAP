#!/usr/bin/env Rscript

gene <- "EIF4G1"
cellLine <- "HepG2"

library(cTRAP)
ENCODEmetadata <- downloadENCODEknockdownMetadata(cellLine, gene)
table(ENCODEmetadata$`Experiment target`)
length(unique(ENCODEmetadata$`Experiment target`))

ENCODEsamples <- loadENCODEsamples(ENCODEmetadata)[[1]]
counts <- prepareENCODEgeneExpression(ENCODEsamples)

# Remove low coverage (at least 10 counts shared across two samples)
minReads   <- 10
minSamples <- 2
filter <- rowSums(counts[ , -c(1, 2)] >= minReads) >= minSamples
counts <- counts[filter, ]

# Convert ENSEMBL identifier to gene symbol
counts$gene_id <- convertENSEMBLtoGeneSymbols(counts$gene_id)

# Perform differential gene expression analysis
diffExpr <- performDifferentialExpression(counts)

diffExprStat <- diffExpr$t
names(diffExprStat) <- diffExpr$Gene_symbol
saveRDS(diffExprStat, "../input/diffExprStat.rds")
