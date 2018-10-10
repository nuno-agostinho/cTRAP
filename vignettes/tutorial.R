## ----load package--------------------------------------------------------
library(cTRAP)

## ----load ENCODE samples, eval=FALSE-------------------------------------
#  gene <- "EIF4G1"
#  cellLine <- "HepG2"
#  
#  ENCODEmetadata <- downloadENCODEknockdownMetadata(cellLine, gene)
#  table(ENCODEmetadata$`Experiment target`)
#  length(unique(ENCODEmetadata$`Experiment target`))
#  
#  ENCODEsamples <- downloadENCODEsamples(ENCODEmetadata)
#  counts <- prepareENCODEgeneExpression(ENCODEsamples)

## ----load data, include=FALSE--------------------------------------------
gene <- "EIF4G1"
cellLine <- "HepG2"

data("ENCODEmetadata")
data("counts")

## ----differential gene expression----------------------------------------
# Remove low coverage (at least 10 counts shared across two samples)
minReads   <- 10
minSamples <- 2
filter <- rowSums(counts[ , -c(1, 2)] >= minReads) >= minSamples
counts <- counts[filter, ]

# Perform differential gene expression analysis
diffExpr <- performDifferentialExpression(counts)

## ----t-statistics--------------------------------------------------------
# Get t-statistics of differential expression with respective gene names
diffExprStat <- diffExpr$t
names(diffExprStat) <- diffExpr$Gene_symbol

## ----L1000 metadata conditions-------------------------------------------
# Check available conditions for L1000 perturbations
l1000metadata <- downloadL1000data("l1000metadata.txt", "metadata")
getL1000conditions(l1000metadata)

## ----L1000 knockdown perturbations, eval=FALSE---------------------------
#  # Code for loading CMap gene KD HepG2 data
#  l1000metadataKnockdown <- filterL1000metadata(
#      l1000metadata, cellLine="HepG2",
#      perturbationType="Consensus signature from shRNAs targeting the same gene")
#  l1000zscores  <- downloadL1000data("l1000zscores.gctx", "zscores",
#                                     l1000metadataKnockdown$sig_id)
#  l1000geneInfo <- downloadL1000data("l1000geneInfo.txt", "geneInfo")
#  
#  l1000perturbationsKnockdown <- loadL1000perturbations(
#      l1000metadataKnockdown, l1000zscores, l1000geneInfo)

## ----L1000 small molecule perturbations, eval=FALSE----------------------
#  l1000metadataSmallMolecules <- filterL1000metadata(
#      l1000metadata, cellLine="HepG2", timepoint="24 h",
#      dosage="5 \U00B5M", # \U00B5 is the unicode code for the micro symbol
#      perturbationType="Compound")
#  l1000zscores  <- downloadL1000data("l1000zscores.gctx", "zscores",
#                                     l1000metadataSmallMolecules$sig_id)
#  l1000geneInfo <- downloadL1000data("l1000geneInfo.txt")
#  
#  l1000perturbationsSmallMolecules <- loadL1000perturbations(
#      l1000metadataSmallMolecules, l1000zscores, l1000geneInfo,
#      sanitizeCompoundNames=TRUE)

## ----load L1000 perturbations, include=FALSE-----------------------------
data("l1000perturbationsKnockdown")
data("l1000perturbationsSmallMolecules")

## ----compare knockdown, include=FALSE------------------------------------
compareKnockdown <- list()

# Compare against L1000 using Spearman correlation
compareKnockdown$spearman <- compareAgainstL1000(
    diffExprStat, l1000perturbationsKnockdown, cellLine, method="spearman")

# Compare against L1000 using Pearson correlation
compareKnockdown$pearson <- compareAgainstL1000(
    diffExprStat, l1000perturbationsKnockdown, cellLine, method="pearson")

# Compare against L1000 using gene set enrichment analysis (GSEA) with the top
# and bottom 150 genes
compareKnockdown$gsea <- compareAgainstL1000(
    diffExprStat, l1000perturbationsKnockdown, cellLine, method="gsea",
    geneSize=150)

## ----compare small molecules---------------------------------------------
compareSmallMolecule <- list()
# Compare against L1000 using Spearman correlation
compareSmallMolecule$spearman <- compareAgainstL1000(
    diffExprStat, l1000perturbationsSmallMolecules, cellLine, method="spearman")

# Compare against L1000 using Pearson correlation
compareSmallMolecule$pearson <- compareAgainstL1000(
    diffExprStat, l1000perturbationsSmallMolecules, cellLine, method="pearson")

# Compare against L1000 using gene set enrichment analysis (GSEA) with the top
# and bottom 150 genes
compareSmallMolecule$gsea <- compareAgainstL1000(
    diffExprStat, l1000perturbationsSmallMolecules, cellLine, method="gsea",
    geneSize=150)

## ----order knockdown perturbation comparison-----------------------------
# Order knockdown perturbations according to similarity
compareKnockdown$spearman_ordered <- compareKnockdown$spearman[
    order(compareKnockdown$spearman$HepG2_t_spearman_coef, decreasing=TRUE)]
compareKnockdown$pearson_ordered <- compareKnockdown$pearson[
    order(compareKnockdown$pearson$HepG2_t_pearson_coef, decreasing=TRUE)]
compareKnockdown$gsea_ordered <- compareKnockdown$gsea[
    order(compareKnockdown$gsea$HepG2_WTCS, decreasing=FALSE)]

# Most positively associated perturbations (note that EIF4G1 knockdown is the
# 6th, 1st and 2nd most positively associated perturbation based on Spearman
# correlation, Pearson correlation and GSEA, respectively)
head(compareKnockdown$spearman_ordered)
head(compareKnockdown$pearson_ordered)
head(compareKnockdown$gsea_ordered)

# Most negatively associated perturbations
tail(compareKnockdown$spearman_ordered)
tail(compareKnockdown$pearson_ordered)
tail(compareKnockdown$gsea_ordered)

## ----order small molecule perturbation comparison------------------------
# Order small molecule perturbations according to similarity
compareSmallMolecule$spearman_ordered <- compareSmallMolecule$spearman[
    order(compareSmallMolecule$spearman$HepG2_t_spearman_coef, decreasing=TRUE)]
compareSmallMolecule$pearson_ordered <- compareSmallMolecule$pearson[
    order(compareSmallMolecule$pearson$HepG2_t_pearson_coef, decreasing=TRUE)]
compareSmallMolecule$gsea_ordered <- compareSmallMolecule$gsea[
    order(compareSmallMolecule$gsea$HepG2_WTCS, decreasing=FALSE)]

# Most positively associated perturbations
head(compareSmallMolecule$spearman_ordered)
head(compareSmallMolecule$pearson_ordered)
head(compareSmallMolecule$gsea_ordered)

# Most negatively associated perturbations
tail(compareSmallMolecule$spearman_ordered)
tail(compareSmallMolecule$pearson_ordered)
tail(compareSmallMolecule$gsea_ordered)

## ----plot knockdown perturbation comparison------------------------------
# Plot relationship with EIF4G1 knockdown
EIF4G1knockdown <- grep("EIF4G1", compareKnockdown$gsea_ordered$genes,
                        value=TRUE)
plotL1000comparison(compareKnockdown$spearman, EIF4G1knockdown)
plotL1000comparison(compareKnockdown$pearson, EIF4G1knockdown)
plotL1000comparison(compareKnockdown$gsea, EIF4G1knockdown, topGenes=TRUE)
plotL1000comparison(compareKnockdown$gsea, EIF4G1knockdown, topGenes=FALSE)

# Plot relationship with most negatively associated perturbation
plotL1000comparison(compareKnockdown$spearman,
                    tail(compareKnockdown$spearman_ordered, 1)$genes)
plotL1000comparison(compareKnockdown$pearson,
                    tail(compareKnockdown$pearson_ordered, 1)$genes)
plotL1000comparison(compareKnockdown$gsea,
                    tail(compareKnockdown$gsea_ordered, 1)$genes, 
                    topGenes=TRUE)
plotL1000comparison(compareKnockdown$gsea, 
                    tail(compareKnockdown$gsea_ordered, 1)$genes, 
                    topGenes=FALSE)

## ----plot small molecule perturbation comparison-------------------------
# Plot relationship with most positively associated perturbation
plotL1000comparison(compareSmallMolecule$spearman,
                    head(compareSmallMolecule$spearman_ordered, 1)$genes)
plotL1000comparison(compareSmallMolecule$pearson,
                    head(compareSmallMolecule$pearson_ordered, 1)$genes)
plotL1000comparison(compareSmallMolecule$gsea,
                    head(compareSmallMolecule$gsea_ordered, 1)$genes)

# Plot relationship with most negatively associated perturbation
plotL1000comparison(compareSmallMolecule$spearman,
                    tail(compareSmallMolecule$spearman_ordered, 1)$genes)
plotL1000comparison(compareSmallMolecule$pearson,
                    tail(compareSmallMolecule$pearson_ordered, 1)$genes)
plotL1000comparison(compareSmallMolecule$gsea,
                    tail(compareSmallMolecule$gsea_ordered, 1)$genes)

