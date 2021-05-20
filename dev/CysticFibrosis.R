# SRP219679: Cystic Fibrosis in whole blood
folder <- "../cTRAP_projects/CysticFibrosis/SRP219679/ReadsPerGene/"
metadata <- loadLocalFiles(file.path(folder, ".."), ignore=".gz")
metadata <- metadata[[1]][["Sample metadata"]]

# Prepare gene expression ------------------------------------------------------
library(psichomics)
countFiles <- list.files(folder, pattern="ReadsPerGene", full.names=TRUE)
geneExpr <- as.data.frame(prepareGeneQuant(countFiles, output=NULL))
colnames(geneExpr) <- gsub("_", "", colnames(geneExpr))

library(org.Hs.eg.db)
rownames(geneExpr) <- convertGeneIdentifiers(org.Hs.eg.db, geneExpr[[1]])
geneExpr <- geneExpr[ , -1]

summary(sort(unlist(geneExpr["CFTR", ])))
# No more than 16 reads for CFTR... no need to continue down this path...

# Normalise gene expression ----------------------------------------------------
plotGeneExprPerSample(geneExpr)
plotDistribution(log10(colSums(geneExpr)))
plotDistribution(log10(colSums(geneExpr)), metadata$source_name)
# outliers <- c("SRR10037639", "SRR10037672", "SRR10037643")

filter <- filterGeneExpr(geneExpr, minCounts=10, minTotalCounts=15)
table(filter)

geneExprNorm <- normaliseGeneExpression(geneExpr, geneFilter=filter,
                                        method="TMM", log2transform=TRUE)

plotGeneExprPerSample(geneExprNorm)
plotDistribution(log10(colSums(geneExprNorm)))
plotDistribution(log10(colSums(geneExprNorm)), metadata$source_name)
# Maybe remove SRR10037664?

# Create groups ----------------------------------------------------------------
groups <- NULL
groups$type   <- createGroupByAttribute("source_name", metadata)
groups$gender <- createGroupByAttribute("gender", metadata)
groups$age    <- createGroupByAttribute("Age", metadata)

# PCA --------------------------------------------------------------------------
pca <- performPCA(t(geneExprNorm))
psichomics::plotPCA(pca, groups=groups$type)

# Differentially expressed genes -----------------------------------------------
design <- model.matrix(~ source_name + gender + Age, metadata)
library(limma)
fit <- lmFit(geneExprNorm, design)
fit <- eBayes(fit)

################################################################################
################################################################################
################################################################################

# SRP273075
folder <- "../cTRAP_projects/CysticFibrosis/SRP273075/ReadsPerGene/"
metadata <- data.table::fread(file.path(folder, "..", "SraRunTable.txt"),
                              data.table=FALSE)
rownames(metadata) <- metadata$Run

# Prepare gene expression ------------------------------------------------------
library(psichomics)
countFiles <- list.files(folder, pattern="ReadsPerGene", full.names=TRUE)
geneExpr <- as.data.frame(prepareGeneQuant(countFiles, output=NULL))
colnames(geneExpr) <- gsub("_", "", colnames(geneExpr))

library(org.Hs.eg.db)
rownames(geneExpr) <- convertGeneIdentifiers(org.Hs.eg.db, geneExpr[[1]])
geneExpr <- geneExpr[ , -1]

# Check if there are enough CFTR reads to continue analysis
summary(sort(unlist(geneExpr["CFTR", ])))

# Normalise gene expression ----------------------------------------------------
plotGeneExprPerSample(geneExpr)
plotDistribution(log10(colSums(geneExpr)))
plotDistribution(log10(colSums(geneExpr)), metadata$Genotype)
plotDistribution(log10(colSums(geneExpr)), metadata$growth_medium)
# outliers <- c("SRR12287872", "SRR12287870", "SRR12287871", "SRR12287869")

filter <- filterGeneExpr(geneExpr, minCounts=10, minTotalCounts=15)
table(filter)

geneExprNorm <- normaliseGeneExpression(geneExpr, geneFilter=filter,
                                        method="TMM", log2transform=TRUE)

plotGeneExprPerSample(geneExprNorm)
plotDistribution(log10(colSums(geneExprNorm)))
plotDistribution(log10(colSums(geneExprNorm)), metadata$source_name)
# Maybe remove SRR10037664?

# Create groups ----------------------------------------------------------------
groups <- NULL
groups$genotype <- createGroupByAttribute("Genotype", metadata)
groups$medium   <- createGroupByAttribute("growth_medium", metadata)

# PCA --------------------------------------------------------------------------
pca <- performPCA(t(geneExprNorm))
psichomics::plotPCA(pca, groups=groups$medium)
psichomics::plotPCA(pca, groups=groups$genotype)

# Differentially expressed genes -----------------------------------------------
design <- model.matrix(~ growth_medium + Genotype, metadata)
library(limma)
fit <- lmFit(geneExprNorm, design)
fit <- eBayes(fit)

table <- topTable(fit, number=nrow(fit), coef=3)
table$genes <- rownames(table)

logFCthreshold  <- abs(table$logFC) > 1
pvalueThreshold <- table$adj.P.Val < 0.05

library(ggplot2)
library(ggrepel)
ggplot(table, aes(logFC, -log10(adj.P.Val))) +
    geom_point(data=table[logFCthreshold & pvalueThreshold, ],
               colour="orange", alpha=0.5, size=3) +
    geom_point(data=table[!logFCthreshold | !pvalueThreshold, ],
               colour="gray", alpha=0.5, size=3) +
    geom_point(data=table["CFTR", ], colour="red") +
    geom_text_repel(data=table["CFTR", ],
                    aes(label=genes), box.padding=0.4, size=5) +
    theme_light(16) +
    ylab("-log10(q-value)")

stat <- data.frame("Genes"=table$genes, "logFC"=table$logFC,
                   stringsAsFactors=FALSE)
saveRDS(stat, file.path(folder, "..", "SRP273075_DGE.rds"))
