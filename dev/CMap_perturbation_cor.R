# Nuno Agostinho (NMorais Lab), Instituto de Medicina Molecular, 27 May 2019
# Correlation matrix of CMap perturbations

# Load CMap perturbations
transFile <- "output/transposed.rds"
if (file.exists(transFile)) {
    message("Loading transposed CMap perturbation matrix...")
    trans <- readRDS(transFile)
} else {
    library(cTRAP)
    perts <- loadCMapPerturbations("cmapMetadata.txt", "cmapZscores.gctx",
                                   "cmapGeneInfo.txt", "cmapCompoundInfo.txt")

    # Transpose data to have genes as columns (required for correlation matrix)
    message("Transposing CMap perturbation matrix...")
    library(WGCNA) # Other libraries that we tried: ccaPP, fastCor, ppcor
    time <- Sys.time()
    trans <- transposeBigData(perts)
    rm(perts)
    gc()
    saveRDS(trans, transFile)
    print(Sys.time() - time)
}

# Creating Very Large Correlation/Covariance Matrices via the base cor function
#
# bigcor uses the framework of the 'ff' package to store the correlation matrix
# in a file. The complete matrix is created by filling a large preallocated
# empty matrix with sub-matrices at the corresponding positions.
#
# Run time: 20s for a 10000 x 100 matrix
# Memory: 18.63 GiB for a 50000 x 100 matrix (50000 ^ 2 x 8 Byte = 18.63 GiB)
message("Performing correlation...")
library("propagate")

time <- Sys.time()
cor_ff <- bigcor(trans, method="spearman", size=2000, verbose=TRUE)
colnames(cor_ff) <- rownames(cor_ff) <- colnames(trans)
print(Sys.time() - time)

message("Saving file as ff...")
library(ff)
time <- Sys.time()
ffsave(cor_ff,
       file=sprintf("output/correlation_ff_%s", format(Sys.time(), "%Y-%m-%d")))
print(Sys.time() - time)

# Save as matrix
message("Converting to matrix...")
time <- Sys.time()
corMat_mat <- cor_ff[1:nrow(cor_ff), 1:ncol(cor_ff)]
colnames(corMat_mat) <- rownames(corMat_mat) <- colnames(trans)
print(Sys.time() - time)

message("Saving file as matrix...")
time <- Sys.time()

saveRDS(corMat_mat,
        sprintf("output/correlation_%s.rds", format(Sys.time(), "%Y-%m-%d")))
print(Sys.time() - time)

# Correlation matrices were performed on:
#     1) CMap compound perturbations only
#     2) CMap knockdown perturbations only
#     3) CMap overexpression perturbations only
#     4) all CMap perturbations
