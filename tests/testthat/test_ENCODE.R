context("Loading ENCODE data")

cellLine <- "HepG2"
gene     <- "EIF4G1"

ENCODEmetadata <- downloadENCODEknockdownMetadata(cellLine, gene)
counts <- loadENCODEgeneExpression(ENCODEmetadata)

# Remove low coverage (at least 10 counts shared across two samples)
minReads   <- 10
minSamples <- 2
filter <- rowSums(counts[ , 3:6] >= minReads) >= minSamples
counts <- counts[filter, ]

# Perform differential gene expression analysis
diffExpr <- performDifferentialExpression(counts)

test_that("Load ENCODE gene expression for EIF4G1", {
    expect_is(ENCODEmetadata, "tbl_df")
    expect_identical(ENCODEmetadata$`Experiment target`, rep("EIF4G1", 2))

    expect_identical(colnames(counts),
                     c("gene_id", "transcript_id(s)",
                       paste0("shRNA", 1:2), paste0("control", 1:2)))
})

test_that("Perform differential gene expression", {
    expect_is(diffExpr, "data.frame")
    expect_identical(colnames(diffExpr), c("Gene_symbol", "logFC", "AveExpr",
                                           "t", "P.Value", "adj.P.Val", "B"))
})
