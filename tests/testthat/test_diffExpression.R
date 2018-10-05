context("Differential gene expression")

test_that("Perform differential gene expression", {
    gene <- "EIF4G1"
    cellLine <- "HepG2"

    ENCODEmetadata <- downloadENCODEknockdownMetadata(cellLine, gene)
    counts <- loadENCODEgeneExpression(ENCODEmetadata)

    # Remove low coverage (at least 10 counts shared across two samples)
    minReads   <- 10
    minSamples <- 2
    filter <- rowSums(counts[ , 3:6] >= minReads) >= minSamples
    counts <- counts[filter, ]

    # Perform differential gene expression analysis
    diffExpr <- performDifferentialExpression(counts)
    expect_is(diffExpr, "data.frame")
    expect_identical(colnames(diffExpr), c("Gene_symbol", "logFC", "AveExpr",
                                           "t", "P.Value", "adj.P.Val", "B"))
})
