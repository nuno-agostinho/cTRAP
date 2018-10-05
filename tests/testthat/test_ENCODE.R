context("Loading ENCODE data")

test_that("Load ENCODE gene expression for EIF4G1", {
    cellLine <- "HepG2"
    gene     <- "EIF4G1"

    ENCODEmetadata <- downloadENCODEknockdownMetadata(cellLine, gene)
    expect_is(ENCODEmetadata, "tbl_df")
    expect_identical(ENCODEmetadata$`Experiment target`, rep("EIF4G1", 2))

    counts <- loadENCODEgeneExpression(ENCODEmetadata)
    expect_identical(colnames(counts),
                     c("gene_id", "transcript_id(s)",
                       paste0("shRNA", 1:2), paste0("control", 1:2)))
})
