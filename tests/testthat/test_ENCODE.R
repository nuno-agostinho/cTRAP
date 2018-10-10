context("Loading ENCODE data")

data("ENCODEmetadata")
data("counts")

test_that("Load ENCODE gene expression for EIF4G1", {
    expect_is(ENCODEmetadata, "tbl_df")
    expect_identical(ENCODEmetadata$`Experiment target`, rep("EIF4G1", 2))

    expect_identical(colnames(counts),
                     c("gene_id", "transcript_id(s)",
                       paste0("shRNA", 1:2), paste0("control", 1:2)))
})

test_that("Perform differential gene expression", {
    diffExpr <- performDifferentialExpression(counts)
    expect_is(diffExpr, "data.frame")
    expect_identical(colnames(diffExpr), c("Gene_symbol", "logFC", "AveExpr",
                                           "t", "P.Value", "adj.P.Val", "B"))
})
