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

test_that("Download ENCODE knockdown metadata", {
    skip_on_bioc()
    metadata <- downloadENCODEknockdownMetadata("HepG2", "EIF4G1")
    expect_equal(nrow(metadata), 2)
    expect_equal(unique(metadata$`File format`), "tsv")
    expect_equal(unique(metadata$`Output type`), "gene quantifications")
    expect_equal(unique(metadata$Assay), "shRNA RNA-seq")
    expect_equal(unique(metadata$`Biosample term name`),  "HepG2")
    expect_equal(unique(metadata$`Biosample genetic modifications categories`),
                 "interference")
    expect_equal(unique(metadata$`Biosample genetic modifications methods`),
                 "RNAi")
    expect_equal(unique(metadata$`Experiment target`), "EIF4G1")
    expect_equal(unique(metadata$`Library made from`), "polyadenylated mRNA")
    expect_equal(unique(metadata$Project), "ENCODE")
    expect_equal(unique(metadata$`Biological replicate(s)`), c("1", "2"))
    expect_equal(unique(metadata$`File assembly`), "hg19")
    expect_equal(unique(metadata$`File Status`), "released")
})

# test_that("Perform differential gene expression", {
#     data("ENCODEsamples")
#     counts <- prepareENCODEgeneExpression(ENCODEsamples)
#
#     diffExpr <- performDifferentialExpression(counts)
#     expect_is(diffExpr, "data.frame")
#     expect_identical(colnames(diffExpr), c("Gene_symbol", "logFC", "AveExpr",
#                                            "t", "P.Value", "adj.P.Val", "B"))
# })
