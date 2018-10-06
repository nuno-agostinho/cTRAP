context("Plot results")

compareKnockdown <- loadInternalData("compareKnockdown")
EIF4G1knockdown <- grep("EIF4G1", compareKnockdown$gsea_ordered$genes,
                        value=TRUE)

test_that("Plot Spearman correlation", {
    plot <- plotL1000comparison(compareKnockdown$spearman, EIF4G1knockdown)
    expect_is(plot, "NULL")
})

test_that("Plot Pearson correlation", {
    plot <- plotL1000comparison(compareKnockdown$pearson, EIF4G1knockdown)
    expect_is(plot, "NULL")
})

test_that("Plot GSEA", {
    plot <- plotL1000comparison(compareKnockdown$gsea, EIF4G1knockdown,
                                topGenes=150)
    expect_is(plot, "ggplot")

    plot <- plotL1000comparison(compareKnockdown$gsea, EIF4G1knockdown,
                                topGenes=20)
    expect_is(plot, "ggplot")
})

