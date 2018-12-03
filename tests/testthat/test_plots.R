context("Plot results")

# data(compareKnockdown)
data("l1000perturbationsSmallMolecules")
perturbations <- l1000perturbationsSmallMolecules
data("diffExprStat")
cellLine <- "HEPG2"

compareKnockdown <- list()
compareKnockdown$pearson  <- compareAgainstL1000(diffExprStat, perturbations,
                                                 cellLine, method="pearson")
compareKnockdown$spearman <- compareAgainstL1000(diffExprStat, perturbations,
                                                 cellLine, method="spearman")
compareKnockdown$gsea     <- compareAgainstL1000(diffExprStat, perturbations,
                                                 cellLine, method="gsea")
condition <- compareKnockdown$gsea$compounds[[1]]

test_that("Plot Spearman correlation", {
    plot <- plotL1000comparison(compareKnockdown$spearman, condition)
    expect_is(plot, "NULL")
})

test_that("Plot Pearson correlation", {
    plot <- plotL1000comparison(compareKnockdown$pearson, condition)
    expect_is(plot, "NULL")
})

test_that("Plot GSEA", {
    plot <- plotL1000comparison(compareKnockdown$gsea, condition)
    expect_is(plot, "ggplot")
})
