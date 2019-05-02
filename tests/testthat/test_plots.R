context("Plot results")

# data(compareKnockdown)
data("cmapPerturbationsCompounds")
perturbations <- cmapPerturbationsCompounds
data("diffExprStat")
cellLine <- "HEPG2"

compareKnockdown <- list()
compareKnockdown$pearson  <- compareAgainstCMap(diffExprStat, perturbations,
                                                cellLine, method="pearson")
compareKnockdown$spearman <- compareAgainstCMap(diffExprStat, perturbations,
                                                cellLine, method="spearman")
compareKnockdown$gsea     <- compareAgainstCMap(diffExprStat, perturbations,
                                                cellLine, method="gsea")
condition <- compareKnockdown$gsea$compounds[[1]]

test_that("Plot Spearman correlation", {
    plot <- plotCMapComparison(compareKnockdown$spearman, condition)
    expect_is(plot, "NULL")
})

test_that("Plot Pearson correlation", {
    plot <- plotCMapComparison(compareKnockdown$pearson, condition)
    expect_is(plot, "NULL")
})

test_that("Plot GSEA", {
    plot <- plotCMapComparison(compareKnockdown$gsea, condition)
    expect_is(plot, "ggplot")
})
