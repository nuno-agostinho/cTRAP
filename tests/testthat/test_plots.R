context("Plot results")

# data(compareKnockdown)
data("cmapPerturbationsCompounds")
perturbations <- cmapPerturbationsCompounds
data("diffExprStat")

compareKnockdown <- compareAgainstCMap(diffExprStat, perturbations)
condition <- compareKnockdown[1, 1]

test_that("Plot Spearman's correlation", {
    plot <- plotCMapComparison(compareKnockdown, condition, "spearman")
    expect_is(plot, "NULL")
})

test_that("Plot Pearson's correlation", {
    plot <- plotCMapComparison(compareKnockdow, condition, "pearson")
    expect_is(plot, "NULL")
})

test_that("Plot GSEA", {
    plot <- plotCMapComparison(compareKnockdown, condition, "gsea")
    expect_is(plot, "ggplot")
})

