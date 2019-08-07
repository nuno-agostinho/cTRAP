context("Working with drug sensitivity data")

data("diffExprStat")
gdsc <- loadExpressionDrugSensitivityAssociation("GDSC")

test_that("Correctly load expression and drug sensitivity association", {
    expect_is(gdsc, "matrix")
    attrs <- c("cellLines", "cellLineInfo", "compoundInfo", "source")
    expect_true(all(attrs %in% names(attributes(gdsc))))
    expect_equal(attr(gdsc, "source"), "GDSC release-7.0")
    expect_equal(length(attr(gdsc, "cellLines")), 983)
})

predicted <- predictTargetingDrugs(diffExprStat, gdsc)

test_that("Predict targeting drug", {
    expect_is(predicted, "referenceComparison")
    expect_is(predicted, "targetingDrugs")
    expect_equal(colnames(predicted)[1], "compound")
})

test_that("Plot predicted targeting drugs", {
    plot <- plot(predicted, method="spearman")
    expect_is(plot, "ggplot")
    expect_equal(plot$labels$y, "Spearman's correlation coefficient")

    plot <- plot(predicted, method="pearson")
    expect_equal(plot$labels$y, "Pearson's correlation coefficient")

    plot <- plot(predicted, method="gsea")
    expect_equal(plot$labels$y, "Weighted Connectivity Score (WCTS)")
})
