context("Working with drug sensitivity data")

data("diffExprStat")
gdsc <- loadExpressionDrugSensitivityAssociation("GDSC")

test_that("Correctly load expression and drug sensitivity association", {
    expect_is(gdsc, "matrix")
    attrs <- c("cellLines", "cellLineInfo", "compoundInfo", "source")
    expect_true(all(attrs %in% names(attributes(gdsc))))
    expect_equal(attr(gdsc, "source"), "GDSC")
    expect_equal(length(attr(gdsc, "cellLines")), 983)
})

test_that("Predict targeting drug", {
    predicted <- predictTargetingDrugs(diffExprStat, gdsc)
    expect_is(predicted, "dataComparison")
    expect_is(predicted, "targetingDrugs")
    expect_equal(colnames(predicted)[1], "compound")

    plot <- plot(predicted, method="spearman")
    expect_is(plot, "ggplot")
    expect_equal(plot$labels$y, "Spearman's correlation coefficient")

    plot <- plot(predicted, method="pearson")
    expect_equal(plot$labels$y, "Pearson's correlation coefficient")

    plot <- plot(predicted, method="gsea")
    expect_equal(plot$labels$y, "Weighted Connectivity Score (WCTS)")
})
