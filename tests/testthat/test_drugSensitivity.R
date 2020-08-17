context("Working with drug sensitivity data")

data("diffExprStat")
gdsc <- loadExpressionDrugSensitivityAssociation("GDSC 7")

test_that("Correctly load expression and drug sensitivity association", {
    expect_is(gdsc, "matrix")
    attrs <- c("cellLines", "cellLineInfo", "compoundInfo", "source")
    expect_true(all(attrs %in% names(attributes(gdsc))))
    expect_equal(attr(gdsc, "source"), "GDSC 7")
    expect_equal(attr(gdsc, "drugActivityMetric"), "-log(IC50)")
    expect_equal(attr(gdsc, "method"), "spearman")
    expect_equal(attr(gdsc, "type"), "compounds")
    expect_true(attr(gdsc, "isDrugActivityDirectlyProportionalToSensitivity"))
    expect_equal(length(attr(gdsc, "cellLines")), 983)
    expect_true(!is.null(attr(gdsc, "cellLineInfo")))
})

test_that("Predict targeting drug", {
    predicted <- predictTargetingDrugs(diffExprStat, gdsc)
    expect_is(predicted, "referenceComparison")
    expect_is(predicted, "targetingDrugs")
    expect_equal(colnames(predicted)[1], "compound")
})

test_that("Predict targeting drug based on gene set", {
    geneset <- sample(names(diffExprStat), 200)
    predicted <- predictTargetingDrugs(geneset, gdsc)

    expect_is(predicted, "referenceComparison")
    expect_is(predicted, "targetingDrugs")
    expect_equal(colnames(predicted)[1], "compound")

    input <- attr(predicted, "input")
    expect_equivalent(input, geneset)
    expect_true(attr(input, "isGeneset"))

    # Automatically perform only GSEA based on a gene set
    expect_warning(predictTargetingDrugs(geneset, gdsc, method="spearman"),
                   "gsea")
})

test_that("Plot predicted targeting drugs", {
    predicted <- predictTargetingDrugs(diffExprStat, gdsc)
    plot <- plot(predicted, method="spearman")
    expect_is(plot, "ggplot")
    expect_equal(plot$labels$y, "Spearman's correlation coefficient")

    plot <- plot(predicted, method="pearson")
    expect_equal(plot$labels$y, "Pearson's correlation coefficient")

    plot <- plot(predicted, method="gsea")
    expect_equal(plot$labels$y, "Weighted Connectivity Score (WTCS)")
})
