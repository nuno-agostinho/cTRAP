context("Plot results")

data("cmapPerturbationsCompounds")
perturbations <- cmapPerturbationsCompounds
data("diffExprStat")

cmp <- compareAgainstCMap(diffExprStat, perturbations)
condition <- cmp[1, 1]

# Plot cmapComparison object ---------------------------------------------------

test_that("Plot similar perturbations according to selected method", {
    plot <- plot(cmp, method="pearson")
    expect_s3_class(plot, "ggplot")
    expect_match(plot$labels$y, "Pearson")

    plot <- plot(cmp, method="spearman")
    expect_s3_class(plot, "ggplot")
    expect_match(plot$labels$y, "Spearman")

    plot <- plot(cmp, method="gsea")
    expect_s3_class(plot, "ggplot")
    expect_match(plot$labels$y, "WCTS")

    # Avoid plotting based on methods not available in the cmapComparison object
    cmpPearson <- compareAgainstCMap(diffExprStat, perturbations,
                                     method="pearson")
    expect_error(plot(cmpPearson, method="gsea"))
})

test_that("Plot custom top and bottom similar perturbations", {
    # Top 7 most similar perturbations
    plotTopBottomSimilarPerts <- function(comp, top, bottom) {
        plot <- plot(comp, n=c(top, bottom))
        expect_s3_class(plot, "ggplot")
        expect_equal(length(unique(plot$data$label)) - 1, sum(top, bottom))
    }
    plotTopBottomSimilarPerts(cmp, 7, 0)

    # Bottom 4 least similar perturbations
    plotTopBottomSimilarPerts(cmp, 0, 4)

    # Top 5 most and bottom 2 least similar perturbations
    plotTopBottomSimilarPerts(cmp, 5, 2)
})

test_that("Plot similar perturbations with/without metadata", {
    plot1 <- plot(cmp, showMetadata=TRUE)
    plot2 <- plot(cmp, showMetadata=FALSE)

    expect_s3_class(plot1, "ggplot")
    expect_s3_class(plot2, "ggplot")

    metadataLabel <- unique(plot1$data$label)
    originalLabel <- unique(plot2$data$label)
    expect_false(all(metadataLabel != originalLabel))
})

test_that("Plot only ranked perturbations", {
    plot <- plot(cmp, "spearman", plotNonRankedPerturbations=FALSE)
    expect_s3_class(plot, "ggplot")
    expect_false(plot$guides$colour)
})

test_that("Plot non-ranked perturbations", {
    plot <- plot(cmp, "spearman", plotNonRankedPerturbations=TRUE)
    expect_s3_class(plot, "ggplot")
    expect_identical(plot$labels$colour, "ranked")
    expect_null(plot$guides$colour)
})

# Plot cmapPerturbations object ------------------------------------------------

test_that("Plot Spearman's correlation for a single perturbation", {
    plot <- plot(perturbations, colnames(perturbations)[[1]], diffExprStat,
                 "spearman")
    expect_s3_class(plot, "ggplot")
    expect_equal(plot$labels$x, "Differentially expressed genes (input)")
})

test_that("Plot Pearson's correlation for a single perturbation", {
    plot <- plot(perturbations, colnames(perturbations)[[1]], diffExprStat,
                 "pearson")
    expect_s3_class(plot, "ggplot")
    expect_equal(plot$labels$x, "Differentially expressed genes (input)")
})

test_that("Plot GSEA for a single perturbation", {
    plot <- plot(perturbations, colnames(perturbations)[[1]], diffExprStat,
                 "gsea")
    expect_s3_class(plot, "ggplot")
})
