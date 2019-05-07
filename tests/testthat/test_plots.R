context("Plot results")

# data(compareKnockdown)
data("cmapPerturbationsCompounds")
perturbations <- cmapPerturbationsCompounds
data("diffExprStat")

compareKnockdown <- compareAgainstCMap(diffExprStat, perturbations)
condition <- compareKnockdown[1, 1]

test_that("Plot similar perturbations according to selected method", {
    plot <- plot(compareKnockdown, method="pearson")
    expect_s3_class(plot, "ggplot")
    expect_match(plot$labels$y, "Pearson")

    plot <- plot(compareKnockdown, method="spearman")
    expect_s3_class(plot, "ggplot")
    expect_match(plot$labels$y, "Spearman")

    plot <- plot(compareKnockdown, method="gsea")
    expect_s3_class(plot, "ggplot")
    expect_match(plot$labels$y, "WCTS")

    # Avoid plotting based on methods not available in the cmapComparison object
    compareKnockdownPearson <- compareAgainstCMap(diffExprStat, perturbations,
                                                  method="pearson")
    expect_error(plot(compareKnockdownPearson, method="gsea"))
})

test_that("Plot custom top and bottom similar perturbations", {
    # Top 7 most similar perturbations
    plotTopBottomSimilarPerts <- function(comp, top, bottom) {
        plot <- plot(comp, n=c(top, bottom))
        expect_s3_class(plot, "ggplot")
        expect_equal(length(unique(plot$data$label)) - 1, sum(top, bottom))
    }
    plotTopBottomSimilarPerts(compareKnockdown, 7, 0)

    # Bottom 4 least similar perturbations
    plotTopBottomSimilarPerts(compareKnockdown, 0, 4)

    # Top 5 most and bottom 2 least similar perturbations
    plotTopBottomSimilarPerts(compareKnockdown, 5, 2)
})

test_that("Plot similar perturbations with/without metadata", {
    plot1 <- plot(compareKnockdown, showMetadata=TRUE)
    plot2 <- plot(compareKnockdown, showMetadata=FALSE)

    expect_s3_class(plot1, "ggplot")
    expect_s3_class(plot2, "ggplot")

    metadataLabel <- unique(plot1$data$label)
    originalLabel <- unique(plot2$data$label)
    expect_false(all(metadataLabel != originalLabel))
})

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
