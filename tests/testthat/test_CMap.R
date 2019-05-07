context("Working with CMap data")

data("cmapPerturbationsCompounds")
perturbations <- cmapPerturbationsCompounds
data("diffExprStat")

test_that("CMap metadata and conditions are loaded", {
    data("cmapMetadata")
    # metadata <- loadCMapData("cmapMetadata.txt", "metadata")
    # metadata <- filterCMapMetadata(metadata, cellLine="HepG2",
    #                                timepoint="2 h")

    expect_is(cmapMetadata, "data.table")
    expect_true(all(c("sig_id", "pert_id", "cell_id",
                      "pert_idose", "pert_itime") %in% colnames(cmapMetadata)))

    conditions <- getCMapConditions(cmapMetadata)
    expect_is(conditions, "list")
    expect_identical(names(conditions), c("perturbationType", "cellLine",
                                          "dosage", "timepoint"))

    expect_true(!is.null(attr(cmapPerturbationsCompounds, "metadata")))
})

test_that("Perturbation types are retrievable", {
    pertTypes <- getCMapPerturbationTypes()
    expect_is(pertTypes, "character")
    expect_true(all(c("trt_cp", "trt_sh.cgs") %in% pertTypes))
})

test_that("Compare using Spearman correlation", {
    data <- compareAgainstCMap(diffExprStat, perturbations, method="spearman")
    expect_is(data, "cmapComparison")
    expect_identical(colnames(data), c(
        "compound_perturbation", "spearman_coef", "spearman_pvalue",
        "spearman_qvalue", "spearman_rank"))

    expect_equal(head(data$compound_perturbation),
                 c("CVD001_24H:BRD-A14014306-001-01-1:4.1",
                   "CVD001_24H:BRD-A65142661-034-01-8:5.35",
                   "CVD001_24H:BRD-K31030218-001-01-1:4.25",
                   "CVD001_24H:BRD-K41172353-001-01-4:4.7",
                   "CVD001_24H:BRD-K60476892-001-02-1:4.1072",
                   "CVD001_24H:BRD-K62810658-001-05-6:4.6768"))

    expect_equal(head(data[order(data$spearman_rank)]$compound_perturbation),
                 c("CVD001_24H:BRD-A14014306-001-01-1:4.1",
                   "CVD001_24H:BRD-K84595254-001-03-0:4.9444",
                   "CVD001_24H:BRD-K41172353-001-01-4:4.7",
                   "CVD001_24H:BRD-K84389640-001-01-5:4.225",
                   "CVD001_24H:BRD-K96188950-001-04-5:4.3967",
                   "CVD001_24H:BRD-K60476892-001-02-1:4.1072"))
})

test_that("Compare using Pearson correlation", {
    data <- compareAgainstCMap(diffExprStat, perturbations, method="pearson")
    expect_is(data, "cmapComparison")
    expect_identical(colnames(data), c(
        "compound_perturbation", "pearson_coef", "pearson_pvalue",
        "pearson_qvalue", "pearson_rank"))

    expect_equal(head(data$compound_perturbation),
                 c("CVD001_24H:BRD-A14014306-001-01-1:4.1",
                   "CVD001_24H:BRD-A65142661-034-01-8:5.35",
                   "CVD001_24H:BRD-K31030218-001-01-1:4.25",
                   "CVD001_24H:BRD-K41172353-001-01-4:4.7",
                   "CVD001_24H:BRD-K60476892-001-02-1:4.1072",
                   "CVD001_24H:BRD-K62810658-001-05-6:4.6768"))

    expect_equal(head(data[order(data$pearson_rank), ]$compound_perturbation),
                 c("CVD001_24H:BRD-A14014306-001-01-1:4.1",
                   "CVD001_24H:BRD-K84595254-001-03-0:4.9444",
                   "CVD001_24H:BRD-K41172353-001-01-4:4.7",
                   "CVD001_24H:BRD-K84389640-001-01-5:4.225",
                   "CVD001_24H:BRD-K96188950-001-04-5:4.3967",
                   "CVD001_24H:BRD-K77508012-001-01-9:6.025"))
})

test_that("Compare using GSEA", {
    data <- compareAgainstCMap(diffExprStat, perturbations, method="gsea")
    expect_is(data, "cmapComparison")
    expect_identical(colnames(data),
                     c("compound_perturbation", "WTCS", "gsea_rank"))
    expect_equal(head(data$compound_perturbation),
                 c("CVD001_24H:BRD-A14014306-001-01-1:4.1",
                   "CVD001_24H:BRD-A65142661-034-01-8:5.35",
                   "CVD001_24H:BRD-K31030218-001-01-1:4.25",
                   "CVD001_24H:BRD-K41172353-001-01-4:4.7",
                   "CVD001_24H:BRD-K60476892-001-02-1:4.1072",
                   "CVD001_24H:BRD-K62810658-001-05-6:4.6768"))
    expect_equal(head(data[order(data$gsea_rank), ]$compound_perturbation),
                 c("CVD001_24H:BRD-A14014306-001-01-1:4.1",
                   "CVD001_24H:BRD-K84595254-001-03-0:4.9444",
                   "CVD001_24H:BRD-K31030218-001-01-1:4.25",
                   "CVD001_24H:BRD-K41172353-001-01-4:4.7",
                   "CVD001_24H:BRD-K84389640-001-01-5:4.225",
                   "CVD001_24H:BRD-K62810658-001-05-6:4.6768"))
})

test_that("Compare against CMap using multiple methods simultaneously", {
    method <- c("pearson", "spearman", "gsea")
    data <- compareAgainstCMap(diffExprStat, perturbations, method=method)

    checkIfAnyColsHaveMethodsName <- function(col, data) {
        any(startsWith(colnames(data), col))
    }

    expect_true(all(vapply(
        method, checkIfAnyColsHaveMethodsName, logical(1), data)))
})



test_that("Compare against CMap by also ranking individual cell lines", {
    areNAsIncludedInRanks <- function(data) {
        anyNAs  <- function(i) any(is.na(data[[i]]))
        rankNAs <- all(vapply(colnames(data)[endsWith(colnames(data), "_rank")],
                              anyNAs, logical(1)))
        return(rankNAs)
    }

    # Expect missing values in rankings when individual cell lines are not
    # ranked and if means across cell lines are calculated
    data <- compareAgainstCMap(diffExprStat, perturbations,
                               rankIndividualCellLines=FALSE,
                               cellLineMean=TRUE)
    expect_true(areNAsIncludedInRanks(data))

    # Expect NO missing values in rankings when individual cell lines are
    # ranked as well
    data <- compareAgainstCMap(diffExprStat, perturbations,
                               rankIndividualCellLines=TRUE)
    expect_false(areNAsIncludedInRanks(data))

    # Expect NO missing values in rankings if cell line means are NOT calculated
    data <- compareAgainstCMap(diffExprStat, perturbations,
                               rankIndividualCellLines=FALSE,
                               cellLineMean=FALSE)
    expect_false(areNAsIncludedInRanks(data))

    data <- compareAgainstCMap(diffExprStat, perturbations,
                               rankIndividualCellLines=TRUE,
                               cellLineMean=FALSE)
    expect_false(areNAsIncludedInRanks(data))
})
