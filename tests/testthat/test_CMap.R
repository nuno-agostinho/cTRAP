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

    expect_true(!is.null(attr(perturbations, "metadata")))
})

test_that("Perturbation types are retrievable", {
    pertTypes <- getCMapPerturbationTypes()
    expect_is(pertTypes, "character")
    expect_true(all(c("trt_cp", "trt_sh.cgs") %in% pertTypes))
})

test_that("Subsetting a perturbationChanges object", {
    ncol <- 22
    nrow <- 12328
    expect_equal(dim(perturbations), c(nrow, ncol))
    expect_equal(nrow(perturbations), nrow)
    expect_equal(ncol(perturbations), ncol)

    expect_identical(head(rownames(perturbations)),
                     c("PSME1", "ATF1", "RHEB", "FOXO3", "RHOA", "IL1B"))
    expect_identical(head(colnames(perturbations)),
                     c("CVD001_HEPG2_24H:BRD-A14014306-001-01-1:4.1",
                       "CVD001_HEPG2_24H:BRD-A65142661-034-01-8:5.35",
                       "CVD001_HEPG2_24H:BRD-K31030218-001-01-1:4.25",
                       "CVD001_HEPG2_24H:BRD-K41172353-001-01-4:4.7",
                       "CVD001_HEPG2_24H:BRD-K60476892-001-02-1:4.1072",
                       "CVD001_HEPG2_24H:BRD-K62810658-001-05-6:4.6768"))
    expect_identical(dimnames(perturbations),
                     list(rownames(perturbations), colnames(perturbations)))

    # Missing i and implicitly missing j
    perts <- perturbations[]
    expect_identical(rownames(perts), rownames(perturbations))
    expect_identical(colnames(perts), colnames(perturbations))
    expect_identical(perts, perturbations)

    # Missing i and explicitly missing j
    perts <- perturbations[ , ]
    expect_identical(rownames(perts), rownames(perturbations))
    expect_identical(colnames(perts), colnames(perturbations))
    expect_identical(perts, perturbations)

    # # Given i and implicitly missing j
    # perts <- perturbations[5:8]
    # expect_identical(rownames(perts), rownames(perturbations))
    # expect_identical(colnames(perts), colnames(perturbations)[5:8])

    # Given i and explicitly missing j
    perts <- perturbations[13:19, ]
    expect_identical(rownames(perts), rownames(perturbations)[13:19])
    expect_identical(colnames(perts), colnames(perturbations))

    # Explicitly missing i and given j
    perts <- perturbations[ , 4:7]
    expect_identical(rownames(perts), rownames(perturbations))
    expect_identical(colnames(perts), colnames(perturbations)[4:7])

    # Given i and j
    perts <- perturbations[21:22, 11:22]
    expect_identical(rownames(perts), rownames(perturbations)[21:22])
    expect_identical(colnames(perts), colnames(perturbations)[11:22])

    # Given i and j
    perts <- perturbations[6:14, 3:9]
    expect_identical(rownames(perts), rownames(perturbations)[6:14])
    expect_identical(colnames(perts), colnames(perturbations)[3:9])

    # Subscript out of bounds - giving both excessive i and j
    expect_error(perturbations[seq(nrow(perturbations) + 10),
                               seq(ncol(perturbations) + 10)])
    # Subscript out of bounds - only giving excessive i
    expect_error(perturbations[seq(nrow(perturbations) + 10),
                               seq(ncol(perturbations))])
    # Subscript out of bounds - only giving excessive j
    expect_error(perturbations[seq(nrow(perturbations)),
                               seq(ncol(perturbations) + 10)])
    # Subscript out of bounds - only giving excessive i and no j
    expect_error(perturbations[seq(nrow(perturbations) + 10)])

    # Extract based on character strings
    rows <- c("ABHD4", "SFN", "SRPRB", "DERA")
    cols <- c("CVD001_HEPG2_24H:BRD-K96188950-001-04-5:4.3967",
              "CVD001_HUH7_24H:BRD-A14014306-001-01-1:4.1",
              "CVD001_HUH7_24H:BRD-A65142661-034-01-8:5.35")
    perts <- perturbations[rows, cols]
    expect_identical(rownames(perts), rows)
    expect_identical(colnames(perts), cols)

    # Subscript out of bounds - no i found based on character
    expect_error(perturbations["ddffd"])
    # Subscript out of bounds - no j found based on character
    expect_error(perturbations[ , c(cols, "test")])
})

checkPerturbationIdentifierInformationInMetadata <- function(data) {
    # All information for perturbation identifiers are available in metadata
    metadata <- attr(data, "metadata")
    expect_true(!is.null(metadata))
    expect_true(all(data$compound_perturbation %in% metadata$sig_id))
    expect_equal(length(data$compound_perturbation),
                 length(metadata$sig_id))
}

testSimilarPerturbationsResult <- function(data, col, cols) {
    expect_is(data, "referenceComparison")
    expect_is(data, "similarPerturbations")
    expect_identical(colnames(data), cols)
    checkPerturbationIdentifierInformationInMetadata(data)

    # Data is ordered based on Spearman's rank
    expect_equal(data[order(data[[col]])]$compound_perturbation,
                 data$compound_perturbation)

    # Expect average perturbations across cell lines
    expect_true(all(parseCMapID(data$compound_perturbation) %in%
                        data$compound_perturbation))
}

test_that("Compare using Spearman correlation", {
    data <- rankSimilarPerturbations(diffExprStat, perturbations,
                                     method="spearman")
    cols <- c("compound_perturbation", "spearman_coef", "spearman_pvalue",
              "spearman_qvalue", "spearman_rank")
    testSimilarPerturbationsResult(data, col="spearman_rank", cols=cols)
})

test_that("Compare using Pearson correlation", {
    data <- rankSimilarPerturbations(diffExprStat, perturbations,
                                     method="pearson")
    cols <- c("compound_perturbation", "pearson_coef", "pearson_pvalue",
              "pearson_qvalue", "pearson_rank")
    testSimilarPerturbationsResult(data, col="pearson_rank", cols=cols)
})

test_that("Compare using GSEA", {
    data <- rankSimilarPerturbations(diffExprStat, perturbations, method="gsea")
    cols <- c("compound_perturbation", "GSEA", "GSEA_rank")
    testSimilarPerturbationsResult(data, col="GSEA_rank", cols=cols)
})

test_that("Compare against CMap using multiple methods simultaneously", {
    method <- c("pearson", "spearman", "gsea")
    data <- rankSimilarPerturbations(diffExprStat, perturbations, method=method)
    checkPerturbationIdentifierInformationInMetadata(data)

    checkIfAnyColsHaveMethodsName <- function(col, data) {
        any(startsWith(tolower(colnames(data)), col))
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
    data <- rankSimilarPerturbations(diffExprStat, perturbations,
                                     rankPerCellLine=FALSE, cellLineMean=TRUE)
    checkPerturbationIdentifierInformationInMetadata(data)
    expect_true(areNAsIncludedInRanks(data))

    # Expect NO missing values in rankings when individual cell lines are
    # ranked as well
    data <- rankSimilarPerturbations(diffExprStat, perturbations,
                                     rankPerCellLine=TRUE)
    checkPerturbationIdentifierInformationInMetadata(data)
    expect_false(areNAsIncludedInRanks(data))

    # Expect NO missing values in rankings if cell line means are NOT calculated
    data <- rankSimilarPerturbations(diffExprStat, perturbations,
                                     rankPerCellLine=FALSE, cellLineMean=FALSE)
    checkPerturbationIdentifierInformationInMetadata(data)
    expect_false(areNAsIncludedInRanks(data))

    data <- rankSimilarPerturbations(diffExprStat, perturbations,
                                     rankPerCellLine=TRUE,
                                     cellLineMean=FALSE)
    checkPerturbationIdentifierInformationInMetadata(data)
    expect_false(areNAsIncludedInRanks(data))
})

test_that("Compare against CMap based on a gene set", {
    geneset <- sample(names(diffExprStat), 200)
    rank <- rankSimilarPerturbations(geneset, perturbations)
    checkPerturbationIdentifierInformationInMetadata(rank)

    cols <- c("compound_perturbation", "GSEA", "GSEA_rank")
    testSimilarPerturbationsResult(rank, col="GSEA_rank", cols=cols)

    input <- attr(rank, "input")
    expect_equivalent(input, geneset)
    expect_true(attr(input, "isGeneset"))

    # Automatically perform only GSEA based on a gene set
    expect_warning(
        rankSimilarPerturbations(geneset, perturbations, method="spearman"),
        "gsea")
})
