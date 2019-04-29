context("Working with CMap data")

cellLine <- "HepG2"
data("cmapPerturbationsSmallMolecules")
perturbations <- cmapPerturbationsSmallMolecules
data("diffExprStat")

test_that("CMap metadata and conditions are loaded", {
    data("cmapMetadata")
    # metadata <- loadCMapData("cmapMetadata.txt", "metadata")
    # metadata <- filterCMapMetadata(metadata, cellLine = cellLine,
    #                                timepoint = "2 h")

    expect_is(cmapMetadata, "data.table")
    expect_true(all(c("sig_id", "pert_id", "cell_id",
                      "pert_idose", "pert_itime") %in% colnames(cmapMetadata)))

    conditions <- getCMapConditions(cmapMetadata)
    expect_is(conditions, "list")
    expect_identical(names(conditions), c("Perturbation type", "Cell line",
                                          "Dosage", "Time points"))

    expect_true(!is.null(attr(cmapPerturbationsSmallMolecules, "metadata")))
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
        "compounds", "HepG2_spearman_coef", "HepG2_spearman_pvalue",
        "HepG2_spearman_qvalue", "Average_spearman_coef"))

    expect_equal(head(data$compounds),
                 c("CVD001_HEPG2_24H:BRD-A14014306-001-01-1:4.1",
                   "CVD001_HEPG2_24H:BRD-A65142661-034-01-8:5.35",
                   "CVD001_HEPG2_24H:BRD-K31030218-001-01-1:4.25",
                   "CVD001_HEPG2_24H:BRD-K41172353-001-01-4:4.7",
                   "CVD001_HEPG2_24H:BRD-K60476892-001-02-1:4.1072",
                   "CVD001_HEPG2_24H:BRD-K62810658-001-05-6:4.6768"))

    expect_equal(head(data[order(data$HepG2_spearman_coef), ]$compounds),
                 c("CVD001_HEPG2_24H:BRD-A65142661-034-01-8:5.35",
                   "CVD001_HEPG2_24H:BRD-K31030218-001-01-1:4.25",
                   "CVD001_HEPG2_24H:BRD-K94818765-001-01-0:4.8",
                   "CVD001_HEPG2_24H:BRD-K77508012-001-01-9:6.025",
                   "CVD001_HEPG2_24H:BRD-K60476892-001-02-1:4.1072",
                   "CVD001_HEPG2_24H:BRD-K84389640-001-01-5:4.225"))
})

test_that("Compare using Pearson correlation", {
    data <- compareAgainstCMap(diffExprStat, perturbations, method="pearson")
    expect_is(data, "cmapComparison")
    expect_identical(colnames(data), c(
        "compounds", "HepG2_pearson_coef", "HepG2_pearson_pvalue",
        "HepG2_pearson_qvalue", "Average_pearson_coef"))

    expect_equal(head(data$compounds),
                 c("CVD001_HEPG2_24H:BRD-A14014306-001-01-1:4.1",
                   "CVD001_HEPG2_24H:BRD-A65142661-034-01-8:5.35",
                   "CVD001_HEPG2_24H:BRD-K31030218-001-01-1:4.25",
                   "CVD001_HEPG2_24H:BRD-K41172353-001-01-4:4.7",
                   "CVD001_HEPG2_24H:BRD-K60476892-001-02-1:4.1072",
                   "CVD001_HEPG2_24H:BRD-K62810658-001-05-6:4.6768"))

    expect_equal(head(data[order(data$HepG2_pearson_coef), ]$compounds),
                 c("CVD001_HEPG2_24H:BRD-A65142661-034-01-8:5.35",
                   "CVD001_HEPG2_24H:BRD-K31030218-001-01-1:4.25",
                   "CVD001_HEPG2_24H:BRD-K94818765-001-01-0:4.8",
                   "CVD001_HEPG2_24H:BRD-K77508012-001-01-9:6.025",
                   "CVD001_HEPG2_24H:BRD-K84389640-001-01-5:4.225",
                   "CVD001_HEPG2_24H:BRD-K60476892-001-02-1:4.1072"))
})

test_that("Compare using GSEA", {
    data <- compareAgainstCMap(diffExprStat, perturbations, method="gsea")
    expect_is(data, "cmapComparison")
    expect_identical(colnames(data),
                     c("compounds", "HepG2_WTCS", "Average_WTCS"))

    expect_equal(head(data$compounds),
                 c("CVD001_HEPG2_24H:BRD-A14014306-001-01-1:4.1",
                   "CVD001_HEPG2_24H:BRD-A65142661-034-01-8:5.35",
                   "CVD001_HEPG2_24H:BRD-K31030218-001-01-1:4.25",
                   "CVD001_HEPG2_24H:BRD-K41172353-001-01-4:4.7",
                   "CVD001_HEPG2_24H:BRD-K60476892-001-02-1:4.1072",
                   "CVD001_HEPG2_24H:BRD-K62810658-001-05-6:4.6768"))

    expect_equal(head(data[order(data$HepG2_WTCS), ]$compounds),
                 c("CVD001_HEPG2_24H:BRD-A65142661-034-01-8:5.35",
                   "CVD001_HEPG2_24H:BRD-K96188950-001-04-5:4.3967",
                   "CVD001_HEPG2_24H:BRD-K94818765-001-01-0:4.8",
                   "CVD001_HEPG2_24H:BRD-K77508012-001-01-9:6.025",
                   "CVD001_HEPG2_24H:BRD-K31030218-001-01-1:4.25",
                   "CVD001_HEPG2_24H:BRD-K41172353-001-01-4:4.7"))
})
