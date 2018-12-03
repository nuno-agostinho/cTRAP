context("Working with L1000 data")

cellLine <- "HepG2"
data("l1000perturbationsSmallMolecules")
perturbations <- l1000perturbationsSmallMolecules
data("diffExprStat")

test_that("L1000 metadata and conditions are loaded", {
    data("l1000metadata")
    # metadata <- downloadL1000data("L1000metadata.txt", "metadata")
    # metadata <- filterL1000metadata(metadata, cellLine = cellLine,
    #                                 timepoint = "2 h")

    expect_is(l1000metadata, "data.table")
    expect_true(all(c("sig_id", "pert_id", "cell_id",
                      "pert_idose", "pert_itime") %in% colnames(l1000metadata)))

    conditions <- getL1000conditions(l1000metadata)
    expect_is(conditions, "list")
    expect_identical(names(conditions), c("Perturbation type", "Cell line",
                                          "Dosage", "Time points"))

    expect_true(!is.null(attr(l1000perturbationsSmallMolecules, "metadata")))
})

test_that("Perturbation types are retrievable", {
    pertTypes <- getL1000perturbationTypes()
    expect_is(pertTypes, "character")
    expect_true(all(c("trt_cp", "trt_sh.cgs") %in% pertTypes))
})

test_that("Compare using Spearman correlation", {
    data <- compareAgainstL1000(diffExprStat, perturbations, cellLine,
                                method="spearman")
    expect_is(data, "l1000comparison")
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
    data <- compareAgainstL1000(diffExprStat, perturbations, cellLine,
                                method="pearson")
    expect_is(data, "l1000comparison")
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
    data <- compareAgainstL1000(diffExprStat, perturbations, cellLine,
                                method="gsea")
    expect_is(data, "l1000comparison")
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
