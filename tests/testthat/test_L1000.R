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
        "genes", "HepG2_t_spearman_coef", "HepG2_t_spearman_pvalue",
        "HepG2_t_spearman_qvalue", "Average_t_spearman_coef"))
})

test_that("Compare using Pearson correlation", {
    data <- compareAgainstL1000(diffExprStat, perturbations, cellLine,
                                method="pearson")
    expect_is(data, "l1000comparison")
    expect_identical(colnames(data), c(
        "genes", "HepG2_t_pearson_coef", "HepG2_t_pearson_pvalue",
        "HepG2_t_pearson_qvalue", "Average_t_pearson_coef"))
})

test_that("Compare using GSEA", {
    data <- compareAgainstL1000(diffExprStat, perturbations, cellLine,
                                method="gsea")
    expect_is(data, "l1000comparison")
    expect_identical(colnames(data), c("genes", "HepG2_WTCS", "Average_WTCS"))
})
