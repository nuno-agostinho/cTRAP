context("Working with L1000 data")

cellLine <- "HepG2"
perturbations <- l1000:::l1000perturbationsSmallMolecules
diffExprGenes <- l1000:::diffExprStat

test_that("L1000 metadata and conditions are loaded", {
    metadata <- loadL1000metadata("L1000metadata.txt")
    expect_is(metadata, "data.table")
    expect_true(c("sig_id", "pert_id", "cell_id",
                  "pert_idose", "pert_itime") %in% colnames(metadata))

    conditions <- getL1000Conditions("L1000metadata.txt")
    expect_is(conditions, "list")
    expect_is(names(conditions), c("Perturbation type", "Cell line", "Dosage",
                                   "Time points"))
})

test_that("Perturbation types are retrievable", {
    pertTypes <- getL1000PerturbationTypes()
    expect_is(pertTypes, "character")
    expect_true(all(c("trt_cp", "trt_sh.cgs") %in% pertTypes))
})

test_that("Compare using Spearman correlation", {
    data <- compareAgainstL1000(diffExprGenes, perturbations, cellLine,
                                method="spearman")
    expect_is(data, "L1000comparison")
    expect_is(colnames(data), c(
        "genes", "HepG2_t_spearman_coef", "HepG2_t_spearman_pvalue",
        "HepG2_t_spearman_qvalue", "Average_t_spearman_coef"))
})

test_that("Compare using Pearson correlation", {
    data <- compareAgainstL1000(diffExprGenes, perturbations, cellLine,
                                method="pearson")
    expect_is(data, "L1000comparison")
    expect_is(colnames(data), c(
        "genes", "HepG2_t_pearson_coef", "HepG2_t_pearson_pvalue",
        "HepG2_t_pearson_qvalue", "Average_t_pearson_coef"))
})

test_that("Compare using GSEA", {
    data <- compareAgainstL1000(diffExprGenes, perturbations, cellLine,
                                method="gsea")
    expect_is(data, "L1000comparison")
    expect_is(colnames(data), c("genes", "HepG2_WTCS", "Average_WTCS"))
})
