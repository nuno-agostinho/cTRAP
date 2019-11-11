context("Drug set enrichment")

library(data.table)

# Load and prepare drug descriptors
descriptors <- loadDrugDescriptors("NCI60", "2D")
drugSets    <- prepareDrugSets(descriptors)

test_that("Sets are prepared from a table of drug descriptors", {
    expect_is(drugSets, "list")
    expect_equal(attr(drugSets, "source"), "NCI60")
    expect_equal(attr(drugSets, "type"),   "2D")

    compoundInfo <- attr(drugSets, "compoundInfo")
    expect_true(!is.null(compoundInfo))
    expect_is(compoundInfo, "tbl")
})

# Match stats with drug sets identifier ----------------------------------------

selectDrugSet              <- drugSets$`FDA_status: Clinical trial`
drugSets_compoundInfo      <- attr(drugSets, "compoundInfo")
selectDrugSet_compoundInfo <- drugSets_compoundInfo[
    drugSets_compoundInfo$id %in% selectDrugSet, ]

test_that("Match ranked CMap perturbations to drug set identifiers", {
    # Match data based on PubChem CID (information found within compound info)
    sequence  <- 1:100
    data <- data.table("compound_perturbation"=paste("sig_id", sequence),
                       "spearman_rank"=sequence, "GSEA_rank"=sequence,
                       "rankProduct_rank"=sequence)
    # Metadata
    metadata <- data.table("sig_id"=paste("sig_id", sequence),
                           "pert_id"=paste("BRD", sequence),
                           "pert_iname"=paste("compound", sequence))
    attr(data, "metadata") <- metadata

    # Compound information
    compoundInfo <- data.table("pert_iname"=paste("compound", sequence),
                               "moa"=paste("moa", sequence),
                               "target"=paste("target", sequence),
                               "pubchem_cid"=paste("cid", sequence))
    compoundInfo[seq(nrow(drugSets_compoundInfo)), ][["pubchem_cid"]] <-
        as.character(drugSets_compoundInfo[["PubChem CID"]])
    attr(data, "compoundInfo") <- compoundInfo

    class(data) <- c("similarPerturbations", "referenceComparison", class(data))
    matched <- matchStatsWithDrugSetsID(drugSets, data)
    expect_identical(attr(matched, "valuesCol"), "rankProduct_rank")

    # Expected result
    expected <- setNames(data[["rankProduct_rank"]],
                         data[["compound_perturbation"]])
    names(expected) <- metadata[["pert_iname"]][match(
        names(expected), metadata[["sig_id"]])]
    names(expected) <- compoundInfo[["pubchem_cid"]][match(
        names(expected), compoundInfo[["pert_iname"]])]

    matches <- lapply(paste0("^", stripStr(names(expected)), "$"), grep,
                      stripStr(drugSets_compoundInfo$`PubChem CID`))
    matchesRep <- vapply(matches, length, numeric(1))

    matches[matchesRep == 0] <- NA
    valuesRep <- matchesRep
    valuesRep[valuesRep == 0] <- 1
    expected <- rep(expected, valuesRep)
    names(expected) <- drugSets_compoundInfo[["id"]][
        as.numeric(unlist(matches))]
    names(expected)[is.na(names(expected))] <- rep(
        data[["compound_perturbation"]], valuesRep)[is.na(names(expected))]

    expect_identical(sort(matched), sort(expected))

    # Match data based on name (information found within metadata)
    metadata[seq(nrow(drugSets_compoundInfo)), "pert_iname"] <-
        drugSets_compoundInfo["name"]
    attr(data, "metadata")     <- metadata
    compoundInfo[seq(nrow(drugSets_compoundInfo)), "pert_iname"] <-
        drugSets_compoundInfo["name"]
    attr(data, "compoundInfo") <- compoundInfo

    expect_identical(sort(matched), sort(expected))
})

test_that("Match putative targeting drugs to drug set identifiers", {
    sequence  <- 1:100
    predicted <- data.table("compound"=paste("id", sequence),
                            "spearman_rank"=sequence, "GSEA_rank"=sequence,
                            "rankProduct_rank"=sequence)
    compoundInfo <- data.table(id=paste("id", sequence),
                               name=paste("compound", sequence))
    compoundInfo[seq(nrow(selectDrugSet_compoundInfo)), "name"] <-
        selectDrugSet_compoundInfo["name"]

    predicted[ , "compound"] <- make.unique(compoundInfo[["id"]])
    attr(predicted, "compoundInfo") <- compoundInfo
    class(predicted) <- c("targetingDrugs", "referenceComparison",
                          class(predicted))

    matched <- matchStatsWithDrugSetsID(drugSets, predicted)
    expect_identical(attr(matched, "valuesCol"), "rankProduct_rank")

    # Use compound names (if matches more than one identifier, return them all)
    res <- setNames(predicted[["rankProduct_rank"]], predicted[["compound"]])
    compoundNames <- compoundInfo$name[match(names(res), compoundInfo$id)]

    # # Avoid repeating identifiers
    # res <- res[!duplicated(compoundNames)]
    # compoundNames <- compoundNames[!duplicated(compoundNames)]
    #
    # matches <- lapply(stripStr(compoundNames), grep,
    #                   stripStr(selectDrugSet_compoundInfo[["name"]]))
    # matchesRep <- vapply(matches, length, numeric(1))
    # matches[matchesRep == 0] <- names(res)[matchesRep == 0]
    # matchedID <- selectDrugSet_compoundInfo[["id"]][
    #     suppressWarnings(as.numeric(unlist(matches)))]
    # matchedID[is.na(matchedID)] <- unlist(matches)[is.na(matchedID)]
    #
    # matchesRep[matchesRep == 0] <- 1
    # res <- rep(res, matchesRep)
    # names(res) <- matchedID
    #
    # matched <- matchStatsWithDrugSetsID(drugSets, predicted)
    # expect_is(matched, "integer")
    # expect_identical(sort(matched), sort(res))
})

customDrugStat <- runif(1000, -5, 5)

test_that("Match a named numeric vector to drug set identifiers", {
    # Use same identifiers between drug sets and stats
    customDrugStat        <- sort(customDrugStat)
    names(customDrugStat) <- selectDrugSet
    names(customDrugStat)[is.na(names(customDrugStat))] <-
        seq(sum(is.na(names(customDrugStat))))
    names(customDrugStat) <- make.unique(names(customDrugStat))

    matched <- matchStatsWithDrugSetsID(drugSets, customDrugStat)
    expect_identical(sort(matched), sort(customDrugStat))
    expect_identical(attr(matched, "valuesCol"), "values")

    # Use different identifiers between drug sets and stats
    compoundInfo <- attr(drugSets, "compoundInfo")
    compoundIDs  <- setNames(compoundInfo[["PubChem SID"]],
                             compoundInfo[["id"]])

    customDrugStat2 <- customDrugStat
    names(customDrugStat2) <- compoundIDs[names(customDrugStat2)]
    names(customDrugStat2)[is.na(names(customDrugStat2))] <- names(
        customDrugStat)[is.na(names(customDrugStat2))]

    matched <- matchStatsWithDrugSetsID(drugSets, customDrugStat2)
    expect_identical(sort(matched), sort(customDrugStat))
    expect_identical(attr(matched, "valuesCol"), "values")
})

test_that("Match a custom list of drug sets to drug set identifiers", {
    names(customDrugStat) <- seq(customDrugStat)
    compounds  <- names(sort(customDrugStat))
    customSets <- list("head"=head(compounds, 20), "tail"=tail(compounds, 20),
                       "random"=sample(compounds, 20))
    matched <- matchStatsWithDrugSetsID(customSets, customDrugStat)
    expect_identical(sort(matched), sort(customDrugStat))
    expect_identical(attr(matched, "valuesCol"), "values")
})

# Analyse drug set enrichment --------------------------------------------------

data("diffExprStat")

test_that("Analyse drug set enrichment for ranked CMap perturbations", {
    data("cmapPerturbationsCompounds")
    perturbations <- cmapPerturbationsCompounds
    data <- rankSimilarPerturbations(diffExprStat, perturbations,
                                     method="spearman")
    # data <- readRDS("../cTRAP_projects/CMapResults_cTRAPvignette.rds")
    dsea <- analyseDrugSetEnrichment(drugSets, data)
    expect_is(dsea, "data.table")
    expect_identical(dsea$padj, sort(dsea$padj))
    expect_identical(colnames(dsea), c("pathway", "pval", "padj", "ES", "NES",
                                       "nMoreExtreme", "size", "leadingEdge"))
    expect_equal(head(dsea$pathway),
                 c("FDA_status: Clinical trial", "Mutagenic: none",
                   "Tumorigenic: high", "Tumorigenic: none",
                   "Reproductive Effective: high",
                   "Reproductive Effective: none"))
})

test_that("Analyse drug set enrichment for putative targeting drugs", {
    gdsc      <- loadExpressionDrugSensitivityAssociation("GDSC 7")
    predicted <- predictTargetingDrugs(diffExprStat, gdsc)
    dsea      <- analyseDrugSetEnrichment(drugSets, predicted)
    expect_is(dsea, "data.table")
    expect_identical(dsea$padj, sort(dsea$padj))
    expect_identical(colnames(dsea), c("pathway", "pval", "padj", "ES", "NES",
                                       "nMoreExtreme", "size", "leadingEdge"))
    expect_true(all(c("Aromatic Nitrogens: 6", "FDA_status: Clinical trial",
                      "Reproductive Effective: high", "Irritant: high",
                      "Aromatic Nitrogens: 2", "Amines: 0") %in% dsea$pathway))
})

customDrugStat <- runif(100, -5, 5)
names(customDrugStat) <- seq(customDrugStat)

test_that("Analyse drug set enrichment for a named numeric vector", {
    dsea <- analyseDrugSetEnrichment(drugSets, customDrugStat)
    expect_is(dsea, "data.table")
    expect_identical(dsea$padj, sort(dsea$padj))
    expect_identical(colnames(dsea), c("pathway", "pval", "padj", "ES", "NES",
                                       "nMoreExtreme", "size", "leadingEdge"))
    expect_true(all(c("Mutagenic: high", "Mutagenic: none", "Tumorigenic: none",
                      "Reproductive Effective: none", "Irritant: none",
                      "Fragments: 2") %in% dsea$pathway))
})

test_that("Analyse drug set enrichment using a custom list of drug sets", {
    compounds  <- names(sort(customDrugStat))
    customSets <- list("head"=head(compounds, 20), "tail"=tail(compounds, 20),
                       "random"=sample(compounds, 20))
    dsea <- analyseDrugSetEnrichment(customSets, customDrugStat)
    expect_is(dsea, "data.table")
    expect_identical(dsea$padj, sort(dsea$padj))
    expect_identical(colnames(dsea), c("pathway", "pval", "padj", "ES", "NES",
                                       "nMoreExtreme", "size", "leadingEdge"))
    expect_equal(sort(dsea$pathway), sort(names(customSets)))
})

