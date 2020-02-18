# Internal GSEA plot -----------------------------------------------------------

#' Perform GSEA
#'
#' @importFrom fgsea calcGseaStat
#'
#' @keywords internal
#' @return List with results of running GSEA
performGSEA <- function(pathways, stats) {
    pathway <- unname(as.vector(na.omit(match(pathways, names(stats)))))
    pathway <- sort(pathway)
    gseaStat <- calcGseaStat(stats, selectedStats=pathway,
                             returnAllExtremes=TRUE)
    enrichment <- data.frame(
        rank =c(0, c(pathway - 1, pathway), length(stats) + 1),
        score=c(0, c(gseaStat$bottoms, gseaStat$tops), 0))
    res <- list(pathway=pathway, enrichmentScore=enrichment, stat=gseaStat)
    return(res)
}

#' Render GSEA enrichment plot
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_line annotate
#' scale_x_continuous scale_y_continuous theme theme_bw ggtitle labs
#' element_text element_blank element_rect expand_scale
#'
#' @keywords internal
#' @return GSEA enrichment plot
plotESplot <- function(enrichmentScore, gseaStat) {
    amp <- range(enrichmentScore$score)
    ES  <- amp[which.max(abs(amp))]
    enrichmentPlot <- ggplot(enrichmentScore, aes(x=enrichmentScore$rank,
                                                  y=enrichmentScore$score)) +
        geom_rug(alpha=0.2, sides="b", length = unit(0.1, "npc")) +
        geom_line(colour="orange", na.rm=TRUE, size=0.7) +
        geom_hline(yintercept=0, colour="darkgrey", linetype="longdash") +
        geom_hline(yintercept=ES, colour="#3895D3") +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=expand_scale(mult = c(0.2, 0.1))) +
                           # breaks=c(round(ES, 2), seq(-20, 20, .1))) +
        labs(y="Enrichment score") +
        theme_bw() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              panel.grid.major.x=element_blank(),
              panel.grid.minor=element_blank(),
              axis.text=element_text(size=10),
              axis.title=element_text(size=12))

    if (ES > 0) {
        label_x     <- max(enrichmentScore$rank)
        label_hjust <- 1.2
        label_vjust <- 1.5
    } else {
        label_x     <- 0
        label_hjust <- -0.2
        label_vjust <- -0.5
    }
    enrichmentPlot <- enrichmentPlot +
        annotate("text", y=ES, colour="#3895D3", label=round(ES, 3),
                 x=label_x, hjust=label_hjust, vjust=label_vjust)
    return(enrichmentPlot)
}

#' Plot metric distribution
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous theme
#' theme_bw labs guides element_text geom_area scale_fill_gradient2
#'
#' @keywords internal
#' @return Metric distribution plot
plotMetricDistribution <- function(stat) {
    # Scale number of breaks according to number of ranked elements
    breaks       <- max(round(-120 + 55 * log10(length(stat))), 10)
    quantile     <- cut(stat, breaks=breaks, labels=FALSE)
    quantile     <- seq(min(stat), max(stat), length.out=breaks)[quantile]
    rankedMetric <- data.frame(sort=seq(stat), stat=stat, quantile=quantile)

    metricPlot <- ggplot(rankedMetric, aes_string("sort", "stat")) +
        scale_x_continuous(expand=c(0, 0)) +
        scale_y_continuous(expand=c(0, 0)) +
        labs(x="Rank", y="Ranked metric") +
        guides(colour=FALSE, fill=FALSE) +
        geom_area(aes_string(group="quantile", fill="quantile"), na.rm=TRUE,
                  position="identity") +
        scale_fill_gradient2(low="dodgerblue", mid="grey95", high="orangered",
                             midpoint=0) +
        theme_bw() +
        theme(plot.margin=unit(c(0, 0, 5.5, 0), "pt"),
              panel.grid.major.x=element_blank(),
              panel.grid.minor=element_blank(),
              axis.text=element_text(size=10),
              axis.title=element_text(size=12))
    return(metricPlot)
}

#' Plot gene set enrichment analysis (GSEA)
#'
#' @param stats Named numeric vector: statistics
#' @param title Character: title
#' @inheritParams plot.referenceComparison
#' @param gseaParam Numeric: GSEA-like parameter
#'
#' @importFrom ggplot2 ggtitle theme unit
#' @importFrom cowplot plot_grid ggdraw draw_label
#' @importFrom stats na.omit
#'
#' @return Grid of plots illustrating a GSEA plot
#' @keywords internal
plotGSEA <- function(stats, geneset, genes=c("both", "top", "bottom"),
                     title="GSEA plot", gseaParam=1) {
    genes <- match.arg(genes)
    if (is.list(geneset)) {
        topGenes    <- c(geneset[["top"]], geneset[["custom"]])
        bottomGenes <- geneset[["bottom"]]
    } else {
        topGenes    <- geneset
        bottomGenes <- NULL
    }
    if (!any(genes %in% c("top", "both")))    topGenes    <- NULL
    if (!any(genes %in% c("bottom", "both"))) bottomGenes <- NULL

    statsOrd <- sort(stats, decreasing=TRUE)
    statsAdj <- abs(statsOrd ^ gseaParam)
    statsAdj <- sign(statsOrd) * statsAdj / max(statsAdj)

    plotOneESplot <- function(end, genes, stats) {
        if (is.null(genes)) return(NULL)
        gseaRes        <- suppressWarnings(performGSEA(genes, stats))
        enrichmentPlot <- plotESplot(gseaRes$enrichmentScore, gseaRes$stat)
        return(enrichmentPlot)
    }
    title        <- ggdraw() + draw_label(title, size=14)
    topESplot    <- plotOneESplot("top",    topGenes,    statsAdj)
    bottomESplot <- plotOneESplot("bottom", bottomGenes, statsAdj)
    metricPlot   <- plotMetricDistribution(statsOrd)

    plots <- list(topESplot, bottomESplot, metricPlot)
    plots <- Filter(length, plots)
    plotHeights <- c(1, 6,
                     if (!is.null(topGenes))    6,
                     if (!is.null(bottomGenes)) 6)

    grid <- plot_grid(title, plotlist=plots, align="v", ncol=1,
                      rel_heights=plotHeights)
    return(grid)
}

# Internal scatter plot for correlations ---------------------------------------

#' Render scatter plot to show a single relationship
#'
#' @param perturbation List of named numeric vectors containing the differential
#' expression profile score per gene for a perturbation; each perturbation of
#' the list will be rendered with a different colour
#' @param ylabel Character: Y axis label
#' @param diffExprGenes Named numeric vector
#'
#' @importFrom ggplot2 ggplot geom_point geom_rug geom_abline xlab ylab
#' geom_density_2d theme guides guide_legend theme_bw
#'
#' @keywords internal
#' @return Scatter plot
plotSingleCorr <- function(perturbation, ylabel, diffExprGenes) {
    prepareDFperPert <- function(perturbation, diffExprGenes) {
        # Intersect common genes
        genes         <- intersect(names(perturbation), names(diffExprGenes))
        perturbation  <- perturbation[genes]
        diffExprGenes <- diffExprGenes[genes]

        df <- data.frame(diffExprGenes, zscores=perturbation)
        return(df)
    }
    dfPerPert  <- lapply(perturbation, prepareDFperPert, diffExprGenes)
    dfAllPerts <- bind_rows(dfPerPert)
    dfAllPerts <- cbind(dfAllPerts, perturbation=rep(names(dfPerPert),
                                                     lapply(dfPerPert, nrow)))
    dfAllPerts$perturbation <- parseCMapID(dfAllPerts$perturbation,
                                           cellLine=TRUE)

    multipleCellLines <- length(perturbation) > 1
    if (multipleCellLines) {
        aesMap <- aes_string("diffExprGenes", "zscores", colour="perturbation")
        guide  <- guides(colour=guide_legend(title="Cell line"))
    } else {
        aesMap <- aes_string("diffExprGenes", "zscores")
        guide  <- NULL
    }

    plot <- ggplot(dfAllPerts, aesMap) +
        geom_point(alpha=0.1) +
        geom_rug(alpha=0.1) +
        geom_density_2d() +
        xlab("Differentially expressed genes (input)") +
        ylab(ylabel) +
        theme_bw() +
        theme(legend.position="bottom") +
        guide
    return(plot)
}

# Exported plot functions ------------------------------------------------------

prepareLabel <- function(data) {
    prepareLabel_similarPerturbations <- function(k, pert) {
        info <- print(pert, pert[[1]][[k]])
        collapse <- function(var) paste(unique(var), collapse="/")

        name <- collapse(info$metadata$pert_iname)
        type <- collapse(info$metadata$pert_type)

        isOverexpression <- startsWith(type, "trt_oe")
        isLossOfFunction <- startsWith(type, "trt_sh") ||
            startsWith(type, "trt_xpr")

        if (isOverexpression) {
            name <- paste(name, "OE")
        } else if (isLossOfFunction) {
            name <- paste(name, "KD")
        }

        cellLine <- unique(info$metadata$cell_id)
        if (length(unique(cellLine)) > 1) cellLine <- "mean"

        dose       <- collapse(info$metadata$pert_idose)
        timepoint  <- collapse(info$metadata$pert_itime)
        res <- sprintf("%s (%s, %s at %s)", name, cellLine, dose, timepoint)
        return(res)
    }

    prepareLabel_targetingDrugs <- function(k, data) {
        compoundInfo    <- attr(data, "compoundInfo")
        data$compound   <- as.character(data$compound)
        compoundInfo$id <- as.character(compoundInfo$id)

        merged          <- merge(data, compoundInfo, by.x="compound", by.y="id")
        compound        <- merged[k]
        name            <- compound[["name"]]
        if (is.null(name) || is.na(name) || name == "") {
            name <- compound[["compound"]]
            if (is.na(name)) name <- "NA"
        }
        target <- compound[["target"]]
        target <- gsub(", |;", "/", target)
        if (is.null(target) || target == "") target <- "?"

        targetPathway <- compound[["target pathway"]]
        if (is.null(targetPathway)) {
            targetPathway <- NA
            res <- sprintf("%s (%s)", name, target)
        } else {
            res <- sprintf("%s (%s: %s)", name, target, targetPathway)
        }
        return(res)
    }

    if (is(data, "similarPerturbations")) {
        FUN <- prepareLabel_similarPerturbations
    } else if (is(data, "targetingDrugs")) {
        FUN <- prepareLabel_targetingDrugs
    }

    res <- sapply(seq(nrow(data)), FUN, data)
    return(res)
}

plotComparison <- function(x, method, n, showMetadata,
                           plotNonRankedPerturbations, alpha) {
    if (method == "gsea") {
        stat     <- "GSEA"
        statRank <- "GSEA_rank"
        yLabel   <- "Weighted Connectivity Score (WCTS)"
    } else if (method == "rankProduct") {
        stat     <- "rankProduct_rank"
        statRank <- stat
        yLabel   <- "Rank product's rank"
    } else {
        stat     <- paste0(method, "_coef")
        statRank <- paste0(method, "_rank")
        yLabel   <- sprintf("%s's correlation coefficient", capitalize(method))
    }
    if (!stat %in% colnames(x)) {
        stop(sprintf("'%s' was not run with method '%s'",
                     deparse(substitute(x)), method))
    }

    # Label perturbations
    if (!plotNonRankedPerturbations) {
        x <- x[order(x[[statRank]], na.last=NA)]
        x$ranked <- "Ranked"
    } else {
        x$ranked <- ifelse(!is.na(x[[statRank]]), "Ranked", "Non-ranked")
    }
    sortedPert <- order(x[[stat]], decreasing=TRUE)

    top <- bottom <- n
    if (length(n) == 2) {
        top    <- n[[1]]
        bottom <- n[[2]]
    }
    index   <- unique(c(head(sortedPert, top), tail(sortedPert, bottom)))
    x$label <- ""

    if (showMetadata && is(x, "referenceComparison")){
        x$label[index] <- prepareLabel(x[index])
    } else {
        x$label[index] <- x[[1]][index]
    }

    # Correlation coefficient with labels for top and bottom perturbations
    vars <- aes_string(x=1, y=stat, label="label", colour="ranked")
    rug  <- geom_rug(alpha=alpha, sides="l")

    plot <- ggplot(x, vars) +
        rug +
        scale_colour_manual(values=c("Ranked"="#56B4E9",
                                     "Non-ranked"="#999999")) +
        geom_text_repel(nudge_x=0.15, direction="y", hjust=0, segment.size=0.2,
                        show.legend=FALSE, size=6) +
        xlim(1, 1.375) +
        ylab(yLabel) +
        theme_classic(base_size=18) +
        theme(legend.title=element_blank(), legend.position="bottom",
              axis.line.x=element_blank(), axis.ticks.x=element_blank(),
              axis.text.x=element_blank(), axis.title.x=element_blank())

    if (length(unique(x$ranked)) == 1) plot <- plot + guides(colour=FALSE)
    if (method == "rankProduct") plot <- plot + scale_y_reverse()
    return(plot)
}

#' Plot data comparison
#'
#' If \code{element = NULL}, comparison is plotted based on all elements.
#' Otherwise, show scatter or GSEA plots for a single element compared with
#' previously given differential expression results.
#'
#' @param x \code{referenceComparison} object: obtained after running
#'   \code{\link{rankSimilarPerturbations}()} or
#'   \code{\link{predictTargetingDrugs}()}
#' @param element Character: identifier in the first column of \code{x}
#' @param method Character: method to plot results; choose between
#'   \code{spearman}, \code{pearson}, \code{gsea} or \code{rankProduct} (the
#'   last one is only available if \code{element = NULL})
#' @param n Numeric: number of top and bottom genes to label (if a vector of two
#'   numbers is given, the first and second numbers will be used as the number
#'   of top and bottom genes to label, respectively); only used if
#'   \code{element = NULL}
#' @param showMetadata Boolean: show available metadata information instead of
#'   identifiers (if available)? Only used if \code{element = NULL}
#' @param alpha Numeric: transparency; only used if \code{element = NULL}
#' @param plotNonRankedPerturbations Boolean: plot non-ranked data in grey? Only
#'   used if \code{element = NULL}
#' @inheritParams compareAgainstReferencePerMethod
#' @param genes Character: when plotting gene set enrichment analysis (GSEA),
#'   plot most up-regulated genes (\code{genes = "top"}), most down-regulated
#'   genes (\code{genes = "bottom"}) or both (\code{genes = "both"}); only used
#'   if \code{method = "gsea"} and \code{geneset = NULL}
#'
#' @param ... Extra arguments currently not used
#'
#' @importFrom graphics plot
#' @importFrom R.utils capitalize
#' @importFrom ggplot2 ggplot aes_string geom_point geom_hline ylab theme
#'   element_blank scale_colour_manual xlim theme_classic guides scale_y_reverse
#' @importFrom ggrepel geom_text_repel
#'
#' @family functions related with the ranking of CMap perturbations
#' @family functions related with the prediction of targeting drugs
#' @return Plot illustrating the reference comparison
#' @export
#'
#' @examples
#' # Example of a differential expression profile
#' data("diffExprStat")
#'
#' \dontrun{
#' # Download and load CMap perturbations to compare with
#' cellLine <- "HepG2"
#' cmapMetadataKD <- filterCMapMetadata(
#'   "cmapMetadata.txt", cellLine=cellLine,
#'   perturbationType="Consensus signature from shRNAs targeting the same gene")
#'
#' cmapPerturbationsKD <- prepareCMapPerturbations(
#'   cmapMetadataKD, "cmapZscores.gctx", "cmapGeneInfo.txt", loadZscores=TRUE)
#' }
#'
#' # Rank similar CMap perturbations
#' compareKD <- rankSimilarPerturbations(diffExprStat, cmapPerturbationsKD)
#'
#' plot(compareKD, method="spearman", n=c(7, 3))
#' plot(compareKD, method="pearson")
#' plot(compareKD, method="gsea")
#'
#' plot(compareKD, compareKD[[1, 1]])
plot.referenceComparison <- function(x, element=NULL,
                                     method=c("spearman", "pearson", "gsea",
                                              "rankProduct"),
                                     n=c(3, 3), showMetadata=TRUE,
                                     plotNonRankedPerturbations=FALSE,
                                     alpha=0.3, geneSize=150,
                                     genes=c("both", "top", "bottom"), ...) {
    if (!is.null(element)) {
        if (length(element) > 1) {
            stop("argument 'element' should be a character of length one")
        } else if (element %in% method) {
            # Legacy: 'method' as the second argument
            method  <- element
            element <- NULL
        }
    }
    cols   <- tolower(colnames(x))
    method <- method[tolower(method) %in% gsub("_rank", "", cols)]
    if (length(method) == 0) stop("no supported methods selected")
    method <- method[[1]]
    method <- match.arg(method)

    if (is.null(element)) {
        plot <- plotComparison(
            x, method=method, n=n, showMetadata=showMetadata,
            plotNonRankedPerturbations=plotNonRankedPerturbations, alpha=alpha)
    } else if (is(x, "similarPerturbations")) {
        metadata     <- attr(x, "metadata")
        zscores      <- attr(x, "zscoresFilename")
        geneInfo     <- attr(x, "geneInfo")
        compoundInfo <- attr(x, "compoundInfo")
        cmapPerturbations <- suppressMessages(
            prepareCMapPerturbations(metadata, zscores, geneInfo, compoundInfo))
        input <- c(attr(compareKD, "diffExprGenes"), # legacy
                   attr(compareKD, "input"))
        plot  <- plot(cmapPerturbations, element, input=input, method=method,
                      geneSize=geneSize, genes=genes, ...)
    } else if (is(x, "targetingDrugs")) {
        plot <- plotTargetingDrug(x, element, method=method, geneSize=geneSize,
                                  genes=genes, ...)
    }
    return(plot)
}

#' Compare vector against its quantile
#'
#' Check which elements of the vector are lower/greater than or equal to the
#' quantile of a given vector.
#'
#' @param vec Numeric vector
#' @param prob Numeric: probability value between [0,1] to produce sample
#'   quantiles
#' @param lte Boolean: check if values are <= quantile? If \code{FALSE}, checks
#'   if values are >= quantile
#'
#' @importFrom stats quantile
#'
#' @return Boolean vector regarding compared elements
#' @keywords internal
compareQuantile <- function(vec, prob, lte=FALSE) {
    if (lte) {
        cmp  <- `<=`
    } else {
        cmp  <- `>=`
        prob <- 1 - prob
    }
    res <- cmp(vec, quantile(vec, prob))
    return(res)
}

#' Plot similar perturbations against predicted targeting drugs
#'
#' @param targetingDrugs \code{targetingDrugs} object
#' @param similarPerturbations \code{similarPerturbations} object
#' @param column Character: column to plot (must be available in both databases)
#' @param labelBy Character: column in \code{similarPerturbations}, its metadata
#'   or compound information to be used for labelling
#' @param showAllScores Boolean: show all scores? If \code{FALSE}, only the best
#'   score per compound will be plotted
#' @param quantileThreshold Numeric: quantile to use for highlight values within
#'   [0, 1]
#'
#' @importFrom ggplot2 ggplot geom_point xlab ylab theme_bw
#' @importFrom ggrepel geom_text_repel
#' @importFrom data.table data.table
#'
#' @family functions related with the ranking of CMap perturbations
#' @family functions related with the prediction of targeting drugs
#' @return \code{ggplot2} plot
#' @export
plotTargetingDrugsVSsimilarPerturbations <- function(
    targetingDrugs, similarPerturbations, column, labelBy="pert_iname",
    quantileThreshold=0.25, showAllScores=FALSE) {

    if (!column %in% colnames(targetingDrugs) &&
        !column %in% colnames(similarPerturbations)) {
        error <- paste0("Selected column ('%s') not available in the columns",
                        "names of provided datasets")
        stop(sprintf(error, column))
    }

    metadata     <- attr(similarPerturbations, "metadata")
    compoundInfo <- attr(similarPerturbations, "compoundInfo")
    merged <- merge(similarPerturbations, metadata,
                    by.x=colnames(similarPerturbations)[[1]],
                    by.y=colnames(metadata)[[1]], all.x=TRUE)
    merged <- merge(merged, compoundInfo, by="pert_iname", all.x=TRUE)

    # Check for intersecting compounds across identifying columns
    keys <- c("pert_iname", "pert_id", "smiles", "InChIKey", "pubchem_cid")
    commonValues <- lapply(keys, function(k, data, comp) {
        values <- stripStr(data[[k]])
        return(data[[k]][which(values %in% comp)])
    }, data=merged, comp=stripStr(targetingDrugs[[1]]))
    mostCommonValuesKey <- which.max(vapply(commonValues, length, numeric(1)))
    commonCompounds <- unique(commonValues[[mostCommonValuesKey]])
    keyCol <- keys[[mostCommonValuesKey]]

    # Prepare data to plot
    yData <- merged[stripStr(merged[[keyCol]]) %in% stripStr(commonCompounds)]
    yAxis <- setNames(yData[[column]],       yData[[keyCol]])
    id    <- setNames(yData[["pert_iname"]], yData[[keyCol]])

    if (!is.null(labelBy)) {
        label <- setNames(yData[[labelBy]], yData[[keyCol]])
    } else {
        label <- NULL
    }

    xData <- targetingDrugs[match(stripStr(commonCompounds),
                                  stripStr(targetingDrugs[[1]]))]
    xAxis <- setNames(xData[[column]], stripStr(xData[[1]]))
    xAxis <- xAxis[stripStr(names(yAxis))]

    # Discard missing values
    df <- data.table(id=id, x=xAxis, y=yAxis, label=label)
    df <- df[!is.na(df$x) & !is.na(df$y), ]

    # Only show the best (instead of all) scores
    isRank <- endsWith(column, "rank")
    if (!showAllScores) {
        df <- df[order(df$y, decreasing=isRank)]
        df <- df[unique(match(df$id, df$id))]
    }

    df$highlight <- compareQuantile(df$x, quantileThreshold, lte=TRUE) &
        compareQuantile(df$y, quantileThreshold, lte=!isRank)
    if (!is.null(label)) df$label[!df$highlight] <- NA

    xlabel <- "Gene expression and drug sensitivity association"
    source <- attr(targetingDrugs, "expressionDrugSensitivityCor")$source
    if (!is.null(source)) xlabel <- sprintf("%s (%s)", xlabel, source)
    xlabel <- sprintf("%s: %s", xlabel, column)
    ylabel <- paste("CMap comparison:", column)

    plot <- ggplot(df, aes_string("x", "y")) +
        geom_point(data=df[df$highlight],  colour="orange", show.legend=FALSE) +
        geom_point(data=df[!df$highlight], colour="grey",   show.legend=FALSE,
                   alpha=0.7) +
        xlab(xlabel) +
        ylab(ylabel) +
        theme_bw()

    if (!is.null(label))
        plot <- plot + geom_text_repel(label=df$label, na.rm=TRUE)
    return(plot)
}
