#' Perform GSEA
#'
#' @importFrom fgsea calcGseaStat
#'
#' @keywords internal
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
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_line
#' scale_x_continuous scale_y_continuous theme theme_bw ggtitle labs
#' element_text element_blank element_rect
#'
#' @keywords internal
plotGSEAenrichment <- function(enrichmentScore, gseaStat, titleSize=14,
                               axisTitleSize=12, axisTextSize=10,
                               pointSize=0.1) {
    enrichmentPlot <- ggplot(enrichmentScore, aes(x=enrichmentScore$rank,
                                                  y=enrichmentScore$score)) +
        geom_point(colour="green", size=pointSize, na.rm=TRUE) +
        geom_hline(yintercept=0, colour="black", size=0.3) +
        geom_line(colour="green", size=1, na.rm=TRUE) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(limits=range(c(gseaStat$bottoms, gseaStat$tops)),
                           expand=c(0,0)) +
        theme_bw(base_size=18) +
        theme(panel.border=element_rect(colour="black"),
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_text(size=axisTitleSize, colour="black"),
              axis.text.y=element_text(size=axisTextSize, colour="black"),
              plot.title=element_text(size=titleSize, colour="black",
                                      hjust=0.5)) +
        labs(y="Enrichment score")

    # enrichmentPlot <- enrichmentPlot +
    #     geom_hline(yintercept=max(gseaStat$tops), colour="red",
    #                linetype="dashed")
    return(enrichmentPlot)
}

#' Render gene hit plot
#'
#' @importFrom ggplot2 ggplot geom_line aes geom_segment scale_x_continuous
#' scale_y_continuous theme theme_bw element_rect element_blank
#'
#' @keywords internal
plotGeneHits <- function(enrichmentScore, gseaStat, pathway) {
    geneHitsPlot <- ggplot(enrichmentScore, aes(x=enrichmentScore$rank,
                                                y=enrichmentScore$score)) +
        geom_line(colour=NA, na.rm=TRUE) +
        geom_segment(data=data.frame(pathway), size=0.4, na.rm=TRUE,
                     mapping=aes(
                         x=pathway, xend=pathway,
                         y=min(gseaStat$bottoms), yend=max(gseaStat$tops)))+
        scale_y_continuous(expand=c(0, 0)) +
        scale_x_continuous(expand=c(0, 0), breaks=seq(0, 50000, 2500)) +
        theme_bw() +
        theme(panel.border=element_rect(colour="black"),
              panel.grid=element_blank(), axis.title=element_blank(),
              axis.text=element_blank(), axis.ticks=element_blank())
    return(geneHitsPlot)
}

#' Plot metric distribution
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous theme
#' theme_bw labs guides element_text geom_area scale_fill_gradient2
#'
#' @keywords internal
plotMetricDistribution <- function(statsOrd, breaks=50, axisTitleSize=12,
                                   axisTextSize=10) {
    rankedMetric <- data.frame(sort=seq(statsOrd), stat=statsOrd)
    rankedMetric$quantile <- cut(rankedMetric$stat, breaks=breaks, labels=FALSE)
    rankedMetric$quantile <- seq(min(rankedMetric$stat), max(rankedMetric$stat),
                                 length.out=breaks)[rankedMetric$quantile]

    metricPlot <- ggplot(rankedMetric, aes(rankedMetric$sort,
                                           rankedMetric$stat)) +
        scale_x_continuous(expand=c(0, 0), breaks = seq(0, 50000, 2500)) +
        scale_y_continuous(expand=c(0, 0)) +
        theme_bw() +
        theme(panel.border=element_rect(colour="black"),
              axis.title=element_text(size=axisTitleSize, colour="black"),
              axis.text=element_text(size=axisTextSize, colour="black")) +
        labs(y="Ranked metric", x="Rank") +
        guides(colour=FALSE, fill=FALSE)

    # metricPlot <- metricPlot +
    #     geom_bar(stat="identity", aes(colour=y), size=0.1) +
    #     scale_colour_gradient2(low="dodgerblue", mid="white",
    #                            high="orangered", midpoint=0)
    metricPlot <- metricPlot +
        geom_area(position="identity",
                  aes_string(group="quantile", fill="quantile"), na.rm=TRUE) +
        scale_fill_gradient2(low="dodgerblue", mid="grey95", high="orangered",
                             midpoint=0)
    return(metricPlot)
}

#' Plot gene set enrichment analysis (GSEA)
#'
#' @inheritParams plotSingleCorr
#' @inheritParams fgsea::fgsea
#' @param title Character: plot title
#' @param titleSize Integer: plot title size
#' @param gseaParam Numeric: GSEA-like parameter
#'
#' @importFrom ggplot2 ggtitle theme unit
#' @importFrom cowplot plot_grid ggdraw draw_label
#' @importFrom stats na.omit
#'
#' @return Grid of plots illustrating a GSEA plot
#' @keywords internal
plotGSEA <- function(perturbation, pathways, genes=c("both", "top", "bottom"),
                     titleSize=14, axisTitleSize=12, axisTextSize=10,
                     enrichmentPlotPointSize=0.1, gseaParam=1) {
    if (length(perturbation) > 1) {
        stop("Plotting GSEA of multiple perturbations is currently unsupported")
    }
    stats <- perturbation[[1]]

    stats <- unclass(stats)
    statsOrd <- stats[order(rank(-stats))]
    statsAdj <- sign(statsOrd) * (abs(statsOrd) ^ gseaParam)
    statsAdj <- statsAdj / max(abs(statsAdj))

    prepareOneEndGseaPlot <- function(pathways, stats, ..., id="") {
        gseaRes <- performGSEA(pathways, stats)

        if (id != "") id <- sprintf(" for %s", id)
        message(sprintf("Preparing enrichment plot%s...", id))
        enrichmentPlot <- plotGSEAenrichment(gseaRes$enrichmentScore,
                                             gseaRes$stat, ...) +
            theme(plot.margin=unit(c(5.5, 5.5, 0, 5.5), "pt"))

        message(sprintf("Preparing gene hits plot%s...", id))
        geneHitsPlot <- plotGeneHits(gseaRes$enrichmentScore, gseaRes$stat,
                                     gseaRes$pathway) +
            theme(plot.margin=unit(c(0, 5.5, 0, 5.5), "pt"))
        return(list(enrichmentPlot=enrichmentPlot, geneHitsPlot=geneHitsPlot))
    }

    message("Preparing metric distribution plot...")
    metricPlot <- plotMetricDistribution(statsOrd, breaks=100,
                                         axisTitleSize=axisTitleSize,
                                         axisTextSize=axisTextSize) +
        theme(plot.margin=unit(c(0, 5.5, 5.5, 5.5), "pt"))

    if (genes %in% c("both", "top")) {
        topGenes <- pathways[[grep("top", names(pathways), value=TRUE)]]
        topPlot  <- prepareOneEndGseaPlot(topGenes, statsAdj, id="top genes")
    }

    if (genes %in% c("both", "bottom")) {
        bottomGenes <- pathways[[grep("bottom", names(pathways), value=TRUE)]]
        bottomPlot  <- prepareOneEndGseaPlot(bottomGenes, statsAdj,
                                             id="bottom genes")
    }

    title <- ggdraw() + draw_label("GSEA plot")
    if (genes == "both") {
        plotHeights <- c(1, 6, 1, 6, 1, 6)
        grid <- plot_grid(title, topPlot$enrichmentPlot,
                          topPlot$geneHitsPlot, bottomPlot$enrichmentPlot,
                          bottomPlot$geneHitsPlot, metricPlot,
                          align="v", ncol=1, rel_heights=plotHeights)
    } else if (genes == "top") {
        plotHeights <- c(1, 6, 1, 6)
        grid <- plot_grid(title, topPlot$enrichmentPlot,
                          topPlot$geneHitsPlot, metricPlot,
                          align="v", ncol=1, rel_heights=plotHeights)
    } else if (genes == "bottom") {
        plotHeights <- c(1, 6, 1, 6)
        grid <- plot_grid(title, bottomPlot$enrichmentPlot,
                          bottomPlot$geneHitsPlot, metricPlot,
                          align="v", ncol=1, rel_heights=plotHeights)
    }
    return(grid)
}

#' Render scatter plot to show a single relationship
#'
#' @param perturbation List of named numeric vectors containing the differential
#' expression profile score per gene for a perturbation; each perturbation of
#' the list will be rendered with a different colour
#' @param ylabel Character: Y axis label
#' @param diffExprGenes Named vector
#'
#' @importFrom ggplot2 ggplot geom_point geom_rug geom_abline xlab ylab theme_bw
#' geom_density_2d theme guides guide_legend
#'
#' @keywords internal
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

#' Plot CMap data comparison
#'
#' @param x \code{cmapComparison} object
#' @param ... Extra arguments currently not used
#' @param method Character: method to plot results (\code{spearman},
#'   \code{pearson} or \code{gsea})
#' @param n Numeric: number of top and bottom genes to label (if a vector of two
#'   numbers is given, the first and second numbers will be used as the number
#'   of top and bottom genes to label, respectively)
#' @param showMetadata Boolean: show available metadata information instead of
#'   perturbation identifiers?
#' @param alpha Numeric: transparency
#' @param plotNonRankedPerturbations Boolean: plot non-ranked perturbations in
#'   grey?
#'
#' @importFrom graphics plot
#' @importFrom R.utils capitalize
#' @importFrom ggplot2 ggplot aes_string geom_point geom_hline ylab theme
#'   element_blank scale_colour_manual xlim theme_classic guides
#' @importFrom ggrepel geom_text_repel
#'
#' @return Plot illustrating the comparison with CMap data
#' @export
#'
#' @examples
#' data("cmapPerturbationsKD")
#'
#' # Compare against CMap using Spearman's correlation, Pearson's correlation
#' # and gene set enrichment analysis (GSEA) with the top and bottom 150 genes
#' # as gene sets
#' compareKD <- compareAgainstCMap(diffExprStat, cmapPerturbationsKD)
#'
#' plot(compareKD, "spearman", c(7, 3))
#' plot(compareKD, "pearson")
#' plot(compareKD, "gsea")
plot.cmapComparison <- function(x, method=c("spearman", "pearson", "gsea"),
                                n=c(3, 3), showMetadata=TRUE,
                                plotNonRankedPerturbations=FALSE,
                                alpha=0.3, ...) {
    method <- match.arg(method)

    if (method == "gsea") {
        stat     <- "GSEA"
        statRank <- "gsea_rank"
        yLabel   <- "Weighted Connectivity Score (WCTS)"
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

    prepareLabel <- function(k, perturbations) {
        info <- print(perturbations, perturbations[[1]][[k]])

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

        # cellLine   <- collapse(info$metadata$cell_id)
        cellLine <- unique(info$metadata$cell_id)
        if (length(unique(cellLine)) > 1) cellLine <- "mean"

        dose       <- collapse(info$metadata$pert_idose)
        timepoint  <- collapse(info$metadata$pert_itime)
        res <- sprintf("%s (%s, %s at %s)", name, cellLine, dose, timepoint)
        return(res)
    }

    if (showMetadata){
        x$label[index] <- sapply(index, prepareLabel, x)
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
    return(plot)
}

#' Plot CMap data comparison
#'
#' @param x \code{cmapPerturbations} object
#' @param ... Extra arguments currently not used
#' @param perturbation Character (perturbation identifier) or a
#'   \code{cmapComparison} table (from which the respective perturbation
#'   identifiers are retrieved)
#' @inheritParams prepareGSEApathways
#' @param method Character: method to plot results (\code{spearman},
#'   \code{pearson} or \code{gsea})
#' @param genes Character: when plotting gene set enrichment analysis (GSEA),
#'   plot top genes (\code{genes = "top"}), bottom genes
#'   (\code{genes = "bottom"}) or both (\code{genes = "both"}); only used if
#'   \code{method = "gsea"}
#'
#' @importFrom methods is
#' @importFrom stats setNames
#'
#' @export
#' @examples
#' data("diffExprStat")
#' data("cmapPerturbationsKD")
#'
#' compareKD <- compareAgainstCMap(diffExprStat, cmapPerturbationsKD)
#' EIF4G1knockdown <- grep("EIF4G1", compareKD[[1]], value=TRUE)
#' plot(cmapPerturbationsKD, EIF4G1knockdown, diffExprStat, method="spearman")
#' plot(cmapPerturbationsKD, EIF4G1knockdown, diffExprStat, method="pearson")
#' plot(cmapPerturbationsKD, EIF4G1knockdown, diffExprStat, method="gsea")
#'
#' data("cmapPerturbationsCompounds")
#' pert <- "CVD001_HEPG2_24H:BRD-A14014306-001-01-1:4.1"
#' plot(cmapPerturbationsCompounds, pert, diffExprStat, method="spearman")
#' plot(cmapPerturbationsCompounds, pert, diffExprStat, method="pearson")
#' plot(cmapPerturbationsCompounds, pert, diffExprStat, method="gsea")
#'
#' # Multiple cell line perturbations
#' pert <- "CVD001_24H:BRD-A14014306-001-01-1:4.1"
#' plot(cmapPerturbationsCompounds, pert, diffExprStat, method="spearman")
#' plot(cmapPerturbationsCompounds, pert, diffExprStat, method="pearson")
#'
#' # Currently unsupported!
#' # plot(cmapPerturbationsCompounds, pert, diffExprStat, method="gsea")
plot.cmapPerturbations <- function(x, perturbation, diffExprGenes,
                                   method=c("spearman", "pearson", "gsea"),
                                   geneSize=150,
                                   genes=c("both", "top", "bottom"), ...) {
    method <- match.arg(method)

    if (is(perturbation, "cmapComparison")) perturbation <- perturbation[[1]]

    cellLinePerts <- colnames(x)[
        parseCMapID(colnames(x), cellLine=FALSE) %in% perturbation]
    isSummaryPert <- length(cellLinePerts) > 0

    if (length(perturbation) == 0) {
        stop("One perturbation ID must be provided")
    } else if (length(perturbation) > 1) {
        stop("Only one perturbation ID is currently supported")
    } else if (!perturbation %in% colnames(x) && !isSummaryPert) {
        stop("Perturbation not found in the columns of the given dataset")
    }

    if (!isSummaryPert) cellLinePerts <- perturbation
    names(cellLinePerts) <- cellLinePerts
    if (is.character(x)) {
        zscores <- loadCMapZscores(x[cellLinePerts])
    } else {
        zscores <- unclass(x)
    }

    data <- lapply(cellLinePerts, function(pert, zscores) {
        sub <- zscores[ , pert, drop=FALSE]
        setNames(as.numeric(sub), rownames(sub))
    }, zscores)

    if (method != "gsea") {
        plotSingleCorr(data, perturbation, diffExprGenes)
    } else {
        pathways <- prepareGSEApathways(diffExprGenes=diffExprGenes,
                                        geneSize=geneSize)
        genes <- match.arg(genes)
        plotGSEA(data, pathways, genes)
    }
}
