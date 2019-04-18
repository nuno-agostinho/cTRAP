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
        labs(y="Enrichment score (ES)")

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
        geom_area(position="identity", aes(group=quantile, fill=quantile),
                  na.rm=TRUE) +
        scale_fill_gradient2(low="dodgerblue", mid="grey95", high="orangered",
                             midpoint=0)
    return(metricPlot)
}

#' Plot gene set enrichment analysis (GSEA)
#'
#' @inheritParams fgsea::fgsea
#' @param title Character: plot title
#' @param titleSize Integer: plot title size
#' @param gseaParam Numeric: GSEA-like parameter
#'
#' @importFrom ggplot2 ggtitle theme unit
#' @importFrom cowplot plot_grid
#' @importFrom stats na.omit
#'
#' @return Grid of plots illustrating a GSEA plot
#' @keywords internal
plotGSEA <- function(pathways, stats, genes=c("both", "top", "bottom"),
                     titleSize=14, axisTitleSize=12, axisTextSize=10,
                     enrichmentPlotPointSize=0.1, gseaParam=1) {
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

    if (genes == "both") {
        plotHeights <- c(6, 1, 6, 1, 6)
        grid <- plot_grid(topPlot$enrichmentPlot + ggtitle("GSEA plot"),
                          topPlot$geneHitsPlot, bottomPlot$enrichmentPlot,
                          bottomPlot$geneHitsPlot, metricPlot,
                          align="v", ncol=1, rel_heights=plotHeights)
    } else if (genes == "top") {
        plotHeights <- c(6, 1, 6)
        grid <- plot_grid(topPlot$enrichmentPlot + ggtitle("GSEA plot"),
                          topPlot$geneHitsPlot, metricPlot,
                          align="v", ncol=1, rel_heights=plotHeights)
    } else if (genes == "bottom") {
        plotHeights <- c(6, 1, 6)
        grid <- plot_grid(bottomPlot$enrichmentPlot + ggtitle("GSEA plot"),
                          bottomPlot$geneHitsPlot, metricPlot,
                          align="v", ncol=1, rel_heights=plotHeights)
    }
    return(grid)
}

#' Render scatter plot to show a single relationship
#'
#' @importFrom ggplot2 ggplot geom_point geom_rug geom_abline xlab ylab theme_bw
#' geom_density_2d
#'
#' @keywords internal
plotSingleCorr <- function(diffExprGenes, perturbation, perturbationID) {
    # Intersect common genes
    genes <- intersect(names(perturbation), names(diffExprGenes))
    perturbation  <- perturbation[genes]
    diffExprGenes <- diffExprGenes[genes]

    df <- data.frame(diffExprGenes, perturbation)
    plot <- ggplot(df, aes(diffExprGenes, perturbation)) +
        geom_point(alpha=0.1) +
        geom_rug(alpha=0.1) +
        geom_density_2d() +
        xlab("Differentially expressed genes (input)") +
        ylab(perturbationID) +
        theme_bw()
    return(plot)
}

#' Plot perturbations according to a given column representation
#'
#' @importFrom ggplot2 ggplot aes_string geom_point geom_hline xlab ylab
#'  theme element_blank scale_colour_manual xlim theme_classic guides
#' @importFrom ggrepel geom_text_repel
#'
#' @keywords internal
plotRes <- function(x, cellLine=NULL, n=5, signifThreshold=0.05, alpha=0.3) {
    if (is.null(cellLine)) {
        cellLines <- unique(attr(attr(x, "perturbations"), "metadata")$cell_id)
        summary <- grepl("Average", colnames(x))
        if (any(summary))
            cellLine <- colnames(x)[summary]
        else
            cellLine <- cellLines[[1]]
    }

    if (attr(x, "method") == "gsea")
        queryStr <- "WTCS"
    else
        queryStr <- c("coef", "qvalue")

    query  <- sprintf("^%s.*\\_%s$", cellLine, queryStr)
    stat   <- grep(query[[1]], colnames(x), value=TRUE, ignore.case=TRUE)
    statLabel <- sprintf("`%s`", stat)

    if (length(query) == 2) {
        signif <- grep(query[[2]], colnames(x), value=TRUE, ignore.case=TRUE)
        legend <- sprintf("p %s %s", c("<", ">="), signifThreshold)
        signifLabel <- sprintf("-log10(`%s`)", signif)

        x$thresh <- ifelse(x[[signif]] < signifThreshold,
                           legend[[1]], legend[[2]])
        colour <- TRUE
    } else {
        colour <- FALSE
    }

    plotType <- 2
    if (plotType == 1) {
        # Correlation coefficient vs q-value
        ggplot(x, aes_string(x=statLabel, y=signifLabel, colour="thresh")) +
            geom_point(alpha=alpha) +
            geom_rug(alpha=alpha) +
            scale_colour_manual(values=c("#56B4E9", "#999999")) +
            theme_classic() +
            theme(legend.title=element_blank(), legend.position="bottom")
    } else if (plotType == 2) {
        # Label perturbations
        sortedPert <- order(x[[stat]], decreasing=TRUE)
        index <- unique(c(head(sortedPert, n), tail(sortedPert, n)))
        x$label  <- ""
        x$label[index] <- x[[1]][index]

        # Correlation coefficient with labels for top and bottom perturbations
        if (colour) {
            vars <- aes_string(x=1, y=statLabel, label="label", colour="thresh")
            rug  <- geom_rug(alpha=alpha, sides="l")
        } else {
            vars <- aes_string(x=1, y=statLabel, label="label")
            rug  <- geom_rug(alpha=alpha, sides="l", colour="#56B4E9")
        }

        plot <- ggplot(x, vars) +
            rug +
            scale_colour_manual(values=c("#56B4E9", "#999999")) +
            geom_text_repel(nudge_x=0.15, direction="y", hjust=0,
                            segment.size=0.2, show.legend=FALSE, size=6) +
            xlim(1, 1.375) +
            theme_classic(base_size=18) +
            theme(legend.title=element_blank(), legend.position="bottom",
                  axis.line.x=element_blank(), axis.ticks.x=element_blank(),
                  axis.text.x=element_blank(), axis.title.x=element_blank())
        return(plot)
    }
}

#' Plot L1000 data comparison
#'
#' @param x \code{l1000comparison} object
#' @param perturbationID Character: perturbation identifier
#' @param genes Character: when plotting gene set enrichment analysis (GSEA),
#'   plot top genes (\code{genes = "top"}), bottom genes
#'   (\code{genes = "bottom"}) or both (\code{genes = "both"})
#'
#' @importFrom graphics plot
#'
#' @return Plot illustrating the comparison with L1000 data
#' @export
#'
#' @examples
#' data("l1000perturbationsKnockdown")
#' cellLine <- "HepG2"
#' compareKnockdown <- list()
#'
#' # Compare against L1000 using Spearman correlation
#' compareKnockdown$spearman <- compareAgainstL1000(
#'     diffExprStat, l1000perturbationsKnockdown, cellLine, method="spearman")
#'
#' # Compare against L1000 using Pearson correlation
#' compareKnockdown$pearson <- compareAgainstL1000(
#'     diffExprStat, l1000perturbationsKnockdown, cellLine, method="pearson")
#'
#' # Compare against L1000 using gene set enrichment analysis (GSEA) with the
#' # top and bottom 150 genes as gene sets
#' compareKnockdown$gsea <- compareAgainstL1000(
#'     diffExprStat, l1000perturbationsKnockdown, cellLine, method="gsea",
#'     geneSize=150)
#'
#' EIF4G1knockdown <- grep("EIF4G1", compareKnockdown$gsea$genes, value=TRUE)
#' plotL1000comparison(compareKnockdown$spearman, EIF4G1knockdown)
#' plotL1000comparison(compareKnockdown$pearson, EIF4G1knockdown)
#' plotL1000comparison(compareKnockdown$gsea, EIF4G1knockdown)
plotL1000comparison <- function(x, perturbation=NULL, perturbationID=NULL,
                                genes=c("both", "top", "bottom")) {
    method <- attr(x, "method")
    perturbation <- unclass(perturbation)

    hasPertID <- !is.null(perturbationID)
    isGSEA    <- method == "gsea"

    if (hasPertID) {
        if (length(perturbationID) == 0)
            stop("One perturbation ID must be provided")
        else if (length(perturbationID) > 1)
            stop("Only one perturbation ID is currently supported")
        else if (!perturbationID %in% colnames(perturbation))
            stop("Perturbation not found in the columns of the given dataset")

        perturbation <- perturbation[ , perturbationID]
        if (!isGSEA) {
            diffExprGenes <- attr(x, "diffExprGenes")
            plotSingleCorr(diffExprGenes, perturbation, perturbationID)
        } else {
            pathways <- attr(x, "pathways")
            genes <- match.arg(genes)
            plotGSEA(pathways, perturbation, genes=genes)
        }
    } else {
        plotRes(x, cellLine=NULL, n=3, signifThreshold=0.05)
    }
}

#' @rdname plotL1000comparison
#' @export
plot.l1000comparison <- plotL1000comparison
