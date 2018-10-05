#' Plot gene set enrichment analysis (GSEA)
#'
#' @inheritParams fgsea::fgsea
#' @param title Character: plot title
#'
#' @importFrom fgsea calcGseaStat
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_line
#' scale_y_continuous scale_x_continuous theme element_rect element_blank
#' element_text labs ggtitle theme_bw geom_segment geom_bar
#' scale_colour_gradient2 guides element_line unit
#' @importFrom cowplot plot_grid
#' @importFrom stats na.omit
#'
#' @return Grid of plots illustrating a GSEA plot
plotGSEA <- function(pathways, stats, title="GSEA plot") {
    # Custom
    axis_title_size <- 12
    axis_text_size  <- 10
    plot_title_size <- 14
    point_size      <- 0.1
    gseaParam       <- 1

    rnk <- rank(-stats)
    ord <- order(rnk)

    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
    statsAdj <- statsAdj / max(abs(statsAdj))

    pathway <- unname(as.vector(na.omit(match(pathways, names(statsAdj)))))
    pathway <- sort(pathway)

    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                            returnAllExtremes = TRUE)

    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops

    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))

    diff <- (max(tops) - min(bottoms)) / 8

    # Getting rid of NOTEs
    x=y=NULL

    # Enrichment plot
    g1 <- ggplot(toPlot, aes(x=x, y=y)) +
        geom_point(color="green", size=point_size) +
        #geom_hline(yintercept=max(tops), colour="red", linetype="dashed") +
        geom_hline(yintercept=0, colour="black", size=0.3) +
        geom_line(color="green", size=1) + theme_bw() +
        scale_y_continuous(limits = c(min(bottoms), max(tops)), expand=c(0,0)) +
        scale_x_continuous(expand=c(0,0)) +
        theme(panel.border=element_rect(colour="black"),
              axis.title.x=element_blank(), axis.text.x=element_blank(),
              axis.ticks.x =element_blank(), axis.title.y = element_text(
                  size=axis_title_size, colour="black"),
              axis.text.y = element_text(size=axis_text_size, colour="black"),
              plot.title  = element_text(size=plot_title_size, colour="black",
                                         hjust=0.5)) +
        labs(y="Enrichment score (ES)") + ggtitle(title)

    # where are the genes
    g2 <- ggplot(toPlot, aes(x=x, y=y)) +
        geom_line(color=NA) + theme_bw() +
        geom_segment(data=data.frame(x=pathway),
                     mapping=aes(x=x, y=min(bottoms),
                                 xend=x, yend=max(tops)),
                     size=0.4) +
        scale_y_continuous(expand=c(0,0)) +
        scale_x_continuous(expand=c(0,0), breaks = seq(0,50000,2500)) +
        theme(panel.border=element_rect(colour="black"),
              panel.grid.minor=element_blank(),
              panel.grid.major =element_blank(),
              axis.title.y=element_text(colour="white"),
              axis.text.y=element_text(colour="white"),
              axis.ticks.y = element_line(colour="white"),
              axis.title.x=element_blank(), axis.text.x=element_blank(),
              axis.ticks.x =element_blank()) +
        labs(x="Rank", y="Enrichment score (ES)")

    # distribution of metric
    toPlot2 <- data.frame(x=seq(1:length(stats[ord])), y=stats[ord])

    g3 <- ggplot(toPlot2, aes(x,y)) +
        theme_bw() +
        geom_bar(stat="identity", aes(color=y), size=0.1) +
        scale_y_continuous(expand=c(0,0)) +
        scale_x_continuous(expand=c(0,0), breaks = seq(0,50000,2500)) +
        theme(panel.border=element_rect(colour="black"),
              axis.title.y = element_text(size=axis_title_size, colour="black"),
              axis.text.y = element_text(size=axis_text_size, colour="black"),
              axis.title.x = element_text(size=axis_title_size, colour="black"),
              axis.text.x = element_text(size=axis_text_size, colour="black")) +
        labs(y="Ranked metric", x="Rank") +
        scale_colour_gradient2(low = "dodgerblue", mid = "white",
                               high = "orangered", midpoint = 0) +
        guides(color=FALSE)

    plot_grid(g1 + theme(plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")),
              g2 + theme(plot.margin = unit(c(0, 5.5, 0, 5.5), "pt")),
              g3 + theme(plot.margin = unit(c(0, 5.5, 5.5, 5.5), "pt")),
              align = "v", ncol = 1, rel_heights = c(3/5, 0.5/5, 1.5/5))
}

#' Plot L1000 data comparison
#'
#' @param object \code{L1000 comparison} object
#' @param perturbationID Character: perturbation identifier
#' @param topGenes Boolean: plot top (\code{topGenes = TRUE}) or bottom genes
#' (\code{topGenes = FALSE}) in case of plotting gene set enrichment analysis
#' (GSEA)
#'
#' @importFrom graphics plot
#'
#' @return Plot illustrating the comparison with L1000 data
#' @export
#'
#' @examples
#' compareKnockdown <- cTRAP:::compareKnockdown
#' EIF4G1knockdown <- grep("EIF4G1", compareKnockdown$gsea_ordered$genes,
#'                         value=TRUE)
#' plotL1000comparison(compareKnockdown$spearman, EIF4G1knockdown)
#' plotL1000comparison(compareKnockdown$pearson, EIF4G1knockdown)
#' plotL1000comparison(compareKnockdown$gsea, EIF4G1knockdown)
plotL1000comparison <- function(object, perturbationID, topGenes=TRUE) {
    method <- attr(object, "method")
    perturbation <- attr(object, "perturbations")[ , perturbationID]
    if (method %in% c("spearman", "pearson")) {
        diffExprGenes <- attr(object, "diffExprGenes")

        # Intersect common genes
        genes <- intersect(names(perturbation), names(diffExprGenes))
        perturbation  <- perturbation[genes]
        diffExprGenes <- diffExprGenes[genes]

        plot(diffExprGenes[genes], perturbation[genes],
             xlab="Differentially expressed genes (input)",
             ylab=perturbationID)
    } else if (method == "gsea") {
        pathways <- attr(object, "pathways")
        if (topGenes) {
            end <- "top"
            title <- "%s (top %s genes)"
        } else {
            end <- "bottom"
            title <- "%s (bottom %s genes)"
        }

        pathway <- pathways[[grep(end, names(pathways))]]
        suppressWarnings(plotGSEA(pathway, perturbation,
                 title=sprintf(title, perturbationID, length(pathways[[1]]))))
    }
}
