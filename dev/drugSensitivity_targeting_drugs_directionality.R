# Load genomics and drug sensitivity associations for GDSC, CTRP and NCI60 -----

library(cTRAP)
gdsc  <- loadGenomicsDrugSensitivityAssociation("GDSC")
ctrp  <- loadGenomicsDrugSensitivityAssociation("CTRP")
nci60 <- loadGenomicsDrugSensitivityAssociation("NCI60")

compoundInfo <- list()
compoundInfo$gdsc  <- attr(gdsc, "compoundInfo")
compoundInfo$ctrp  <- attr(ctrp, "compoundInfo")
compoundInfo$nci60 <- attr(nci60, "compoundInfo")

# Comparing gene profile across common compounds between datasets --------------
stripStr <- function(str) gsub("[^[:alnum:] ]", "", as.character(str))

library(ggpubr)
library(pbapply)
plotCommonCompounds <- function(data1, data2, compoundInfo1, compoundInfo2) {
    genes <- intersect(rownames(data1), rownames(data2))

    compoundName <- list(data1=stripStr(name1), data2=stripStr(name2))
    compounds    <- intersect(compoundName$data1, compoundName$data2)
    message("Common compounds: ", length(compounds))

    id1   <- compoundInfo1$id
    id2   <- compoundInfo2$id
    name1 <- compoundInfo1$name
    name2 <- compoundInfo2$name

    plotCommonGenes <- function(compound, data1, data2, id1, id2, genes,
                                compoundName) {
        comp1 <- id1[match(compound, compoundName$data1)]
        comp2 <- id2[match(compound, compoundName$data2)]

        df <- data.frame(x=data1[genes, comp1], y=data2[genes, comp2])
        plot <- ggscatter(
            df, "x", "y", title=compound, xlab=comp1, ylab=comp2,
            add="reg.line", add.params=list(color="orange", fill="lightgray"),
            conf.int=TRUE, alpha=0.2) +
            stat_cor(method="pearson", color="orange")
    }
    pblapply(compounds, plotCommonGenes, data1=data1, data2=data2,
             id1=id1, id2=id2, genes=genes, compoundName=compoundName)
}

## GDSC vs NCI60
gdscVSnci60 <- plotCommonCompounds(data1=gdsc, data2=nci60,
                                   compoundInfo1=compoundInfo$gdsc,
                                   compoundInfo2=compoundInfo$nci60)
cowplot::plot_grid(plotlist=gdscVSnci60)

ctrpVSnci60 <- plotCommonCompounds(data1=ctrp, data2=nci60,
                                   compoundInfo1=compoundInfo$ctrp,
                                   compoundInfo2=compoundInfo$nci60)
cowplot::plot_grid(plotlist=ctrpVSnci60)

ctrpVSgdsc <- plotCommonCompounds(data1=ctrp, data2=gdsc,
                                  compoundInfo1=compoundInfo$ctrp,
                                  compoundInfo2=compoundInfo$gdsc)
cowplot::plot_grid(plotlist=ctrpVSgdsc)
