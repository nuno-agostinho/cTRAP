descriptors <- list()
descriptors$cmap2D  <- loadDrugDescriptors("CMap", "2D")
descriptors$cmap3D  <- loadDrugDescriptors("CMap", "3D")
descriptors$nci602D <- loadDrugDescriptors("NCI60", "2D")
descriptors$nci603D <- loadDrugDescriptors("NCI60", "3D")

# How many are categorical or numeric variables?
checkColClasses <- function(i) table(sapply(i, class))
lapply(descriptors, checkColClasses)

# How many and what categories on the categorical variables?
lapply(descriptors, function(i) {
    isCategorical <- sapply(i, class) == "character"
    isCategorical[1] <- FALSE # Ignore first column (identifier)
    data       <- i[ , isCategorical, with=F]
    categories <- lapply(data, table)
    len        <- lapply(categories, length)
    ifelse(len < 5, lapply(data, table)[len < 5], len)
})

# Why does cmap2D, but not nci602D, has duplicated categories whose values are
# not identical?
par(mfrow=c(1, 1))
plot(table(abs(
    descriptors$cmap2D$`H-Acceptors` - descriptors$cmap2D$H.Acceptors)))
plot(table(abs(
    descriptors$cmap2D$`Rotatable Bonds` - descriptors$cmap2D$Rotatable.Bonds)))
plot(table(abs(
    descriptors$cmap2D$Aromatic.Rings - descriptors$cmap2D$`Aromatic Rings`)))

# How many values on the non-categorical variables?
lapply(descriptors[c("cmap2D", "nci602D")], function(i) {
    isNotCategorical <- sapply(i, class) != "character"
    isNotCategorical[1] <- FALSE # Ignore first column (identifier)
    data <- i[ , isNotCategorical, with=F]
    lapply(lapply(data, table), length)
})

# How many sets to create for numeric variables?
calculateBins <- function(numbers, target.bins=15, minpts=NULL, ...) {
    numbers <- na.omit(numbers)
    print(diptest::dip.test(numbers)) # Test if multimodal
    factors <- calculateEvenlyDistributedBins(numbers, target.bins=target.bins,
                                              minpts=minpts, ...)
    par(mfrow=c(1, 2))
    plot(table(factors))
    plot(density(numbers), main = NA)
    return(invisible(factors))
}
# CMap 2D
calculateBins(descriptors$cmap2D$H.Acceptors)
calculateBins(descriptors$cmap2D$H.Donors)
calculateBins(descriptors$cmap2D$`Basic Nitrogens`)
calculateBins(descriptors$cmap2D$`Acidic Oxygens`)
calculateBins(descriptors$cmap2D$Amines)
calculateBins(descriptors$cmap2D$Amides)
calculateBins(descriptors$cmap2D$Polar.Surface.Area)
calculateBins(descriptors$cmap2D$`Symmetric atoms`)
calculateBins(descriptors$cmap2D$XLogP)
calculateBins(descriptors$cmap2D$Rotatable.Bonds)
calculateBins(descriptors$cmap2D$Total.Molweight)

# NCI60 2D
calculateBins(descriptors$nci602D$`Total Molweight`)
calculateBins(descriptors$nci602D$`Monoisotopic Mass`)
calculateBins(descriptors$nci602D$`H-Acceptors`)
calculateBins(descriptors$nci602D$`Polar Surface Area`)
calculateBins(descriptors$nci602D$`Aromatic Nitrogens`)
calculateBins(descriptors$nci602D$`Shape Index`)
calculateBins(descriptors$nci602D$`Acidic Oxygens`)
calculateBins(descriptors$nci602D$`Small Rings`)

# CMap 3D
calculateBins(descriptors$cmap3D$nAcid)
calculateBins(descriptors$cmap3D$nBase)
calculateBins(descriptors$cmap3D$ABC)
calculateBins(descriptors$cmap3D$ABCGG)
calculateBins(descriptors$cmap3D$SpAbs_A)
calculateBins(descriptors$cmap3D$mZagreb1)
calculateBins(descriptors$cmap3D$nBridgehead)
calculateBins(descriptors$cmap3D$nAtom)

# NCI60 3D
calculateBins(descriptors$nci603D$nAcid)
calculateBins(descriptors$nci603D$nBase)
calculateBins(descriptors$nci603D$ABC)
calculateBins(descriptors$nci603D$ABCGG)
calculateBins(descriptors$nci603D$SpAbs_A)
calculateBins(descriptors$nci603D$mZagreb1)
calculateBins(descriptors$nci603D$nBridgehead, minpts = 500)
calculateBins(descriptors$nci603D$nAtom)

# Nasty functions: 73 categories... include or not?
table(descriptors$cmap2D$`Nasty Functions`)

# Prepare drug sets from molecular descriptors
drugSets <- list()
drugSets$cmap2D  <- prepareDrugSets(descriptors$cmap2D)
drugSets$nci602D <- prepareDrugSets(descriptors$nci602D)
drugSets$cmap3D  <- prepareDrugSets(descriptors$cmap3D)
drugSets$nci603D <- prepareDrugSets(descriptors$nci603D)
