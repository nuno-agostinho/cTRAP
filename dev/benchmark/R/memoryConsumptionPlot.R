#!/usr/bin/env Rscript

# Plot memory consumption
# Nuno Agostinho, 27 November 2020

# 1. Run heaptrack
# 2. Convert heaptrack output to massif
# 3. Get time stamp annotation from cTRAP logs
# 4. Create annotated plot of consumed memory

source("R/memoryConsumptionPlot_help.R")

# Load heaptrack data ----------------------------------------------------------
library(data.table)
massif <- list.files("output_heaptrack/cTRAP_1.8.1_3e3720d_1thread",
                     pattern="massif", full.names=TRUE)

heaptrack <- lapply(massif, readLines)
names(heaptrack) <- sapply(heaptrack, getValuesFromItem, "cmd: .*rds ",
                           fixed=FALSE)

# Assume time is displayed in seconds in all files
isTimeInSeconds <- sapply(heaptrack, getValuesFromItem, "time_unit: ") == "s"
stopifnot(all(isTimeInSeconds))
snapshots <- lapply(heaptrack, getMemorySnapshots)

# Peak heap memory consumption
sort(sapply(lapply(snapshots, "[[", "memory"), max))

# Load annotation data ---------------------------------------------------------
annot <- prepareAnnotation("logs/cTRAP_1.8.1_3e3720d/")
snapshots <- appendSnapshotAnnotation(snapshots, annot)
snapshots <- lapply(snapshots, addNAsBetweenCategories)
data <- Reduce(rbind, snapshots)
type <- rep(names(snapshots), sapply(snapshots, nrow))

data$chunking <- grepl("TRUE", type)
data$type     <- gsub("(.*) .* .*", "\\1", type)
data <- data[data$category != "finished!", ]
levels(data$category)[
  levels(data$category) == "starting to load CMap z-scores"] <- "data loading"
levels(data$category)[
  levels(data$category) == "comparing with pearson"] <- "pearson correlation"
levels(data$category)[
  levels(data$category) == "comparing with spearman"] <- "spearman correlation"
levels(data$category)[
  levels(data$category) == "comparing with gsea"] <- "gsea score"

table(data$category)

# Plot memory consumption ------------------------------------------------------
library(gridExtra)
plotMemConsumption(data) + facet_grid(chunking ~ type, scales="free")
