#!/usr/bin/env Rscript

# Plot memory consumption (auxiliary functions)
# Nuno Agostinho, 27 November 2020

getValuesFromItem <- function(raw, str, fixed=TRUE) {
  set <- grep(str, raw, fixed=fixed, value=TRUE)
  res <- gsub(str, "", set, fixed=fixed)
  return(res)
}

getMemorySnapshots <- function(raw) {
  time <- as.numeric(getValuesFromItem(raw, "time="))
  mem  <- as.numeric(getValuesFromItem(raw, "mem_heap_B="))
  df   <- data.frame("time"=time,    # in minutes
                     "memory"=mem / 10e8) # in GB
  return(df)
}

prepareAnnotation <- function(folder) {
  annot <- list.files(folder, pattern="heaptrack.*1-thread", full.names=TRUE)
  names(annot) <- gsub("heaptrack_rankCMap_(.*)_load-(.*)_(.*)-threads.log",
                       "\\1 \\2 \\3", basename(annot))
  annot <- lapply(annot, readLines)
  
  dateRegex <- "\\d{1,4}-\\d{1,2}-\\d{1,2} \\d{1,2}:\\d{1,2}:\\d{1,2}"
  annot <- lapply(annot, grep, pattern=dateRegex, value=TRUE)
  
  convertAnnotationToTimepoints <- function(annot, dateRegex) {
    dateRegex   <- sprintf("(%s) - (.*)", dateRegex)
    info        <- gsub(dateRegex, "\\2", annot)
    time        <- as.POSIXct(gsub(dateRegex, "\\1", annot))
    names(info) <- time - time[[1]]
    
    info <- info[!names(info) == 0]
    info <- info[!info %in% c("starting...", "loaded CMap z-scores")]
    
    # Join result preparation
    res <- grep("preparing .* results", info)
    info[res] <- "preparing results"
    names(info)[res] <- min(as.numeric(names(info)[res]))
    
    info <- info[!duplicated(names(info), fromLast=TRUE)]
    df   <- data.frame(time=names(info), info=info)
    return(df)
  }
  
  annot <- lapply(annot, convertAnnotationToTimepoints, dateRegex)
  return(annot)
}

addCategory <- function(item, snapshots, annot) {
  snap   <- snapshots[[item]]
  time   <- snap$time
  breaks <- unique(c(0, as.numeric(annot[[item]]$time), max(time)))
  bins   <- cut(time, breaks)
  levels(bins) <- c("data loading", annot[[item]]$info, "end")[seq(breaks)]
  snap$category <- bins
  return(snap)
}

appendSnapshotAnnotation <- function(snapshots, annot) {
  stopifnot(all(names(snapshots) %in% names(annot)))
  ns <- setNames(names(snapshots), names(snapshots))
  snapshots <- lapply(ns, addCategory, snapshots, annot)
  return(snapshots)
}

addNAsBetweenCategories <- function(snap) {
  # Add NAs between each category's transition to avoid connecting lines/areas
  transitions   <- cumsum(rle(as.character(snap$category))[[1]])
  naRows        <- snap[transitions, ]
  naRows$time   <- naRows$time + 0.001
  naRows$memory <- NA
  snap          <- rbind(snap, naRows)
  return(snap)
}

library(ggplot2)
plotMemConsumption <- function(data) {
  ggplot(data, aes(.data[["time / 60"]], .data[["memory"]], ymax="memory",
                   colour="category", fill="category")) +
    xlab("Time (minutes)") +
    ylab("Memory (GB)") +
    geom_ribbon(alpha=0.8, ymin=0) +
    # geom_line() +
    # ggtitle(type) +
    theme_light()
}
