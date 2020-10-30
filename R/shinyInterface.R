# Skeleton for common elements -------------------------------------------------

#' Plot packed bubbles
#'
#' @importFrom highcharter hchart hcaes hc_title hc_tooltip %>% highchart JS
#'
#' @param data Data to plot
#' @param title Character: plot title
#' @param colour Character: bubble colour
#'
#' @return \code{highchart} object
#' @keywords internal
.plotBubbles <- function(data, title, colour="orange") {
    df <- as.data.frame(table(data))
    if (nrow(df) == 0) {
        hc <- highchart()
    } else {
        tooltip <- paste0("function() {",
                          "return '<b>' + this.series.name + '</b><br/>'",
                          "+ this.y + ' perturbations'; }")
        Freq <- NULL
        hc <- hchart(df, "packedbubble", hcaes(group=data, value=Freq),
                     showInLegend=FALSE, maxSize=50, color=colour) %>%
            hc_tooltip(formatter=JS(tooltip))
    }
    hc <- hc %>% hc_title(text=title)
    return(hc)
}

#' Prepare Shiny page template
#'
#' @importFrom shiny navbarPage tabPanel sidebarLayout tags
#'
#' @return HTML elements
#' @keywords internal
.prepareNavPage <- function(...) {
    app <- "cTRAP"
    ui  <- navbarPage(app, ...)

    ui[[3]][[1]][[2]]$class <- "navbar navbar-default"
    ui <- tags$div(class="container-fluid", style="padding-top: 15px;", ui)
    ui[[3]][[1]][[3]][[2]] <- ui[[3]][[1]][[3]][[2]][[3]][[1]]
    return(ui)
}


#' @importFrom DT datatable formatSignif
.prepareDT <- function(table, suffixes=NULL, ..., columnDefs=NULL) {
    if (is.null(table)) return(NULL)
    if (nrow(table) == 0) return(table)

    # Prepare columns whose significant digits will be rounded
    cols <- unlist(lapply(c("_coef", "_pvalue", "_qvalue", "GSEA", "WTCS"),
                          paste0, suffixes))
    cols <- lapply(cols, function(i) endsWith(names(table), i))
    cols <- which(Reduce("|", cols))
    cols <- cols[!sapply(table, class)[cols] %in% c("character", "logical")]

    # Convert to factor if there are low number of unique values
    lenUniqValuesPerCol <- sapply(sapply(table, unique), length)
    for (col in names(lenUniqValuesPerCol[lenUniqValuesPerCol < 50])) {
        if (!col %in% names(table)[cols]) { # Ignore columns that are to round
            table[[col]] <- as.factor(table[[col]])
        }
    }
    dt <- datatable(
        table, style="bootstrap", rownames=FALSE, ...,
        selection="single", extensions="Buttons", filter="top",
        options=list(dom="Bfrtip", buttons=list(list(extend="colvis")),
                     columnDefs=columnDefs))

    # Round significant digits
    if (length(cols) > 0) dt <- formatSignif(dt, cols)
    return(dt)
}

.prepareMetadata <- function(x) {
    input <- attr(x, "input")
    if (!is.null(input)) {
        attr(x, "input") <- as.data.frame(
            cbind("Element"=names(input), "Value"=input))
    }
    geneset <- attr(x, "geneset")
    if (!is.null(geneset)) {
        category           <- rep(names(geneset), sapply(geneset, length))
        attr(x, "geneset") <- as.data.frame(
            cbind("Element"=unname(unlist(geneset)), "Category"=category))
    }

    tables <- sapply(attributes(x), is, "data.frame")
    tables <- attributes(x)[tables]

    if (is(x, "perturbationChanges")) {
        data <- rbind(
            cbind("Z-scores filename", x),
            cbind("Source", attr(x, "source")),
            cbind("Type", attr(x, "type")),
            cbind("Gene size", length(attr(x, "genes"))),
            cbind("Perturbation number",
                  length(attr(x, "perturbations"))))
        colnames(data) <- c("Key", "Value")
        dataName <- "Object summary"
    } else {
        data <- tryCatch(as.table(x, clean=FALSE), error=return)
        if (is(data, "error")) data <- data.frame("Value"=x)
        dataName <- "Object values"
    }

    tables <- c(list(data), tables)
    names(tables)[[1]] <- dataName

    # Improve names of common tables
    name <- c("metadata"="Metadata",
              "geneInfo"="Gene information",
              "compoundInfo"="Compound information",
              "input"="Input",
              "rankingInfo"="Ranking information",
              "geneset"="Geneset",
              "extra"="Extra")
    pretty <- name[names(tables)]
    names(tables) <- unname(ifelse(is.na(pretty), names(tables), pretty))
    return(tables)
}

.prepareEllipsis <- function(...) {
    elems <- list(...)
    names(elems) <- sapply(substitute(list(...))[-1], deparse)
    return(elems)
}

# Internal interface server and UI ---------------------------------------------

#' @importFrom shiny NS sidebarPanel selectizeInput mainPanel tabPanel
#' sidebarLayout
#' @importFrom DT DTOutput
.diffExprENCODEloaderUI <- function(id, metadata, cellLine=NULL, gene=NULL,
                              title="Differential Expression Loader") {
    ns <- NS(id)
    cellLines <- sort(unique(metadata$`Biosample term name`))
    genes     <- sort(unique(metadata$`Experiment target`))
    genes     <- genes[genes != ""]

    onInitializeCellLine <- NULL
    if (is.null(cellLine)) {
        onInitializeCellLine <- I('function() { this.setValue(""); }')
    }

    onInitializeGene <- NULL
    if (is.null(gene)) {
        onInitializeGene <- I('function() { this.setValue(""); }')
    }

    sidebar <- sidebarPanel(
        selectizeInput(ns("cellLine"), "Cell line", cellLines, cellLine,
                       options=list(placeholder='Please select a cell line',
                                    onInitialize=onInitializeCellLine)),
        selectizeInput(ns("gene"), "Gene", genes, gene,
                       options=list(placeholder='Please select a gene',
                                    onInitialize=onInitializeGene)),
        actionButton(ns("load"), "Load data", class="btn-primary"))

    mainPanel <- mainPanel(DTOutput(ns("table")))
    ui <- tabPanel(title, sidebarLayout(sidebar, mainPanel))
    return(ui)
}

#' @importFrom shiny moduleServer stopApp
#' @importFrom DT renderDT
.diffExprENCODEloaderServer <- function(id, metadata) {
    moduleServer(
        id, function(input, output, session) {
            output$table <- renderDT({
                cellLine <- input$cellLine
                if (cellLine == "") cellLine <- NULL

                gene <- input$gene
                if (gene == "") gene <- NULL

                .prepareDT(downloadENCODEknockdownMetadata(cellLine, gene))
            })
            observeEvent(input$load, {
                message("Filtering data...")
                ENCODEmetadata <- downloadENCODEknockdownMetadata(
                    input$cellLine, input$gene)

                if (nrow(ENCODEmetadata) == 0) {
                    stopApp()
                    stop("No samples match the selected criteria")
                }
                ENCODEsamples <- loadENCODEsamples(ENCODEmetadata)[[1]]
                counts <- prepareENCODEgeneExpression(ENCODEsamples)

                # Remove low coverage (at least 10 counts shared by 2 samples)
                minReads   <- 10
                minSamples <- 2
                filter <- rowSums(counts[ , -c(1, 2)] >= minReads) >= minSamples
                counts <- counts[filter, ]

                # Convert ENSEMBL identifier to gene symbol
                message("Converting ENSEMBL identifiers to gene symbols...")
                counts$gene_id <- convertENSEMBLtoGeneSymbols(counts$gene_id)

                # Perform differential gene expression analysis
                message("Performing differential gene expression analysis...")
                diffExpr <- performDifferentialExpression(counts)

                diffExprStat <- diffExpr$t
                names(diffExprStat) <- diffExpr$Gene_symbol
                stopApp(diffExprStat)
            })
        }
    )
}

.findMatch <- function(what, where, ignore.case=TRUE) {
    if (is.null(what)) return(NULL)
    unname(sapply(what, grep, where, ignore.case=ignore.case, value=TRUE))
}

#' @importFrom shiny NS selectizeInput checkboxGroupInput sidebarPanel
.cmapDataLoaderUI <- function(id, metadata, zscores, geneInfo, compoundInfo,
                              cellLine, timepoint, dosage, perturbationType,
                              title="CMap Data Loader") {
    ns <- NS(id)
    dataTypes <- c("Perturbation metadata"="metadata",
                   "Perturbation z-scores"="zscores",
                   "Gene information"="geneInfo",
                   "Compound information"="compoundInfo")
    selectizeCondition <- function(id, label, choices, selected, ...) {
        selected <- .findMatch(selected, choices)
        plugins  <- list("remove_button", "restore_on_backspace")
        return(selectizeInput(id, label, choices, selected, ...,
                              options=list(plugins=plugins)))
    }

    conditions <- getCMapConditions(metadata)
    sidebar <- sidebarPanel(
        selectizeCondition(ns("type"), "Perturbation type",
                           conditions$perturbationType, perturbationType),
        selectizeCondition(ns("cellLine"), "Cell lines", multiple=TRUE,
                           conditions$cellLine, cellLine),
        selectizeCondition(ns("dosage"), "Dosages", multiple=TRUE,
                           conditions$dosage, dosage),
        selectizeCondition(ns("timepoint"), "Time points", multiple=TRUE,
                           conditions$timepoint, timepoint),
        checkboxGroupInput(ns("data"), "Data to load", dataTypes, dataTypes),
        actionButton(ns("cancel"), "Cancel"),
        actionButton(ns("load"), "Load data", class="btn-primary"))
    mainPanel <- mainPanel(
        fluidRow(column(4, highchartOutput(ns("cellLinePlot"),
                                           height="200px")),
                 column(4, highchartOutput(ns("dosagePlot"),
                                           height="200px")),
                 column(4, highchartOutput(ns("timepointPlot"),
                                           height="200px")),
        ),
        DTOutput(ns("table")))
    ui <- tabPanel(title, sidebar, mainPanel)
    return(ui)
}

#' @importFrom shiny isolate updateSelectizeInput moduleServer stopApp
#' observeEvent
.cmapDataLoaderServer <- function(id, metadata, zscores, geneInfo, compoundInfo,
                                  cellLine, timepoint, dosage) {
    updateSelectizeCondition <- function(session, id, choices, selected, ...) {
        selected <- .findMatch(selected, choices)
        return(updateSelectizeInput(session, id, choices=choices,
                                    selected=selected, ...))
    }

    moduleServer(
        id,
        function(input, output, session) {
            # Update conditions based on selected perturbation type
            observeEvent(input$type, {
                available <- getCMapConditions(metadata)
                updateSelectizeCondition(
                    session, "cellLine", selected=cellLine,
                    choices=available$cellLine)
                updateSelectizeCondition(
                    session, "dosage", selected=dosage,
                    choices=available$dosage)
                updateSelectizeCondition(
                    session, "timepoint", selected=timepoint,
                    choices=available$timepoint)
            })

            # Filter metadata based on selected inputs
            getFilteredMetadata <- reactive({
                filterCMapMetadata(
                    metadata, perturbationType=input$type,
                    cellLine=input$cellLine,
                    dosage=input$dosage,
                    timepoint=input$timepoint)
            })

            # Show plots
            output$cellLinePlot <- renderHighchart({
                subset <- getFilteredMetadata()
                .plotBubbles(subset$cell_id, "Cell lines", "orange")
            })
            output$dosagePlot <- renderHighchart({
                subset <- getFilteredMetadata()
                .plotBubbles(subset$pert_idose, "Dosages", "green")
            })
            output$timepointPlot <- renderHighchart({
                subset <- getFilteredMetadata()
                .plotBubbles(subset$pert_itime, "Time points", "purple")
            })

            # Show table
            output$table <- renderDT({
                hiddenCols <- c("pert_dose", "pert_dose_unit", "pert_time",
                                "pert_time_unit")
                hiddenCols <- match(hiddenCols, colnames(metadata))
                columnDefs <- list(list(visible=FALSE, targets=hiddenCols - 1))
                .prepareDT(getFilteredMetadata(), columnDefs=columnDefs)
            })

            # Load data
            observeEvent(input$load, {
                types <- isolate(input$data)

                returnIf <- function(bool, val) {
                    res <- NULL
                    if (bool) res <- val
                    return(res)
                }
                metadata     <- returnIf("metadata" %in% types,
                                         getFilteredMetadata())
                zscores      <- returnIf("zscores" %in% types, zscores)
                geneInfo     <- returnIf("geneInfo" %in% types, geneInfo)
                compoundInfo <- returnIf("compoundInfo" %in% types,
                                         compoundInfo)

                perturbations <- prepareCMapPerturbations(
                    metadata=metadata, zscores=zscores,
                    geneInfo=geneInfo, compoundInfo=compoundInfo)
                stopApp(perturbations)
            })

            # Cancel
            observeEvent(input$cancel, stopApp(stop("User cancel",
                                                    call.=FALSE)))
        }
    )
}

#' @importFrom shiny NS sidebarPanel mainPanel tabPanel sidebarLayout
#' selectizeInput
.metadataViewerUI <- function(id, x, title="Metadata Viewer") {
    ns <- NS(id)
    sidebar <- sidebarPanel(
        selectizeInput(ns("object"), "Dataset", names(x)),
        selectizeInput(ns("attr"), "Table", choices=NULL))
    main    <- mainPanel(DTOutput(ns("table")))
    ui      <- tabPanel(title, sidebarLayout(sidebar, main))
    return(ui)
}

#' @importFrom shiny moduleServer observe
#' @importFrom DT renderDT
.metadataViewerServer <- function(id, x) {
    moduleServer(
        id,
        function(input, output, session) {
            getSelectedObject <- reactive(.prepareMetadata(x[[input$object]]))

            observe(updateSelectizeInput(session, "attr",
                                         choices=names(getSelectedObject())))

            output$table <- renderDT(
                .prepareDT(getSelectedObject()[[input$attr]]))
        }
    )
}

#' @importFrom shiny NS sidebarPanel plotOutput selectizeInput mainPanel
#' tabPanel
#' @importFrom DT DTOutput
.dataPlotterUI <- function(id, x, title="Data Plotter") {
    ns <- NS(id)
    sidebar <- sidebarPanel(
        selectizeInput(ns("object"), "Dataset", names(x)),
        selectizeInput(ns("element"), "Element", choices=NULL),
        selectizeInput(ns("method"), "Method", choices=NULL))
    sidebar[[3]][[2]] <- plotOutput(ns("plot"))

    mainPanel <- mainPanel(DTOutput(ns("table")))
    ui <- tabPanel(title, sidebarLayout(sidebar, mainPanel))
    return(ui)
}

#' @importFrom shiny renderPlot observeEvent observe
#' @importFrom DT renderDT
.dataPlotterServer <- function(id, x) {
    moduleServer(
        id,
        function(input, output, session) {
            getSelectedObject <- reactive(x[[input$object]])

            # Update element and methods choices depending on selected object
            observeEvent(input$object, {
                obj <- getSelectedObject()
                updateSelectizeInput(session, "element", choices=obj[[1]])
                # Update methods depending on selected object
                methods <- c("Spearman's correlation"="spearman",
                             "Pearson's correlation"="pearson",
                             "GSEA"="gsea")
                isPresent <- sapply(methods, function(i)
                    any(grepl(i, colnames(obj), ignore.case=TRUE)))
                methods <- methods[isPresent]
                updateSelectizeInput(session, "method", choices=methods)
            })

            # Update selected element
            observe({
                selected <- getSelectedObject()[[1]][input$table_rows_selected]
                updateSelectizeInput(session, "element", selected=selected)
            })

            output$plot <- renderPlot({
                obj <- getSelectedObject()

                method <- input$method
                if (is.null(method) || method == "") return(NULL)

                element <- input$element
                if (element == "") element <- NULL
                if (!is.null(element) && !element %in% obj[[1]]) return(NULL)

                plot(obj, element, method=method, n=6)
            })

            output$table <- renderDT(
                .prepareDT(as.table(getSelectedObject(), clean=FALSE)))
        }
    )
}

#' @importFrom shiny NS sidebarPanel plotOutput selectizeInput mainPanel
#' tabPanel hr
#' @importFrom DT DTOutput
.datasetComparisonUI <- function(id, x, title="Dataset Comparison") {
    ns <- NS(id)
    sidebar <- sidebarPanel(
        selectizeInput(ns("data1"), "Dataset 1", names(x)),
        selectizeInput(ns("col1"), "Column to plot in X axis", choices=NULL),
        hr(),
        selectizeInput(ns("data2"), "Dataset 2", names(x), selected=2),
        selectizeInput(ns("col2"), "Column to plot in Y axis", choices=NULL))
    sidebar[[3]][[2]] <- plotOutput(ns("plot"))

    mainPanel <- mainPanel(DTOutput(ns("table")))
    ui <- tabPanel(title, sidebarLayout(sidebar, mainPanel))
    return(ui)
}


#' @importFrom shiny renderPlot observe
#' @importFrom DT renderDT
.datasetComparisonServer <- function(id, x) {
    moduleServer(
        id,
        function(input, output, session) {
            getSelectedDataset1 <- reactive(x[[input$data1]])
            getSelectedDataset2 <- reactive(x[[input$data2]])

            observe({
                dataset1 <- getSelectedDataset1()
                updateSelectizeInput(session, "col1", choices=names(dataset1))
            })

            observe({
                dataset2 <- getSelectedDataset2()
                updateSelectizeInput(session, "col2", choices=names(dataset2))
            })

            output$plot <- renderPlot({
                isColValid <- function(col, dataset) {
                    !is.null(col) && col != "" && col %in% colnames(dataset)
                }

                dataset1 <- getSelectedDataset1()
                col1 <- input$col1
                if (!isColValid(col1, dataset1)) return(NULL)

                dataset2 <- getSelectedDataset2()
                col2 <- input$col2
                if (!isColValid(col2, dataset2)) return(NULL)

                if (any(dataset1[[1]] %in% dataset2[[1]])) {
                    return(plot(dataset1[[col1]], dataset2[[col2]]))
                } else {
                    stop("no common identifiers between datasets")
                }
            })
            # output$table <- renderDT(
            #     .prepareDT(as.table(getSelectedObject(), clean=FALSE)))
        }
    )
}

.targetingDrugsVSsimilarPerturbationsPlotterUI <- function(
    id, x, elemClasses, title="Dataset Comparison") {
    if (!all(c("targetingDrugs", "similarPerturbations") %in% elemClasses)) {
        return(NULL)
    }
    ns <- NS(id)
    sidebar <- sidebarPanel(
        selectizeInput(
            ns("data1"), "Dataset with predicted targeting drugs",
            names(x)[elemClasses == "targetingDrugs"]),
        # selectizeInput(ns("col1"), "Column to plot in X axis", choices=NULL),
        selectizeInput(
            ns("data2"), "Dataset with similar perturbations",
            names(x)[elemClasses == "similarPerturbations"]))
        # selectizeInput(ns("col2"), "Column to plot in Y axis", choices=NULL))

    mainPanel <- mainPanel(plotOutput(ns("plot"), brush=ns("brush")),
                           DTOutput(ns("table")))
    ui <- tabPanel(title, sidebarLayout(sidebar, mainPanel))
    return(ui)
}

#' @importFrom shiny moduleServer reactive observe renderPlot brushedPoints
#' @importFrom DT renderDT
.targetingDrugsVSsimilarPerturbationsPlotterServer <- function(id, x) {
    moduleServer(
        id,
        function(input, output, session) {
            getSelectedDataset1 <- reactive(x[[input$data1]])
            getSelectedDataset2 <- reactive(x[[input$data2]])

            observe({
                plot <- suppressMessages(
                    plotTargetingDrugsVSsimilarPerturbations(
                        getSelectedDataset1(), getSelectedDataset2(),
                        column="spearman_rank", labelBy=NULL) + theme_bw(16))

                output$plot  <- renderPlot(plot)
                output$table <- renderDT({
                    data  <- attr(plot, "data")
                    brush <- input$brush
                    if (!is.null(brush)) data <- brushedPoints(data, brush)

                    hiddenCols <- grep("^pearson|GSEA|spearman", colnames(data))
                    columnDefs <- list(list(visible=FALSE,
                                            targets=hiddenCols - 1))
                    return(.prepareDT(data, attr(plot, "suffixes"),
                                      columnDefs=columnDefs))
                })
            })
        }
    )
}

# Launch graphical interface ---------------------------------------------------

#' Load differential expression data via a visual interface
#'
#' Currently only supports loading data from ENCODE knockdown experiments
#'
#' @inheritParams downloadENCODEknockdownMetadata
#'
#' @return Differential expression data
#' @family visual interface functions
#' @export
launchDiffExprLoader <- function(cellLine=NULL, gene=NULL) {
    metadata <- downloadENCODEknockdownMetadata()
    id       <- "diffExpr"
    ui       <- .prepareNavPage(
        .diffExprENCODEloaderUI(id, metadata, cellLine, gene))
    server   <- function(input, output, session) {
        .diffExprENCODEloaderServer(id, metadata)
    }
    app    <- runApp(shinyApp(ui, server))
    return(app)
}

#' Load CMap data via a visual interface
#'
#' @importFrom highcharter highchartOutput renderHighchart
#' @importFrom shiny sidebarPanel selectizeInput actionButton mainPanel fluidRow
#' column observeEvent updateSelectizeInput reactive runApp shinyApp tabPanel
#' @importFrom DT DTOutput renderDT
#'
#' @inheritParams getCMapConditions
#' @inheritParams prepareCMapPerturbations
#' @inherit prepareCMapPerturbations return
#'
#' @return CMap data
#' @family visual interface functions
#' @export
launchCMapDataLoader <- function(metadata="cmapMetadata.txt",
                                 zscores="cmapZscores.gctx",
                                 geneInfo="cmapGeneInfo.txt",
                                 compoundInfo="cmapCompoundInfo.txt",
                                 cellLine=NULL, timepoint=NULL, dosage=NULL,
                                 perturbationType=NULL) {
    metadata <- loadCMapData(metadata, type="metadata")
    id <- "cmapDataLoader"
    ui <- .prepareNavPage(
        .cmapDataLoaderUI(id, metadata, zscores, geneInfo, compoundInfo,
                          cellLine, timepoint, dosage, perturbationType))
    server <- function(input, output, session) {
        .cmapDataLoaderServer(id, metadata, zscores, geneInfo, compoundInfo,
                              cellLine, timepoint, dosage)
    }
    app <- runApp(shinyApp(ui, server))
    return(app)
}

#' View metadata via a visual interface
#'
#' @param ... Objects
#'
#' @importFrom shiny runApp shinyApp
#'
#' @return Metadata viewer (retunrs \code{NULL})
#' @family visual interface functions
#' @export
launchMetadataViewer <- function(...) {
    elems  <- .prepareEllipsis(...)
    id     <- "metadataViewer"
    ui     <- .prepareNavPage(.metadataViewerUI(id, elems))
    server <- function(input, output, session) .metadataViewerServer(id, elems)
    app    <- runApp(shinyApp(ui, server))
    return(app)
}

#' View and plot results via a visual interface
#'
#' @param ... Objects
#'
#' @importFrom shiny tagList
#'
#' @return Launches result viewer and plotter (returns \code{NULL})
#' @family visual interface functions
#' @export
launchResultPlotter <- function(...) {
    elems <- .prepareEllipsis(...)
    compareId     <- "datasetComparison"
    comparePlotId <- "comparePlotter"
    dataId        <- "dataPlotter"
    metadataId    <- "metadataViewer"

    elemClasses       <- sapply(lapply(list(...), class), "[[", 1)
    hasSimilarPerts   <- "similarPerturbations" %in% elemClasses
    hasTargetingDrugs <- "targetingDrugs" %in% elemClasses
    showTwoKindPlot   <- hasSimilarPerts && hasTargetingDrugs

    uiList <- tagList(
        # .datasetComparisonUI(compareId, elems),
        .targetingDrugsVSsimilarPerturbationsPlotterUI(
            comparePlotId, elems, elemClasses),
        .dataPlotterUI(dataId, elems),
        .metadataViewerUI(metadataId, elems))
    uiList <- Filter(length, uiList)
    ui     <- do.call(.prepareNavPage, uiList)

    server <- function(input, output, session) {
        # .datasetComparisonServer(compareId, elems)
        if (showTwoKindPlot) {
            .targetingDrugsVSsimilarPerturbationsPlotterServer(
                comparePlotId, elems)
        }
        .dataPlotterServer(dataId, elems)
        .metadataViewerServer(metadataId, elems)
    }
    app <- runApp(shinyApp(ui, server))
    return(app)
}
