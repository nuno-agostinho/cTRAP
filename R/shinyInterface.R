# Skeleton for common elements -------------------------------------------------

#' @importFrom shiny div h3 tags
.panel <- function(
    title=NULL, ..., footer=NULL, collapse=TRUE,
    type=c("default", "primary", "success", "info", "warning", "danger")) {

    type <- match.arg(type)
    if (!is.null(title)) {
        if (collapse) {
            title <- tags$a(role="button", "data-toggle"="collapse",
                            "data-parent"="#accordion", href="#collapseOne",
                            "aria-expanded"="false",
                            "aria-controls"="collapseOne", title)
        }
        title <- div(class="panel-heading", h3(class="panel-title", title))
    }

    body <- list(...)
    if (length(body) > 0) body <- div(class="panel-body", ...)
    if (collapse) {
        body <- div(id="collapseOne", class="panel-collapse collapse in",
                    role="tabpanel", "aria-labelledby"="headingOne", body)
    }
    if (!is.null(footer)) footer <- div(class="panel-footer", footer)
    div(class=paste0("panel panel-", type), title, body, footer)
}

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
    if (!is.table(data)) data <- table(data)
    df <- as.data.frame(data)
    names(df) <- c("data", "Freq")
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
.prepareDT <- function(table, suffixes=NULL, ..., columnDefs=NULL,
                       scrollX=TRUE, pagingType="simple_numbers") {
    if (is.null(table)) return(NULL)

    cols <- NULL
    if (nrow(table) > 0) {
        # Prepare columns whose significant digits will be rounded
        cols <- unlist(lapply(c("_coef", "_pvalue", "_qvalue", "GSEA", "WTCS"),
                              paste0, suffixes))
        cols <- c(cols, "pval", "padj", "ES", "NES") # DSEA results
        cols <- lapply(cols, function(i) endsWith(colnames(table), i))
        cols <- which(Reduce("|", cols))
        cols <- cols[!sapply(table, class)[cols] %in% c("character", "logical")]

        # Convert to factor if there are low number of unique values
        if (!all(c("Key", "Value") %in% colnames(table))) {
            lenUniqValuesPerCol <- sapply(sapply(table, unique), length)
            for (col in names(lenUniqValuesPerCol[lenUniqValuesPerCol < 50])) {
                if (!col %in% names(table)[cols]) { # Ignore columns to round
                    table[[col]] <- as.factor(table[[col]])
                }
            }
        }
    }

    # # Fix specific bug with targetingDrugs columns when using NCI-60 data
    # if (all(c("PubChem SID", "PubChem CID") %in% colnames(table))) {
    #     table$`PubChem CID` <- as.numeric(table$`PubChem CID`)
    #     table$`PubChem SID` <- as.numeric(table$`PubChem SID`)
    # }

    dt <- datatable(
        table, rownames=FALSE, ..., escape=FALSE,
        selection="single", extensions="Buttons", filter="top",
        options=list(dom="Bfrtip", buttons=list(list(extend="colvis")),
                     columnDefs=columnDefs, scrollX=scrollX,
                     pagingType=pagingType))
    # Round significant digits
    if (length(cols) > 0) dt <- formatSignif(dt, cols)
    return(dt)
}

.prepareMetadata <- function(x) {
    input <- attr(x, "input")
    if (!is.null(input)) {
        attr(x, "input") <- data.frame("Element"=names(input), "Value"=input)
    }
    geneset <- attr(x, "geneset")
    if (!is.null(geneset)) {
        category           <- rep(names(geneset), sapply(geneset, length))
        attr(x, "geneset") <- data.frame("Element"=unname(unlist(geneset)),
                                         "Category"=category)
    }

    tables <- sapply(attributes(x), is, "data.frame")
    tables <- attributes(x)[tables]

    if (is(x, "perturbationChanges")) {
        data <- rbind(
            cbind("Z-scores filename", prepareWordBreak(x)),
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
        dataName <- "Data"
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

.getENCODEconditions <- function(metadata) {
    cellLines <- sort(unique(metadata$`Biosample term name`))
    genes     <- sort(unique(metadata$`Experiment target`))
    genes     <- genes[genes != ""]

    res <- list(genes=genes, cellLines=cellLines)
    return(res)
}

#' @importFrom shiny NS sidebarPanel selectizeInput mainPanel tabPanel
#' sidebarLayout
#' @importFrom DT DTOutput
.diffExprENCODEloaderUI <- function(id, metadata, cellLine=NULL, gene=NULL,
                                    title="ENCODE Knockdown Data Loader") {
    ns <- NS(id)

    conditions <- .getENCODEconditions(metadata)
    cellLines  <- conditions$cellLines
    genes      <- conditions$genes

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
                .prepareDT(downloadENCODEknockdownMetadata())
            })
            proxy <- dataTableProxy("table")
            observe({
                cellLine <- input$cellLine
                if (cellLine == "") cellLine <- NULL

                gene <- input$gene
                if (gene == "") gene <- NULL
                data <- downloadENCODEknockdownMetadata(cellLine, gene)
                observe(replaceData(proxy, data, rownames=FALSE))
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

#' @importFrom shiny NS selectizeInput checkboxGroupInput sidebarPanel helpText
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
        plugins  <- list("remove_button")
        return(selectizeInput(id, label, choices, selected, ...,
                              options=list(plugins=plugins)))
    }

    plots <- fluidRow(
        column(3, highchartOutput(ns("perturbationPlot"), height="200px")),
        column(3, highchartOutput(ns("cellLinePlot"), height="200px")),
        column(3, highchartOutput(ns("dosagePlot"), height="200px")),
        column(3, highchartOutput(ns("timepointPlot"), height="200px")))
    conditions <- getCMapConditions(metadata, control=TRUE)
    sidebar <- sidebarPanel(
        selectizeCondition(ns("type"), "Perturbation type", multiple=TRUE,
                           conditions$perturbationType, perturbationType),
        selectizeCondition(ns("cellLine"), "Cell lines", multiple=TRUE,
                           conditions$cellLine, cellLine),
        selectizeCondition(ns("dosage"), "Dosages", multiple=TRUE,
                           conditions$dosage, dosage),
        selectizeCondition(ns("timepoint"), "Time points", multiple=TRUE,
                           conditions$timepoint, timepoint),
        checkboxGroupInput(ns("data"), "Data to load", dataTypes, dataTypes),
        helpText("By default, data will be downloaded if not found."),
        actionButton(ns("cancel"), "Cancel"),
        actionButton(ns("load"), "Load data", class="btn-primary"))
    mainPanel <- mainPanel(DTOutput(ns("table")))

    ui <- tabPanel(title, sidebar, mainPanel)
    ui[[3]] <- c(list(plots), ui[[3]])
    return(ui)
}

#' @importFrom shiny isolate updateSelectizeInput moduleServer stopApp
#' observeEvent
#' @importFrom DT dataTableProxy replaceData
.cmapDataLoaderServer <- function(id, metadata, zscores, geneInfo, compoundInfo,
                                  cellLine, timepoint, dosage) {
    updateSelectizeCondition <- function(session, input, id, choices, ...) {
        selected <- isolate(input[[id]])
        selected <- .findMatch(selected, choices)
        return(updateSelectizeInput(session, id, choices=choices,
                                    selected=selected, ...))
    }

    server <- function(input, output, session) {
        # Update conditions based on selected perturbation type
        observe({
            available <- getCMapConditions(metadata,
                                           perturbationType=input$type)
            updateSelectizeCondition(session, input, "cellLine",
                                     choices=available$cellLine)
            updateSelectizeCondition(session, input, "dosage",
                                     choices=available$dosage)
            updateSelectizeCondition(session, input, "timepoint",
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
        output$perturbationPlot <- renderHighchart({
            subset <- getFilteredMetadata()
            labels <- getCMapPerturbationTypes(control=TRUE)
            types  <- table(subset$pert_type)
            names(types) <- names(labels)[match(names(types), labels)]
            .plotBubbles(types, "Perturbations", "red")
        })

        renderBubbleChart <- function(subset, title, colour) {
            renderHighchart({
                data <- getFilteredMetadata()[[subset]]
                data[is.na(data)] <- "NA"
                .plotBubbles(data, title, colour)
            })
        }
        output$cellLinePlot  <- renderBubbleChart("cell_id", "Cell lines",
                                                  "orange")
        output$dosagePlot    <- renderBubbleChart("pert_idose", "Dosages",
                                                  "green")
        output$timepointPlot <- renderBubbleChart("pert_itime", "Time points",
                                                "purple")

        # Show table
        output$table <- renderDT({
            hiddenCols <- c("pert_dose", "pert_dose_unit", "pert_time",
                            "pert_time_unit")
            hiddenCols <- match(hiddenCols, colnames(metadata))
            columnDefs <- list(list(visible=FALSE, targets=hiddenCols - 1))
            .prepareDT(isolate(getFilteredMetadata()), columnDefs=columnDefs,
                               scrollX=TRUE)
        })
        proxy <- dataTableProxy("table")
        observe(replaceData(proxy, getFilteredMetadata(), rownames=FALSE))

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
            compoundInfo <- returnIf("compoundInfo" %in% types, compoundInfo)

            perturbations <- prepareCMapPerturbations(
                metadata=metadata, zscores=zscores,
                geneInfo=geneInfo, compoundInfo=compoundInfo)
            stopApp(perturbations)
        })

        # Cancel
        observeEvent(input$cancel, stopApp(stop("User cancel", call.=FALSE)))
    }

    moduleServer(id, server)
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
        selectizeInput(ns("method"), "Method", choices=NULL))
    sidebar[[3]][[2]] <- tagList(
        selectizeInput(ns("element"), "Row ID to plot", choices=NULL,
                       width="100%"),
        plotOutput(ns("plot"), brush=ns("brush")),
        DTOutput(ns("pointTable")))

    mainPanel <- mainPanel(DTOutput(ns("table")))
    ui <- tabPanel(title, sidebarLayout(sidebar, mainPanel))
    return(ui)
}

#' @importFrom shiny renderPlot observeEvent observe
#' @importFrom DT renderDT formatSignif dataTableProxy replaceData
.dataPlotterServer <- function(id, x) {
    isValid <- function(e) !is.null(e) && e != ""
    moduleServer(
        id,
        function(input, output, session) {
            getSelectedObject <- reactive(x[[input$object]])

            # Update element and methods choices depending on selected object
            observeEvent(input$object, {
                obj <- getSelectedObject()
                updateSelectizeInput(session, "element", choices=obj[[1]],
                                     selected=list(), server=TRUE)
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

            output$table <- renderDT({
                data <- as.table(getSelectedObject(), clean=FALSE)
                .prepareDT(data)
            })

            # Filter table based on overall plot
            proxy <- dataTableProxy("table")
            observe({
                elem <- input$element
                if (is.null(elem) || elem == "") {
                    obj   <- as.table(getSelectedObject(), clean=FALSE)
                    brush <- input$brush
                    if (!is.null(brush)) {
                        val  <- obj[[brush$mapping$y]]
                        rows <- val >= brush$ymin & val <= brush$ymax
                        obj  <- obj[rows, ]
                    }
                    replaceData(proxy, obj, rownames=FALSE)
                }
            })

            plotData <- reactive({
                obj <- getSelectedObject()
                method <- input$method
                if (!isValid(method)) return(NULL)
                element <- input$element
                if (element == "") element <- NULL
                if (!is.null(element) && !element %in% obj[[1]]) return(NULL)
                plot(obj, element, method=method, n=6)
            })

            output$plot <- renderPlot(plotData())
            output$pointTable <- renderDT({
                data <- attr(plotData(), "data")
                if (is.null(data)) return(NULL)
                brush <- input$brush
                if (!is.null(brush)) data <- brushedPoints(data, brush)
                numericCols <- sapply(data, is.numeric)
                dt <- .prepareDT(data, scrollX=TRUE, pagingType="simple",
                                 class="compact hover stripe")
                return(formatSignif(dt, numericCols))
            })
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
    sidebar[[3]][[2]] <- tagList(
        plotOutput(ns("plot"), brush=ns("brush")),
        helpText("Click-and-drag points in the plot to filter the table."))

    mainPanel <- mainPanel(DTOutput(ns("table")))
    ui <- tabPanel(title, sidebarLayout(sidebar, mainPanel))
    return(ui)
}

#' @importFrom shiny renderPlot observe
#' @importFrom DT renderDT
#' @importFrom ggplot2 theme_bw geom_density_2d geom_rug geom_point
.datasetComparisonServer <- function(id, x) {
    getNumericCols <- function(x) colnames(x)[vapply(x, is.numeric, logical(1))]

    moduleServer(
        id,
        function(input, output, session) {
            getSelectedDataset1 <- reactive(x[[input$data1]])
            getSelectedDataset2 <- reactive(x[[input$data2]])

            observe({
                dataset1 <- getSelectedDataset1()
                numericCols <- getNumericCols(dataset1)
                updateSelectizeInput(session, "col1", choices=numericCols)
            })

            observe({
                dataset2 <- getSelectedDataset2()
                numericCols <- getNumericCols(dataset2)
                updateSelectizeInput(session, "col2", choices=numericCols)
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
                    plot <- ggplot(data=NULL, aes(x=dataset1[[col1]],
                                                  y=dataset2[[col2]])) +
                        geom_rug(alpha=0.1) +
                        geom_abline(slope=1, intercept=0, colour="orange") +
                        geom_point(size=0.1, alpha=0.3) +
                        geom_density_2d(alpha=0.3) +
                        xlab(paste(input$data1, input$col1, sep="_")) +
                        ylab(paste(input$data2, input$col2, sep="_")) +
                        theme_bw()
                    return(plot)
                } else {
                    stop("no common identifiers between datasets")
                }
            })

            output$table <- renderDT({
                isColValid <- function(col, dataset) {
                    !is.null(col) && col != "" && col %in% colnames(dataset)
                }

                dataset1 <- getSelectedDataset1()
                col1 <- input$col1
                if (!isColValid(col1, dataset1)) return(NULL)

                dataset2 <- getSelectedDataset2()
                col2 <- input$col2
                if (!isColValid(col2, dataset2)) return(NULL)
                common <- dataset1[[1]] %in% dataset2[[1]]

                if (any(common)) {
                    df <- data.frame(dataset1[[1]][common],
                                     getSelectedDataset1()[[col1]][common],
                                     getSelectedDataset2()[[col2]][common])
                    colnames(df) <- c("id",
                                      paste(input$data1, input$col1, sep="_"),
                                      paste(input$data2, input$col2, sep="_"))

                    brush <- input$brush
                    if (!is.null(brush)) {
                        df <- brushedPoints(df, brush,
                                            xvar=names(df)[[2]],
                                            yvar=names(df)[[3]])
                    }

                    return(.prepareDT(df))
                } else {
                    stop("no common identifiers between datasets")
                }
            })
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
        selectizeInput(
            ns("data2"), "Dataset with similar CMap perturbations",
            names(x)[elemClasses == "similarPerturbations"]),
        selectizeInput(ns("col"), "Column to plot in both axes", choices=NULL))
    sidebar[[3]][[2]] <- plotOutput(ns("plot"), brush=ns("brush"))

    mainPanel <- mainPanel(DTOutput(ns("table")))
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
                data1 <- getSelectedDataset1()
                data2 <- getSelectedDataset2()
                if (is.null(data1) || is.null(data2)) return(NULL)
                cols <- intersect(colnames(data1), colnames(data2))
                # Select first ranked column
                rankCol <- grep("rank$", cols, value=TRUE)[1]
                if (is.na(rankCol)) rankCol <- NULL

                updateSelectizeInput(session, "col", choices=cols,
                                     selected=rankCol)
            })

            observe({
                data1 <- getSelectedDataset1()
                data2 <- getSelectedDataset2()
                col <- input$col
                isColValid <- !is.null(col) && col != ""
                if (is.null(data1) || is.null(data2) || !isColValid) {
                    return(NULL)
                }

                plot <- suppressMessages(
                    plotTargetingDrugsVSsimilarPerturbations(
                        data1, data2, column=col, labelBy=NULL) + theme_bw(16))

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

#' @importFrom shiny NS sidebarPanel plotOutput selectizeInput mainPanel
#' tabPanel uiOutput
#' @importFrom DT DTOutput
.drugSetEnrichmentAnalyserUI <- function(id, sets, x,
                                         title="Drug Set Enrichment Analyser") {
    ns <- NS(id)
    sidebar <- sidebarPanel(
        selectizeInput(ns("object"), "Dataset", names(x)),
        selectizeInput(ns("sort"), "Sorting metric", choices=NULL), hr(),
        selectizeInput(ns("statsKey"), "Dataset column to match compounds",
                       choices=NULL),
        selectizeInput(ns("setsKey"), "Drug set column to match compounds",
                       choices=NULL),
        uiOutput(ns("msg")),
        actionButton(ns("analyse"), "Analyse", class="btn-primary"))
    sidebar[[3]][[2]] <- tagList(
        selectizeInput(ns("element"), "Row ID to plot", choices=NULL,
                       width="100%"),
        plotOutput(ns("plot")))

    mainPanel <- mainPanel(DTOutput(ns("table")))
    ui <- tabPanel(title, sidebarLayout(sidebar, mainPanel))
    return(ui)
}

#' @importFrom shiny renderPlot observeEvent observe isolate renderUI
#' @importFrom DT renderDT
.drugSetEnrichmentAnalyserServer <- function(id, sets, x) {
    moduleServer(
        id,
        function(input, output, session) {
            getSelectedObject <- reactive(x[[input$object]])
            getDSEAresult     <- reactive({
                obj      <- getSelectedObject()
                sort     <- input$sort
                statsKey <- input$statsKey
                setsKey  <- input$setsKey

                isValid <- function(e) !is.null(e) && e != ""
                if (is.null(obj) || !isValid(sort)) return(NULL)
                if (!isValid(statsKey) || !isValid(setsKey)) return(NULL)

                analyseDrugSetEnrichment(
                    sets, obj, col=sort,
                    keyColSets=setsKey, keyColStats=statsKey)
            })

            observeEvent(input$object, {
                obj <- getSelectedObject()
                numericCols <- names(obj)[vapply(obj, is.numeric, logical(1))]
                updateSelectizeInput(session, "sort", choices=numericCols)
            })

            # Update available keys to select for datasets
            observe({
                obj <- getSelectedObject()
                if (is.null(obj)) return(NULL)

                statsInfo <- prepareStatsCompoundInfo(obj)$statsInfo
                setsInfo  <- prepareSetsCompoundInfo(sets)$setsCompoundInfo

                probableKey <- findIntersectingCompounds(statsInfo, setsInfo)
                statsKey    <- probableKey$key2
                setsKey     <- probableKey$key1

                keyList      <- getCompoundIntersectingKeyList()
                statsOptions <- intersect(names(statsInfo), keyList)
                setsOptions  <- intersect(names(setsInfo), keyList)

                updateSelectizeInput(session, "statsKey", selected=statsKey,
                                     choices=statsOptions)
                updateSelectizeInput(session, "setsKey", selected=setsKey,
                                     choices=setsOptions)
            })

            # Update number of intersecting compounds based on selected keys
            observe({
                obj <- getSelectedObject()
                if (is.null(obj)) return(NULL)

                statsInfo <- prepareStatsCompoundInfo(obj)$statsInfo
                setsInfo  <- prepareSetsCompoundInfo(sets)$setsCompoundInfo

                statsKey <- input$statsKey
                setsKey  <- input$setsKey
                isValid <- function(e) !is.null(e) && e != ""
                if (!isValid(statsKey) || !isValid(setsKey)) return(NULL)

                probableKey <- findIntersectingCompounds(statsInfo, setsInfo,
                                                         statsKey, setsKey)
                num <- length(probableKey[[3]])
                msg <- "cross-matches found using the selected keys"
                output$msg <- renderUI(helpText(paste(num, msg)))
            })

            observeEvent(input$analyse, {
                dsea <- getDSEAresult()
                updateSelectizeInput(session, "element", choices=dsea[[1]])
                output$table <- renderDT(.prepareDT(dsea))

                output$plot <- renderPlot({
                    obj <- getSelectedObject()
                    if (is.null(obj)) return(NULL)

                    element <- input$element
                    if (element == "") element <- NULL
                    if (is.null(element) || !element %in% names(sets)) {
                        return(NULL)
                    }

                    isolate({
                        setsKey  <- input$setsKey
                        statsKey <- input$statsKey
                        sort     <- input$sort
                    })
                    isValid <- function(e) !is.null(e) && e != ""
                    if (!isValid(statsKey) || !isValid(setsKey)) return(NULL)
                    plotDrugSetEnrichment(sets, obj, col=sort,
                                          selectedSets=element,
                                          keyColStats=statsKey,
                                          keyColSets=setsKey)[[1]]
                })
            })

            # Update selected element
            observeEvent(input$table_rows_selected, {
                selected <- getDSEAresult()[[1]][input$table_rows_selected]
                updateSelectizeInput(session, "element", selected=selected)
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
        .dataPlotterUI(dataId, elems),
        .targetingDrugsVSsimilarPerturbationsPlotterUI(
            comparePlotId, elems, elemClasses),
        .datasetComparisonUI(compareId, elems),
        .metadataViewerUI(metadataId, elems))
    uiList <- Filter(length, uiList)
    ui     <- do.call(.prepareNavPage, uiList)

    server <- function(input, output, session) {
        .dataPlotterServer(dataId, elems)
        if (showTwoKindPlot) {
            .targetingDrugsVSsimilarPerturbationsPlotterServer(
                comparePlotId, elems)
        }
        .datasetComparisonServer(compareId, elems)
        .metadataViewerServer(metadataId, elems)
    }
    app <- runApp(shinyApp(ui, server))
    return(app)
}

#' View and plot results via a visual interface
#'
#' @inheritParams analyseDrugSetEnrichment
#' @param ... Objects
#'
#' @importFrom shiny tagList
#'
#' @return Launches result viewer and plotter (returns \code{NULL})
#' @family visual interface functions
#' @export
launchDrugSetEnrichmentAnalyser <- function(sets, ...) {
    dataId     <- "dataPlotter"
    metadataId <- "metadataViewer"
    dseaId     <- "drugSetAnalyser"

    elems <- .prepareEllipsis(...)

    uiList <- tagList(
        .drugSetEnrichmentAnalyserUI(dseaId, sets, elems),
        .dataPlotterUI(dataId, elems),
        .metadataViewerUI(metadataId, elems))
    uiList <- Filter(length, uiList)
    ui     <- do.call(.prepareNavPage, uiList)

    server <- function(input, output, session) {
        .drugSetEnrichmentAnalyserServer(dseaId, sets, elems)
        .dataPlotterServer(dataId, elems)
        .metadataViewerServer(metadataId, elems)
    }
    app <- runApp(shinyApp(ui, server))
    return(app)
}
