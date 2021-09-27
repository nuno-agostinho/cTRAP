# Skeleton for common elements -------------------------------------------------

.convertToFunction <- function(x) {
    # Easier to deal with reactives
    res <- function() x
    if (is.reactive(x)) res <- x
    return(res)
}

.filterDatasetsByClass <- function(dataset, class) {
    selected <- sapply(dataset, is, class)
    if (!any(selected)) return(NULL)
    return(dataset[selected])
}

#' @importFrom shiny tags
.alert <- function(..., type=c("danger", "success", "info", "warning"),
                   condition=NULL, ns=NULL) {
    type <- match.arg(type)
    alert <- tags$div(class=paste0("alert alert-", type), role="alert", ...)
    if (!is.null(condition)) alert <- conditionalPanel(condition, alert, ns=ns)
    return(alert)
}

.prepareDiffExprDataset <- function(data, name) {
    attr(data, "name") <- name
    class(data) <- c("diffExpr", class(data))
    return(data)
}

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
    if (nrow(df) == 0) {
        hc <- highchart()
    } else {
        names(df) <- c("data", "Freq")
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
            lenUniqValuesPerCol <- sapply(lapply(table, unique), length)
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
    if (is.null(x)) return(NULL)
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
        
        filter <- attr(attr(x, "metadata"), "filter")
        if (!is.null(filter)) {
            key <- c("cellLine"="cell line",
                     "timepoint"="time point",
                     "dosage"="dosage",
                     "perturbationType"="perturbation type")
            key <- paste("Filtered by", key[names(filter)])
            filter <- cbind(key, lapply(filter, paste, collapse=", "))
            data <- rbind(data, filter)
        }
        
        colnames(data) <- c("Key", "Value")
        dataName <- "Object summary"
    } else {
        data <- tryCatch(as.table(x, clean=FALSE), error=return)
        if (is(data, "error") || !is.data.frame(data)) {
            if (!is.null(names(x))) {
                data <- data.frame("Element"=names(x), "Value"=x)
            } else {
                data <- data.frame("Value"=x)
            }
        }
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
    if (length(elems) == 0 && length(names(elems)) == 0) elems <- NULL
    return(elems)
}

# Data input -------------------------------------------------------------------

.getENCODEconditions <- function(metadata) {
    cellLines <- sort(unique(metadata$`Biosample term name`))
    genes     <- sort(unique(metadata$`Experiment target`))
    genes     <- genes[genes != ""]
    
    res <- list(genes=genes, cellLines=cellLines)
    return(res)
}

#' @importFrom shiny textAreaInput
.diffExprLoadUI <- function(id, title="User data") {
    ns <- NS(id)
    sidebar <- sidebarPanel(
        textAreaInput(ns("diffExpr"), "Differential gene expression",
                      height="200px"),
        selectizeInput(
            ns("sep"), "Separator",
            c("Automatic"="auto", "Tab (\\t)"="\t", "Space"=" ", ",", ";")),
        textInput(ns("name"), "Dataset name", "Differential expression"),
        actionButton(ns("load"), "Load data", class="btn-primary"))
    mainPanel <- mainPanel(
        .alert("No differential gene expression dataset loaded/selected",
               condition="input.dataset == ''", ns=ns),
        selectizeInput(ns("dataset"), "Differential gene expression dataset",
                       choices=NULL, width="100%"),
        DTOutput(ns("table")))
    ui <- tabPanel(title, sidebarLayout(sidebar, mainPanel))
    return(ui)
}

#' @importFrom data.table fread
#' @importFrom DT renderDT
.diffExprLoadServer <- function(id, x) {
    x <- .convertToFunction(x)
    server <- function(input, output, session) {
        loadData <- eventReactive(input$load, {
            diffExpr <- input$diffExpr
            req(diffExpr)
            
            # Remove comments
            diffExpr <- gsub("[(^|\n)]{0,1}#.*?\n", "\n", diffExpr)
            diffExpr <- gsub("\n+", "\n", diffExpr)
            
            withProgress(message="Loading differential expression", value=2, {
                data <- fread(text=diffExpr, sep=input$sep, data.table=FALSE,
                              header=FALSE)
                incProgress(1)
                
                if (ncol(data) == 1) {
                    data <- data[[1]]
                } else if (ncol(data) == 2) {
                    data <- setNames(data[[2]], data[[1]])
                } else {
                    showNotification(
                        tagList(tags$b("Input requires 1 or 2 columns."),
                                ncol(data), "columns are not supported."),
                        type="error")
                    return(NULL)
                }
                data <- .prepareDiffExprDataset(data, input$name)
                incProgress(2)
            })
            return(data)
        })
        
        observe({
            choices <- names(x())[sapply(x(), is, "diffExpr")]
            if (is.null(choices)) choices <- list()
            selected <- choices[length(choices)]
            updateSelectizeInput(session, "dataset", choices=choices,
                                 selected=selected)
        })
        
        output$table <- renderDT({
            req(input$dataset)
            data <- x()[[input$dataset]]
            if (!is.null(names(data))) {
                df <- data.frame("Genes"=names(data), "Score"=data)
            } else {
                df <- data.frame("Genes"=data)
            }
            return(df)
        }, rownames=FALSE)
        return(loadData)
    }
    moduleServer(id, server)
}

#' @importFrom shiny NS sidebarPanel selectizeInput mainPanel tabPanel
#' sidebarLayout
#' @importFrom DT DTOutput
.diffExprENCODEloaderUI <- function(id,
                                    metadata=downloadENCODEknockdownMetadata(),
                                    cellLine=NULL, gene=NULL,
                                    title="ENCODE knockdown data") {
    ns <- NS(id)
    conditions <- .getENCODEconditions(metadata)
    cellLines  <- conditions$cellLines
    genes      <- conditions$genes
    
    nullStrJS <- I('function() { this.setValue(""); }')
    onInitializeCellLine <- NULL
    if (is.null(cellLine)) onInitializeCellLine <- nullStrJS
    onInitializeGene <- NULL
    if (is.null(gene)) onInitializeGene <- nullStrJS
    
    ui <- tabPanel(
        title,
        div(class="well", style="padding-bottom: 9px;",
            fluidRow(
                column(4, selectizeInput(
                    ns("cellLine"), "Cell line", cellLines, cellLine, 
                    width="100%",
                    options=list(placeholder='Select a cell line',
                                 onInitialize=onInitializeCellLine))),
                column(4, selectizeInput(
                    ns("gene"), "Gene", genes, gene, width="100%",
                    options=list(placeholder='Select a gene',
                                 onInitialize=onInitializeGene))),
                column(4, actionButton(ns("load"), "Load data",
                                       style="position:absolute; top:25px;",
                                       class="btn-primary")))),
        DTOutput(ns("table")))
    return(ui)
}

#' @importFrom shiny moduleServer stopApp withProgress incProgress
#' @importFrom DT renderDT
.diffExprENCODEloaderServer <- function(
    id, metadata=downloadENCODEknockdownMetadata(), shinyproxy=FALSE) {
    server <- function(input, output, session) {
        output$table <- renderDT(.prepareDT(metadata))
        proxy <- dataTableProxy("table")
        
        observe({
            cellLine <- input$cellLine
            if (is.null(cellLine) || cellLine == "") cellLine <- NULL
            
            gene <- input$gene
            if (is.null(gene) || gene == "") gene <- NULL
            data <- downloadENCODEknockdownMetadata(cellLine, gene)
            observe(replaceData(proxy, data, rownames=FALSE))
        })
        
        loadData <- eventReactive(input$load, {
            withProgress(message="Preparing differential expression", {
                steps <- 3
                
                incProgress(1/steps, detail="Downloading")
                message("Downloading data...")
                ENCODEmetadata <- downloadENCODEknockdownMetadata(
                    input$cellLine, input$gene)
                
                if (nrow(ENCODEmetadata) == 0) {
                    showMessage("No samples match the selected criteria",
                                type="error")
                    return(NULL)
                }
                ENCODEsamples <- loadENCODEsamples(ENCODEmetadata)[[1]]
                counts <- prepareENCODEgeneExpression(ENCODEsamples)
                
                # Remove low coverage (>= 10 counts shared by 2 samples)
                minReads   <- 10
                minSamples <- 2
                filter <- rowSums(
                    counts[ , -c(1:2)] >= minReads) >= minSamples
                counts <- counts[filter, ]
                
                # Convert ENSEMBL identifier to gene symbol
                msg <- "Converting ENSEMBL identifiers to gene symbols"
                incProgress(1/steps, detail=msg)
                message(paste0(msg, "..."))
                counts$gene_id <- convertGeneIdentifiers(counts$gene_id)
                
                # Perform differential gene expression analysis
                msg <- "Performing differential gene expression analysis"
                incProgress(1/steps, detail=msg)
                message(paste0(msg, "..."))
                diffExpr <- performDifferentialExpression(counts)
                
                diffExprStat <- diffExpr$t
                names(diffExprStat) <- diffExpr$Gene_symbol
                
                msg <- sprintf("%s vs control DGE in %s (ENCODE)",
                               input$gene, input$cellLine)
                diffExprStat <- .prepareDiffExprDataset(diffExprStat, msg)
                
                message(paste(msg, "data loaded"))
                return(diffExprStat)
            })
        })
        
        if (!shinyproxy) observeEvent(input$load, stopApp(loadData()))
        return(loadData)
    }
    moduleServer(id, server)
}

.findMatch <- function(what, where, ignore.case=TRUE) {
    if (is.null(what)) return(NULL)
    if (is.list(where)) where <- unlist(where)
    what <- paste0("^", what, "$") # exact match
    unname(sapply(what, grep, where, ignore.case=ignore.case, value=TRUE))
}

#' @importFrom shiny NS selectizeInput checkboxGroupInput sidebarPanel helpText
#' @importFrom shinycssloaders withSpinner
#' @importFrom htmltools tagQuery
.cmapDataLoaderUI <- function(id, cellLine=NULL, timepoint=NULL, dosage=NULL,
                              perturbationType=NULL, title="CMap perturbations",
                              shinyproxy=FALSE) {
    ns <- NS(id)
    dataTypes <- c("Perturbation metadata"="metadata",
                   "Perturbation z-scores"="zscores",
                   "Gene information"="geneInfo",
                   "Compound information"="compoundInfo")
    selectizeCondition <- function(id, label, choices, selected, ...) {
        selected    <- .findMatch(selected, choices)
        plugins     <- list("remove_button")
        placeholder <- paste("Select one or more", tolower(label))
        return(selectizeInput(id, label, choices, selected, ...,
                              options=list(plugins=plugins,
                                           placeholder=placeholder)))
    }
    
    colChart <- function(id) column(3, highchartOutput(id, height="200px"))
    plots <- fluidRow(
        colChart(ns("pertPlot")),
        colChart(ns("cellLinePlot")),
        colChart(ns("dosagePlot")),
        colChart(ns("timepointPlot")))
    
    dataToLoad <- checkboxGroupInput(ns("data"), "Data to load",
                                     dataTypes, dataTypes)
    if (shinyproxy) {
        # Disable check boxes
        dataToLoad <- tagQuery(dataToLoad)
        dataToLoad <- dataToLoad$find("input")$addAttrs("disabled"=NA)$allTags()
    }
    
    sidebar <- sidebarPanel(
        selectizeCondition(ns("type"), "Perturbation types", multiple=TRUE,
                           choices=perturbationType, perturbationType),
        selectizeCondition(ns("cellLine"), "Cell lines", multiple=TRUE,
                           choices=cellLine, cellLine),
        selectizeCondition(ns("dosage"), "Dosages", multiple=TRUE,
                           choices=dosage, dosage),
        selectizeCondition(ns("timepoint"), "Time points", multiple=TRUE,
                           choices=timepoint, timepoint),
        dataToLoad,
        if (shinyproxy) textInput(ns("name"), "Dataset name", value="cmapData"),
        if (!shinyproxy)
            helpText("By default, data will be downloaded if not found."),
        if (!shinyproxy) actionButton(ns("cancel"), "Cancel"),
        actionButton(ns("load"), "Load data", class="btn-primary"))
    mainPanel <- mainPanel( withSpinner(DTOutput(ns("table")), type=6) )
    
    ui <- tabPanel(title, sidebar, mainPanel)
    ui[[3]] <- c(list(plots), ui[[3]])
    return(ui)
}

#' @importFrom shiny isolate updateSelectizeInput moduleServer stopApp
#' observeEvent eventReactive
#' @importFrom DT dataTableProxy replaceData
.cmapDataLoaderServer <- function(id, metadata="cmapMetadata.txt",
                                  zscores="cmapZscores.gctx",
                                  geneInfo="cmapGeneInfo.txt",
                                  compoundInfo="cmapCompoundInfo.txt",
                                  cellLine=NULL, timepoint=NULL, dosage=NULL,
                                  shinyproxy=FALSE, tab=NULL) {
    server <- function(input, output, session) {
        updateSelectizeCondition <- function(session, input, id, choices, ...) {
            selected <- isolate(input[[id]])
            selected <- .findMatch(selected, choices)
            return(updateSelectizeInput(session, id, choices=choices,
                                        selected=selected, ...))
        }
        loadCMapMetadata <- reactive({
            withProgress(message="Loading CMap metadata", {
                metadata <- loadCMapData(metadata, "metadata")
                incProgress(1)
                return(metadata)
            })
        })
        
        happened <- reactiveVal(FALSE)
        checkTab <- reactive({
            load <- TRUE
            if (!isolate(happened())) {
                if (is.null(tab)) {
                    load <- TRUE
                } else {
                    tab  <- tab()
                    load <- !is.null(tab) && grepl(
                      "CMap", tab, ignore.case=TRUE)
                }
            }
            happened(load)
            if (!load) return(NULL)
            return(load)
        })
        
        getCMapMetadata <- eventReactive(checkTab(), {
            if (!is.data.frame(metadata)) metadata <- loadCMapMetadata()
            return(metadata)
        })
        
        # Update conditions based on selected perturbation type
        observe({
            isolate({
                pertType  <- input$type
                cellLine  <- input$cellLine
                dosage    <- input$dosage
                timepoint <- input$timepoint
            })
            
            metadata <- getCMapMetadata()
            req(metadata)
            available <- getCMapConditions(metadata, control=TRUE)
            
            # Prepare perturbation types
            perts <- available$perturbationType
            controls <- startsWith(perts, "Controls")
            choices <- list("Perturbations"=perts[!controls],
                            "Controls"=perts[controls])
            updateSelectizeCondition(session, input, "type", choices=choices)
            updateSelectizeCondition(session, input, "cellLine",
                                     choices=available$cellLine)
            updateSelectizeCondition(session, input, "dosage",
                                     choices=available$dosage)
            updateSelectizeCondition(session, input, "timepoint",
                                     choices=available$timepoint)
        })
        
        # Filter metadata based on selected inputs
        getFilteredMetadata <- reactive({
            metadata <- getCMapMetadata()
            req(metadata)
            filterCMapMetadata(
                metadata, perturbationType=input$type, cellLine=input$cellLine,
                dosage=input$dosage, timepoint=input$timepoint)
        })
        
        renderBubbleChart <- function(subset, title, colour) {
            renderHighchart({
                filt <- getFilteredMetadata()
                req(filt)
                if (nrow(filt) != 0) {
                    data <- filt[[subset]]
                    if (subset == "pert_type") {
                        label <- getCMapPerturbationTypes(control=TRUE)
                        data  <- table(data)
                        names(data) <- names(label)[match(names(data), label)]
                    } else {
                        data[is.na(data)] <- "NA"
                    }
                } else {
                    data <- NULL
                }
                .plotBubbles(data, title, colour)
            })
        }
        output$cellLinePlot  <- renderBubbleChart("cell_id", "Cell lines",
                                                  "orange")
        output$dosagePlot    <- renderBubbleChart("pert_idose", "Dosages",
                                                  "green")
        output$timepointPlot <- renderBubbleChart("pert_itime", "Time points",
                                                  "purple")
        output$pertPlot      <- renderBubbleChart("pert_type", "Perturbations",
                                                  "red")
        
        # Show table
        output$table <- renderDT({
            hiddenCols <- c("pert_dose", "pert_dose_unit", "pert_time",
                            "pert_time_unit")
            hiddenCols <- match(hiddenCols, colnames(getCMapMetadata()))
            columnDefs <- list(list(visible=FALSE, targets=hiddenCols - 1))
            .prepareDT(isolate(getFilteredMetadata()), columnDefs=columnDefs,
                       scrollX=TRUE)
        }, server=TRUE)
        proxy <- dataTableProxy("table")
        observe(replaceData(proxy, getFilteredMetadata(), rownames=FALSE))
        
        # Load data
        loadData <- eventReactive(input$load, {
            types <- input$data
            
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
            
            withProgress(message="Loading CMap perturbations", {
                perturbations <- prepareCMapPerturbations(
                    metadata=metadata, zscores=zscores,
                    geneInfo=geneInfo, compoundInfo=compoundInfo)
                attr(perturbations, "name") <- input$name
                incProgress(1)
            })
            return(perturbations)
        })
        
        observe({
            txt <- "CMap"
            
            type <- input$type
            kd <- "Consensus signature from shRNAs targeting the same gene"
            oe <- "cDNA for overexpression of wild-type gene"
            
            if (length(input$cellLine) == 1) txt <- paste(txt, input$cellLine)
            if (length(type) == 1) {
                txtType <- NULL
                if (type == "Compound") txtType <- "compounds"
                if (type == kd) txtType <- "knockdowns"
                if (type == oe) txtType <- "overexpressions"
                txt <- paste(txt, txtType)
            }
            
            if (txt == "cmap") txt <- "CMap data"
            updateTextInput(session, "name", value=txt)
        })
        
        if (!shinyproxy) observeEvent(input$load, stopApp(loadData()))
        observeEvent(input$cancel, stopApp(stop("User cancel", call.=FALSE)))
        return(loadData)
    }
    
    moduleServer(id, server)
}

# Result plotting and viewing --------------------------------------------------

#' @importFrom shiny NS sidebarPanel mainPanel tabPanel sidebarLayout
#' selectizeInput fluidRow column
.metadataViewerUI <- function(id, title="Metadata") {
    ns <- NS(id)
    ui <- tabPanel(
        title,
        div(class="well", style="padding-bottom: 9px;",
            fluidRow(
                column(4, selectizeInput(ns("object"), "Dataset", choices=NULL, 
                                         width="100%")),
                column(4, selectizeInput(ns("attr"), "Table", choices=NULL, 
                                         width="100%")))),
        DTOutput(ns("table")))
    return(ui)
}

#' @importFrom shiny moduleServer observe
#' @importFrom DT renderDT
.metadataViewerServer <- function(id, x) {
    x <- .convertToFunction(x)
    moduleServer(
        id,
        function(input, output, session) {
            getSelectedObject <- reactive(.prepareMetadata(x()[[input$object]]))
            
            observe(updateSelectizeInput(session, "object", choices=names(x())))
            
            observe({
                selected <- isolate(input$attr)
                choices  <- names(getSelectedObject())
                anyChoiceSelected <- !is.null(selected) && !is.null(choices) &&
                    selected %in% choices
                if (!anyChoiceSelected) selected <- NULL
                updateSelectizeInput(session, "attr", selected=selected,
                                     choices=choices)
            })
            
            output$table <- renderDT(
                .prepareDT(getSelectedObject()[[input$attr]]))
        }
    )
}

#' @importFrom shiny NS sidebarPanel plotOutput selectizeInput mainPanel
#' tabPanel
#' @importFrom DT DTOutput
.dataPlotterUI <- function(id, x, title="Plot Data") {
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

.rankSimilarityMethods <- function() {
    c("Spearman's correlation"="spearman",
      "Pearson's correlation"="pearson",
      "GSEA"="gsea")
}

#' @importFrom shiny renderPlot observeEvent observe
#' @importFrom DT renderDT formatSignif dataTableProxy replaceData
.dataPlotterServer <- function(id, x) {
    x <- .convertToFunction(x)
    isValid <- function(e) !is.null(e) && e != ""
    moduleServer(
        id,
        function(input, output, session) {
            getSelectedObject <- reactive(x()[[input$object]])
            
            observe({
                obj <- x()
                if (length(obj) == 0) return(NULL)
                ref <- sapply(obj, is, "referenceComparison")
                if (!any(ref)) return(NULL)
                choices <- names(obj[ref])
                updateSelectizeInput(session, "object", choices=choices,
                                     server=TRUE)
            })
            
            # Update element and methods choices depending on selected object
            observeEvent(input$object, {
                obj <- getSelectedObject()
                if (is.null(obj)) return(NULL)
                
                updateSelectizeInput(session, "element", choices=obj[[1]],
                                     selected=list(), server=TRUE)
                # Update methods depending on selected object
                methods <- .rankSimilarityMethods()
                isPresent <- sapply(methods, function(i)
                    any(grepl(i, colnames(obj), ignore.case=TRUE)))
                methods <- methods[isPresent]
                
                selected <- isolate(input$method)
                if (!selected %in% methods) selected <- NULL
                updateSelectizeInput(session, "method", choices=methods,
                                     selected=selected)
            })
            
            # Update selected element
            observe({
                obj <- getSelectedObject()
                if (is.null(obj) || !is(obj, "referenceComparison")) {
                    return(NULL)
                }
                
                selected <- obj[[1]][input$table_rows_selected]
                updateSelectizeInput(session, "element", selected=selected)
            })
            
            output$table <- renderDT({
                obj <- getSelectedObject()
                if (is.null(obj) || !is(obj, "referenceComparison")) {
                    return(NULL)
                }
                data <- as.table(obj, clean=FALSE)
                .prepareDT(data)
            })
            
            # Filter table based on overall plot
            proxy <- dataTableProxy("table")
            observe({
                elem <- input$element
                if (is.null(elem) || elem == "") {
                    obj <- getSelectedObject()
                    if (is.null(obj) || !is(obj, "referenceComparison")) {
                        return(NULL)
                    }
                    
                    obj   <- as.table(obj, clean=FALSE)
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
                if (is.null(obj) || !is(obj, "referenceComparison")) {
                    return(NULL)
                }
                
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
    x <- .convertToFunction(x)
    getNumericCols <- function(x) colnames(x)[vapply(x, is.numeric, logical(1))]
    
    moduleServer(
        id,
        function(input, output, session) {
            getSelectedDataset1 <- reactive(x()[[input$data1]])
            getSelectedDataset2 <- reactive(x()[[input$data2]])
            
            observe({
                obj <- x()
                req(names(obj))
                choices <- names(obj)[sapply(obj, is, "referenceComparison")]
                req(choices)
                if (length(choices) < 2) return(NULL)
                
                updateSelectizeInput(session, "data1", choices=choices)
                updateSelectizeInput(session, "data2", choices=choices,
                                     selected=choices[[2]])
            })
            
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
    id, title="Dataset Comparison 2") {
    
    if (!all(c("targetingDrugs", "similarPerturbations") %in% elemClasses)) {
        return(NULL)
    }
    ns <- NS(id)
    sidebar <- sidebarPanel(
        selectizeInput(
            ns("data1"), "Dataset with predicted targeting drugs", NULL),
        selectizeInput(
            ns("data2"), "Dataset with similar CMap perturbations", NULL),
        selectizeInput(ns("col"), "Column to plot in both axes", choices=NULL))
    sidebar[[3]][[2]] <- plotOutput(ns("plot"), brush=ns("brush"))
    
    mainPanel <- mainPanel(DTOutput(ns("table")))
    ui <- tabPanel(title, sidebarLayout(sidebar, mainPanel))
    return(ui)
}

#' @importFrom shiny moduleServer reactive observe renderPlot brushedPoints
#' @importFrom DT renderDT
.targetingDrugsVSsimilarPerturbationsPlotterServer <- function(id, x) {
    x <- .convertToFunction(x)
    moduleServer(
        id,
        function(input, output, session) {
            getDataset <- function(obj, item) {
                if (length(obj) > 0 && !is.null(item)) obj[[item]]
            }
            getSelectedDataset1 <- reactive(getDataset(x(), input$data1))
            getSelectedDataset2 <- reactive(getDataset(x(), input$data2))
            
            observe({
                obj <- x()
                elemClasses <- sapply(lapply(obj, class), "[[", 1)
                updateSelectizeInput(
                    session, "object", server=TRUE,
                    names(obj)[elemClasses == "targetingDrugs"])
                updateSelectizeInput(
                    session, "object", server=TRUE,
                    choices=names(obj)[elemClasses == "similarPerturbations"])
            })
            
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
                if (is.null(data1) || is.null(data2)) return(NULL)
                browser()
                
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

# Data analysis ----------------------------------------------------------------

#' @importFrom shiny NS sidebarPanel plotOutput selectizeInput mainPanel
#' tabPanel uiOutput column fluidRow numericInput checkboxInput conditionalPanel
#' @importFrom DT DTOutput
.rankSimilarPerturbationsUI <- function(
    id, diffExpr, perturbations,
    title="Rank CMap perturbations by similarity") {
    
    ns <- NS(id)
    sidebar <- sidebarPanel(
        selectizeInput(ns("diffExpr"), "Differential gene expression",
                       names(diffExpr)),
        selectizeInput(ns("perts"), "CMap perturbations", names(perturbations)),
        selectizeInput(ns("method"), "Method", multiple=TRUE,
                       .rankSimilarityMethods(), .rankSimilarityMethods(),
                       options=list(plugins=list("remove_button"))),
        conditionalPanel(
            "input.method.includes('gsea')",
            fluidRow(
                column(6, numericInput(ns("upGenes"), "Top genes", 150)),
                column(6, numericInput(ns("downGenes"), "Bottom genes", 150))),
            ns=ns),
        selectizeInput(ns("cellLineMean"),
                       "Calculate mean across cell lines",
                       c("For data with ≥ 2 cell lines"="auto",
                         "Always"=TRUE,
                         "Never"=FALSE)),
        conditionalPanel(
          "input.cellLineMean != 'FALSE'",
          selectizeInput(ns("rankPerCellLine"), "Rank results based on",
                         c("Mean scores only"=FALSE,
                           "Mean + individual cell lines' scores"=TRUE)),
          ns=ns),
        textInput(ns("name"), "Dataset name", "cmap_ranking"),
        uiOutput(ns("msg")),
        actionButton(ns("analyse"), "Rank by similarity", class="btn-primary"))
    
    mainPanel <- mainPanel(DTOutput(ns("table")))
    ui <- tabPanel(title, sidebarLayout(sidebar, mainPanel))
    return(ui)
}

#' @importFrom shiny renderPlot observeEvent observe isolate renderUI
#' @importFrom DT renderDT
#' @importFrom data.table rbindlist
.rankSimilarPerturbationsServer <- function(id, diffExpr, perturbations,
                                            shinyproxy=FALSE) {
    server <- function(input, output, session) {
        getSelectedObject <- reactive(perturbations[[input$object]])
        
        updateDatasetChoices <- function(id, var, class=NULL) {
            req(names(var))
            if (!is.null(class)) {
                ref <- sapply(var, is, class)
            } else {
                ref <- TRUE
            }
            if (!any(ref)) return(NULL)
            choices <- names(var[ref])
            updateSelectizeInput(session, id, choices=choices)
        }
        
        observe(updateDatasetChoices(
            "diffExpr", diffExpr(), "diffExpr"))
        observe(updateDatasetChoices(
            "perts", perturbations(), "perturbationChanges"))
        
        rankData <- eventReactive(input$analyse, {
            diffExprDataset <- req(input$diffExpr)
            pertsDataset    <- req(input$perts)
            method          <- input$method
            upGenes         <- input$upGenes
            downGenes       <- input$downGenes
            cellLineMean    <- input$cellLineMean
            rankPerCellLine <- input$rankPerCellLine
            dataset         <- input$name
            
            selectedDiffExpr <- diffExpr()[[req(diffExprDataset)]]
            selectedPerts    <- perturbations()[[req(pertsDataset)]]
            
            withProgress(message="Ranking against CMap perturbations", {
                ranking <- rankSimilarPerturbations(
                    selectedDiffExpr, selectedPerts, method, c(upGenes, downGenes),
                    cellLineMean, rankPerCellLine)
                incProgress(1)
            })
            attr(ranking, "name") <- dataset
            
            # Prepare form input
            cellLineMean <- switch(cellLineMean,
                                   "TRUE"="Always", "FALSE"="Never",
                                   "auto"="For data with ≥ 2 cell lines")
            rankPerCellLine <- ifelse(
                rankPerCellLine,
                "Mean + individual cell lines' scores",
                "Mean scores only")
            attr(ranking, "formInput") <- list(
                "Differential expression dataset"=diffExprDataset,
                "CMap perturbation dataset"=pertsDataset,
                "Methods"=paste(method, collapse=", "),
                "Top genes"=upGenes,
                "Bottom genes"=downGenes,
                "Calculate mean across cell lines"=cellLineMean,
                "Rank results based on"=rankPerCellLine)
            return(ranking)
        })
        
        if (!shinyproxy) observeEvent(input$load, stopApp(rankData()))
        
        output$table <- renderDT({
            data <- .filterDatasetsByClass(diffExpr(), "similarPerturbations")
            req(data)
            
            formInput <- lapply(data, attr, "formInput")
            ns <- unique(unlist(lapply(formInput, names)))
            
            # Add name this way to avoid issues with data from outside the UI
            for (i in seq(formInput)) formInput[[i]]$Dataset <- names(data)[[i]]
            res <- rbindlist(formInput, fill=TRUE)
            res <- cbind("Dataset"=names(data), "Progress"="Loaded", res)
            return(res)
        }, rownames=FALSE)
        
        return(rankData)
    }
    moduleServer(id, server)
}

#' @importFrom shiny NS sidebarPanel plotOutput selectizeInput mainPanel
#' tabPanel uiOutput
#' @importFrom DT DTOutput
.drugSetEnrichmentAnalyserUI <- function(id, sets, x,
                                         title="Drug Set Enrichment") {
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
            
            observe({
                req(names(x))
                ref <- sapply(x, is, "referenceComparison")
                if (!any(ref)) return(NULL)
                choices <- names(x[ref])
                updateSelectizeInput(session, "object", choices=choices)
            })
            
            getDSEAresult <- reactive({
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
                
                if (!is(sets, "drugSets")) return(NULL)
                
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
        .cmapDataLoaderUI(id, cellLine, timepoint, dosage, perturbationType))
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
    ui     <- .prepareNavPage(.metadataViewerUI(id))
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
        .targetingDrugsVSsimilarPerturbationsPlotterUI(comparePlotId),
        .datasetComparisonUI(compareId, elems),
        .metadataViewerUI(metadataId))
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
        .metadataViewerUI(metadataId))
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
