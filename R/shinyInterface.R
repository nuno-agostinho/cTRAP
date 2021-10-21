# Skeleton for common elements -------------------------------------------------

#' @importFrom shiny is.reactive
.convertToFunction <- function(x) {
    # Easier to deal with reactives
    res <- function() x
    if (is.reactive(x)) res <- x
    return(res)
}

.filterDatasetsByClass <- function(data, class, expected=FALSE) {
    selected <- sapply(data, is, class)
    if (expected) {
        # Return values that are expected to turn into the given class
        expected <- sapply(data, is, paste0("expected", capitalize(class)))
        if (is.list(expected) && length(expected) == 0) return(NULL)
        selected <- selected | expected
    }
    if (!any(selected)) return(NULL)
    return(data[selected])
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

# Replace a string with another in a list
.replaceStrInList <- function(tag, old, new) {
    FUN <- function(x) {
        res <- x
        if (grepl(old, x)) res <- gsub(old, new, x, fixed=TRUE)
        return(res)
    }
    rapply(tag, FUN, how="replace", classes="character")
}

#' Prepare Shiny page template
#'
#' @importFrom shiny navbarPage tabPanel sidebarLayout tags
#' @importFrom highcharter %>%
#'
#' @return HTML elements
#' @keywords internal
.prepareNavPage <- function(...) {
    app <- "cTRAP"
    ui  <- navbarPage(app, ...) %>%
        .replaceStrInList("navbar-static-top", "") %>%
        .replaceStrInList("container-fluid", "") %>%
        tags$div(class="container-fluid", style="padding-top: 15px;")
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

.prepareReferenceComparisonDT <- function(data, class) {
    data <- .filterDatasetsByClass(req(data), class, expected=TRUE)
    req(data)
    
    formInput <- lapply(data, attr, "formInput")
    ns <- unique(unlist(lapply(formInput, names)))
    
    # Add name to avoid issues with data from outside the shiny UI
    for (i in seq(formInput)) formInput[[i]]$Dataset <- names(data)[[i]]
    res <- rbindlist(formInput, fill=TRUE)
    
    # Return task state if available
    state     <- sapply(data, getTaskState)
    isLoaded  <- state == "Loaded"
    stateHTML <- sapply(state, convertTaskState2HTML)
    
    # Create link to data results
    dataset   <- names(data)
    js        <- paste("$(\"a[data-value*='Plot']\").click(); ",
                       "$('#dataPlotter-object')[0].selectize.setValue('%s');")
    datasetJS <- as.character(tags$a(href="#", onclick=js, "%s"))
    datasetJS <- sprintf(datasetJS, dataset, dataset)
    dataset   <- ifelse(isLoaded, datasetJS, dataset)
    
    res <- cbind("Dataset"=dataset, "Progress"=stateHTML, res)
    # Reverse order (so newest datasets are on top)
    res <- res[rev(seq(nrow(res))), unique(colnames(res)), with=FALSE]
    return(res)
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
        data <- tryCatch(as.table(x, clean=FALSE), error=function(e) e)
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

# Add tags to datasets in "display" attribute
.addDatasetTags <- function(elems) {
    # Only edit data with no "display" attribute
    noDisplay <- sapply(lapply(elems, attr, "display"), is.null)
    if (!any(noDisplay)) return(elems)
    
    getDatasetTags <- function(name, data) {
        dataset <- data[[name]]
        if (is.null(dataset)) return(dataset)
        class  <- class(dataset)[[1]]
        source <- attr(dataset, "source")
        tags   <- paste0("#", c(class, source), collapse=" ")
        
        html   <- convertTaskState2HTML(getTaskState(dataset), label=TRUE)
        attr(dataset, "display") <- paste(name, tags, html)
        return(dataset)
    }
    dataset          <- names(elems[noDisplay])
    elems[noDisplay] <- lapply(dataset, getDatasetTags, elems[noDisplay])
    return(elems)
}

# Update dataset choices (optionally, filter datasets by class)
.updateDatasetChoices <- function(session, id, data, class=NULL) {
    if (!is.null(class)) data <- .filterDatasetsByClass(data, class)
    if (is.null(data) || length(data) == 0) {
        choices <- list()
    } else {
        data    <- .addDatasetTags(data)
        choices <- setNames(names(data), sapply(data, attr, "display"))
    }
    render  <- I('{ option: renderSelectizeTags, item: renderSelectizeTags }')
    
    # Keep previous selection if possible
    selected <- isolate(session$input[[id]])
    if (is.null(selected) || !selected %in% choices) selected <- NULL
    
    updateSelectizeInput(session, id, choices=choices, selected=selected,
                         options=list(render=render))
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
    
    diffExpr <- tags$span(
        "data-toggle"="tooltip", "data-placement"="right",
        title="Gene symbols and respective differential expression values (e.g. t-statistics)",
        "Differential gene expression", icon("question-circle"))
    
    sidebar <- sidebarPanel(
        textAreaInput(ns("diffExpr"), diffExpr, height="300px"),
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
        
        observe(.updateDatasetChoices(session, "dataset", x(), "diffExpr"))
        
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
.diffExprENCODEloaderUI <- function(id, title="ENCODE knockdown data") {
    ns <- NS(id)
    
    nullStrJS <- I('function() { this.setValue(""); }')
    sidebar <- sidebarPanel(
        selectizeInput(ns("cellLine"), "Cell line", choices=NULL, width="100%",
                       options=list(placeholder='Select a cell line',
                                    onInitialize=nullStrJS)),
        selectizeInput(ns("gene"), "Gene", choices=NULL, width="100%",
                       options=list(placeholder='Select a gene',
                                    onInitialize=nullStrJS)),
        actionButton(ns("load"), "Load data", class="btn-primary"))
    main <- mainPanel(DTOutput(ns("table")))
    ui <- tabPanel(title, sidebarLayout(sidebar, main))
    return(ui)
}

#' @importFrom shiny moduleServer stopApp withProgress incProgress
#' @importFrom DT renderDT
.diffExprENCODEloaderServer <- function(id, metadata, cellLine=NULL, gene=NULL,
                                        path=".", globalUI=FALSE) {
    server <- function(input, output, session) {
        output$table <- renderDT({
            hiddenCols <- c(
                "File format", "File type", "File format type", "Donors",
                "Biosample treatments", "Biosample treatments amount",
                "Biosample treatments duration",
                "Biosample genetic modifications targets",
                "Biosample genetic modifications gene targets",
                "Biosample genetic modifications site coordinates",
                "Biosample genetic modifications zygosity",
                "Library lysis method", "Library crosslinking method",
                "Project", "RBNS protein concentration",
                "Read length", "Mapped read length", "Run type",
                "Paired end", "Paired with", "Index of", "md5sum", "dbxrefs",
                "File download URL", "File analysis status",
                "Platform", "Controlled by", "s3_uri",
                "Audit ERROR", "Audit WARNING")
            hiddenCols <- match(hiddenCols, colnames(metadata))
            columnDefs <- list(list(visible=FALSE, targets=hiddenCols - 1))
            .prepareDT(metadata, columnDefs=columnDefs)
        })
        proxy <- dataTableProxy("table")
        
        observe({
            conditions <- .getENCODEconditions(metadata)
            updateSelectizeInput(
                session, "gene", choices=conditions$genes, selected=gene)
            updateSelectizeInput(
                session, "cellLine", choices=conditions$cellLines,
                selected=cellLine)
        })
        
        observe({
            cellLine <- input$cellLine
            if (is.null(cellLine) || cellLine == "") cellLine <- NULL
            
            gene <- input$gene
            if (is.null(gene) || gene == "") gene <- NULL
            data <- filterENCONDEmetadata(metadata, cellLine, gene)
            observe(replaceData(proxy, data, rownames=FALSE))
        })
        
        loadData <- eventReactive(input$load, {
            withProgress(message="Preparing differential expression", {
                steps <- 3
                
                incProgress(1/steps, detail="Downloading ENCODE data")
                message("Downloading data...")
                ENCODEmetadata <- filterENCONDEmetadata(
                    metadata, input$cellLine, input$gene)
                
                if (nrow(ENCODEmetadata) == 0) {
                    showNotification("No samples match the selected criteria",
                                     type="error")
                    return(NULL)
                }
                if (is.function(path)) path <- path()
                ENCODEsamples <- loadENCODEsamples(ENCODEmetadata, path=path)
                ENCODEsamples <- ENCODEsamples[[1]]
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
        
        if (!globalUI) observeEvent(input$load, stopApp(loadData()))
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
                              globalUI=FALSE) {
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
    if (globalUI) {
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
        if (globalUI) textInput(ns("name"), "Dataset name", value="cmapData"),
        if (!globalUI)
            helpText("By default, data will be downloaded if not found."),
        if (!globalUI) actionButton(ns("cancel"), "Cancel"),
        actionButton(ns("load"), "Load data", class="btn-primary"),
        convertTaskState2HTML("Loaded", toStr=FALSE, id=ns("loaded"),
                              style="margin-left: 10px; opacity: 0;"))
    mainPanel <- mainPanel( withSpinner(DTOutput(ns("table")), type=6) )
    
    ui <- tabPanel(title, sidebarLayout(sidebar, mainPanel))
    ui[[3]] <- c(list(plots), ui[[3]])
    return(ui)
}

#' @importFrom shiny isolate updateSelectizeInput moduleServer stopApp
#' observeEvent eventReactive updateTextInput reactiveVal
#' @importFrom DT dataTableProxy replaceData
.cmapDataLoaderServer <- function(id, metadata="cmapMetadata.txt",
                                  zscores="cmapZscores.gctx",
                                  geneInfo="cmapGeneInfo.txt",
                                  compoundInfo="cmapCompoundInfo.txt",
                                  cellLine=NULL, timepoint=NULL, dosage=NULL,
                                  globalUI=FALSE, tab=NULL) {
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
                    load <- !is.null(tab) && tab == "CMap perturbations"
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
            session$sendCustomMessage("brieflyShowElem", session$ns("loaded"))
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
        
        if (!globalUI) observeEvent(input$load, stopApp(loadData()))
        observeEvent(input$cancel, stopApp(stop("User cancel", call.=FALSE)))
        return(loadData)
    }
    
    moduleServer(id, server)
}

# Result plotting and viewing --------------------------------------------------

#' @importFrom shiny NS sidebarPanel mainPanel tabPanel sidebarLayout
#' selectizeInput fluidRow column
.metadataViewerUI <- function(id, title="Data", icon=NULL) {
    ns <- NS(id)
    sidebar <- sidebarPanel(
        selectizeInput(ns("object"), "Dataset", choices=NULL, width="100%"),
        # DTOutput(ns("selection")),
        selectizeInput(ns("attr"), "Table", choices=NULL, width="100%"))
    main <- mainPanel(
        .alert("No dataset loaded/selected",
               condition="input.object == ''", ns=ns),
        DTOutput(ns("table")))
    ui <- tabPanel(title, icon=icon, sidebarLayout(sidebar, main))
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
            
            observe( .updateDatasetChoices(session, "object", x()) )
            
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
.dataPlotterUI <- function(id, title="Plot Data") {
    ns <- NS(id)
    sidebar <- sidebarPanel(
        selectizeInput(ns("object"), "Dataset", choices=NULL),
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

.comparisonMethods <- function() {
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
                .updateDatasetChoices(session, "object", x(),
                                      "referenceComparison")
            })
            
            # Update element and methods choices depending on selected object
            observeEvent(input$object, {
                obj <- getSelectedObject()
                if (is.null(obj)) return(NULL)
                
                updateSelectizeInput(session, "element", choices=obj[[1]],
                                     selected=list(), server=TRUE)
                # Update methods depending on selected object
                methods <- .comparisonMethods()
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
.datasetComparisonUI <- function(id, title="Dataset Comparison") {
    ns <- NS(id)
    sidebar <- sidebarPanel(
        selectizeInput(ns("data1"), "Dataset 1", choices=NULL),
        selectizeInput(ns("col1"), "Column to plot in X axis", choices=NULL),
        hr(),
        selectizeInput(ns("data2"), "Dataset 2", choices=NULL),
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
                .updateDatasetChoices(session, "data1", x(),
                                      "referenceComparison")
                .updateDatasetChoices(session, "data2", x(),
                                      "referenceComparison")
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
                        xlab(paste(input$data1, input$col1, sep="_")) +
                        ylab(paste(input$data2, input$col2, sep="_")) +
                        theme_bw()
                    
                    # Avoid showing density instead of calculating an error
                    qd1 <- quantile(dataset1[[col1]], c(0.25, 0.75))
                    qd2 <- quantile(dataset2[[col2]], c(0.25, 0.75))
                    if (diff(qd1) != 0 && diff(qd2) != 0) {
                        plot <- plot + geom_density_2d(alpha=0.3)
                    }
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
    id, title="Targeting Drugs vs Similar Perturbations") {
    
    ns <- NS(id)
    sidebar <- sidebarPanel(
        selectizeInput(ns("data1"), "Predicted targeting drugs", NULL),
        selectizeInput(ns("data2"), "Similar CMap perturbations", NULL),
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
                .updateDatasetChoices(session, "data1", x(), "targetingDrugs")
                .updateDatasetChoices(session, "data2", x(),
                                      "similarPerturbations")
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

getTaskState <- function(dataset) {
    if ("state" %in% names(dataset)) {
        state <- capitalize(dataset[["state"]])
    } else {
        state <- "Loaded"
    }
    return(state)
}

convertTaskState2HTML <- function(state, toStr=TRUE, ..., label=FALSE) {
    state <- capitalize(tolower(state))
    if (state %in% c("Failure", "Revoked")) {
        colour <- "red"
        icon   <- icon("times-circle")
        state  <- "Error"
        class  <- "label label-danger"
    } else if (state %in% c("Received", "Pending", "Retry")) {
        colour <- "grey"
        icon   <- icon("pause-circle")
        state  <- "Waiting"
        class  <- "label label-default"
    } else if (state %in% c("Started")) {
        colour <- "orange"
        icon   <- icon("circle-notch", "fa-spin")
        state  <- "Running"
        class  <- "label label-warning"
    } else if (state %in% c("Success", "Loaded")) {
        colour <- "green"
        icon   <- icon("check-circle")
        state  <- "Loaded"
        class  <- "label label-success"
    } else {
        return(state)
    }
    
    if (!label) {
        colour <- sprintf("color: %s;", colour)
        class  <- NULL
    } else {
        colour <- NULL
    }
    html   <- tags$span(style=colour, icon, state, ..., class=class)
    if (toStr) html <- as.character(html)
    return(html)
}

# Run rank similar perturbations in Celery
celery_rankAgainstRef <- function(..., mode, token) {
    # Prepare filenames for input and output
    rand       <- .genRandomString()
    inputFile  <- file.path(token, sprintf("input_%s.Rda",  rand))
    outputFile <- file.path(token, sprintf("output_%s.rds", rand))
    
    # Save variables in Rda file
    save(..., file=inputFile, envir=parent.frame())
    
    # Rank comparisons via Celery/Flower and save as RDS file
    if (mode == "similarPerturbations") {
        cmd <- "cTRAP::rankSimilarPerturbations(
                    selectedDiffExpr, selectedPerts, method,
                    c(upGenes, downGenes), cellLineMean, rankPerCellLine)"
    } else if (mode == "targetingDrugs") {
        cmd <- "cTRAP::predictTargetingDrugs(
                    selectedDiffExpr, selectedCorMatrix, method,
                    c(upGenes, downGenes))"
    }
    cmd <- list(sprintf("load('%s')", inputFile),
                paste("ranking <-", cmd),
                "attr(ranking, 'name') <- dataset",
                sprintf("saveRDS(ranking, '%s')", outputFile),
                sprintf("unlink('%s')", inputFile))
    cmd <- gsub("\n *", "", paste(cmd, collapse="; "))
    taskAsync <- floweRy::taskAsyncApply("tasks.R", cmd)
    
    # Prepare object
    ranking <- taskAsync
    ranking$state <- capitalize(ranking$state)
    ranking[["outputFile"]] <- outputFile
    class(ranking) <- c(paste0("expected", capitalize(mode)), "expected",
                        class(ranking))
    return(ranking)
}

#' @importFrom shiny NS sidebarPanel plotOutput selectizeInput mainPanel
#' tabPanel uiOutput column fluidRow numericInput checkboxInput conditionalPanel
#' @importFrom DT DTOutput
.rankSimilarPerturbationsUI <- function(
    id, title="Rank CMap perturbations by similarity") {
    
    ns <- NS(id)
    sidebar <- sidebarPanel(
        selectizeInput(ns("diffExpr"), "Differential gene expression",
                       choices=NULL),
        selectizeInput(ns("perts"), "CMap perturbations", choices=NULL),
        selectizeInput(ns("method"), "Method", multiple=TRUE,
                       .comparisonMethods(), .comparisonMethods(),
                       options=list(plugins=list("remove_button"))),
        conditionalPanel(
            "input.method.includes('gsea')",
            fluidRow(
                column(6, numericInput(ns("upGenes"), "Top genes", 150)),
                column(6, numericInput(ns("downGenes"), "Bottom genes", 150))),
            ns=ns),
        selectizeInput(ns("cellLineMean"),
                       "Calculate mean across cell lines",
                       c("For data with \u2265 2 cell lines"="auto",
                         "Always"=TRUE,
                         "Never"=FALSE)),
        conditionalPanel(
            "input.cellLineMean != 'FALSE'",
            selectizeInput(ns("rankPerCellLine"), "Rank results based on",
                           c("Mean scores only"=FALSE,
                             "Mean + individual cell lines' scores"=TRUE)),
            ns=ns),
        textInput(ns("name"), "Dataset name", "Ranked CMap perturbations"),
        uiOutput(ns("msg")),
        actionButton(ns("analyse"), "Rank by similarity", class="btn-primary"))
    
    mainPanel <- mainPanel(DTOutput(ns("table")))
    ui <- tabPanel(title, sidebarLayout(sidebar, mainPanel))
    return(ui)
}

#' @importFrom shiny renderPlot observeEvent observe isolate renderUI
#' @importFrom DT renderDT
#' @importFrom data.table rbindlist
.rankSimilarPerturbationsServer <- function(id, x, globalUI=FALSE,
                                            flower=FALSE, token=NULL) {
    server <- function(input, output, session) {
        observe({
            .updateDatasetChoices(session, "diffExpr", x(), "diffExpr")
            .updateDatasetChoices(session, "perts", x(), "perturbationChanges")
        })
        
        rankData <- eventReactive(input$analyse, {
            diffExprDataset <- req(input$diffExpr)
            pertsDataset    <- req(input$perts)
            method          <- input$method
            upGenes         <- input$upGenes
            downGenes       <- input$downGenes
            cellLineMean    <- input$cellLineMean
            rankPerCellLine <- input$rankPerCellLine
            dataset         <- input$name
            
            if (cellLineMean == "TRUE") {
                cellLineMean <- TRUE
                cellLineMeanTxt <- "Always"
            } else if (cellLineMean == "FALSE") {
                cellLineMean <- FALSE
                cellLineMeanTxt <- "Never"
            } else if (cellLineMean == "auto") {
                cellLineMeanTxt <- "For data with \u2265 2 cell lines" 
            }
            
            if (rankPerCellLine == "TRUE") {
                rankPerCellLine <- TRUE
            } else if (rankPerCellLine == "FALSE") {
                rankPerCellLine <- FALSE
            }
            
            selectedDiffExpr <- x()[[req(diffExprDataset)]]
            selectedPerts    <- x()[[req(pertsDataset)]]
            if (!flower) {
                withProgress(message="Ranking against CMap perturbations", {
                    ranking <- rankSimilarPerturbations(
                        selectedDiffExpr, selectedPerts, method,
                        c(upGenes, downGenes), cellLineMean, rankPerCellLine)
                    incProgress(1)
                })
            } else {
                ranking <- celery_rankAgainstRef(
                    selectedDiffExpr, selectedPerts, method, upGenes, downGenes,
                    cellLineMean, rankPerCellLine, dataset,
                    token=isolate(token()), mode="similarPerturbations")
            }
            attr(ranking, "name") <- dataset
            
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
                "Calculate mean across cell lines"=cellLineMeanTxt,
                "Rank results based on"=rankPerCellLine)
            return(ranking)
        })
        
        if (!globalUI) observeEvent(input$load, stopApp(rankData()))
        
        output$table <- renderDT({
            .prepareReferenceComparisonDT(x(), "similarPerturbations")
        }, rownames=FALSE, escape=FALSE, selection="none")
        
        return(rankData)
    }
    moduleServer(id, server)
}

#' @importFrom shiny NS sidebarPanel plotOutput selectizeInput mainPanel
#' tabPanel uiOutput column fluidRow numericInput checkboxInput conditionalPanel
#' @importFrom DT DTOutput
.predictTargetingDrugsUI <- function(id, title="Predict targeting drugs") {
    ns <- NS(id)
    sidebar <- sidebarPanel(
        selectizeInput(ns("diffExpr"), choices=NULL,
                       "Differential gene expression"),
        selectizeInput(ns("corMatrix"),
                       "Gene expression and drug sensitivity association",
                       choices=listExpressionDrugSensitivityAssociation()),
        selectizeInput(ns("method"), "Method", multiple=TRUE,
                       .comparisonMethods(), .comparisonMethods(),
                       options=list(plugins=list("remove_button"))),
        conditionalPanel(
            "input.method.includes('gsea')",
            fluidRow(
                column(6, numericInput(ns("upGenes"), "Top genes", 150)),
                column(6, numericInput(ns("downGenes"), "Bottom genes", 150))),
            ns=ns),
        textInput(ns("name"), "Dataset name", "Targeting drugs"),
        uiOutput(ns("msg")),
        actionButton(ns("analyse"), "Predict targeting drugs",
                     class="btn-primary"))
    
    mainPanel <- mainPanel(DTOutput(ns("table")))
    ui <- tabPanel(title, sidebarLayout(sidebar, mainPanel))
    return(ui)
}

#' @importFrom shiny renderPlot observeEvent observe isolate renderUI
#' @importFrom DT renderDT
#' @importFrom data.table rbindlist
.predictTargetingDrugsServer <- function(id, x, path=".", globalUI=FALSE,
                                         flower=FALSE, token=NULL) {
    server <- function(input, output, session) {
        observe( .updateDatasetChoices(session, "diffExpr", x(), "diffExpr") )
        
        observe({
            dataset <- "Targeting drugs"
            
            corMatrix <- input$corMatrix
            if (!is.null(corMatrix) || corMatrix != "") {
                dataset <- paste(corMatrix, tolower(dataset))
            }
            updateTextInput(session, "name", value=dataset)
        })
        
        rankData <- eventReactive(input$analyse, {
            diffExprDataset <- req(input$diffExpr)
            corMatrix       <- req(input$corMatrix)
            method          <- input$method
            upGenes         <- input$upGenes
            downGenes       <- input$downGenes
            dataset         <- input$name
            
            selectedDiffExpr  <- x()[[diffExprDataset]]
            selectedCorMatrix <- loadExpressionDrugSensitivityAssociation(
                corMatrix, path=path)
            if (!flower) {
                withProgress(message="Predict targeting drugs", {
                    ranking <- predictTargetingDrugs(
                        selectedDiffExpr, selectedCorMatrix, method,
                        c(upGenes, downGenes))
                    incProgress(1)
                })
            } else {
                ranking <- celery_rankAgainstRef(
                    selectedDiffExpr, selectedCorMatrix, method,
                    upGenes, downGenes, dataset, token=isolate(token()),
                    mode="targetingDrugs")
            }
            attr(ranking, "name") <- dataset
            
            # Prepare form input
            attr(ranking, "formInput") <- list(
                "Differential expression dataset"=diffExprDataset,
                "Gene expression and drug sensitivity association"=corMatrix,
                "Methods"=paste(method, collapse=", "),
                "Top genes"=upGenes,
                "Bottom genes"=downGenes)
            return(ranking)
        })
        
        if (!globalUI) observeEvent(input$load, stopApp(rankData()))
        
        output$table <- renderDT({
            .prepareReferenceComparisonDT(x(), "targetingDrugs")
        }, rownames=FALSE, escape=FALSE, selection="none")
        
        return(rankData)
    }
    moduleServer(id, server)
}

#' @importFrom shiny NS sidebarPanel plotOutput selectizeInput mainPanel
#' tabPanel helpText textOutput
#' @importFrom shinycssloaders withSpinner
#' @importFrom DT DTOutput
.drugSetEnrichmentAnalyserUI <- function(id, title="Drug Set Enrichment") {
    ns <- NS(id)
    sidebar <- sidebarPanel(
        selectizeInput(ns("object"), "Dataset", choices=NULL),
        selectizeInput(ns("sort"), "Ranked by", choices=NULL), hr(),
        selectizeInput(ns("drugSet"), "Drug set", choices=NULL,
                       options=list(placeholder="Select a drug set")), hr(),
        tags$h4("Match compounds"),
        helpText("Select key columns to match compounds between datasets"),
        withSpinner(uiOutput(ns("matchCompounds")), type=8,
                    proxy.height="150px", hide.ui=TRUE))
    sidebar[[3]][[2]] <- tagList(
        selectizeInput(ns("element"), "Row ID to plot", choices=NULL,
                       width="100%"),
        plotOutput(ns("plot")))
    
    mainPanel <- mainPanel(DTOutput(ns("table")))
    ui <- tabPanel(title, sidebarLayout(sidebar, mainPanel))
    return(ui)
}

#' @importFrom shiny renderPlot observeEvent observe isolate renderText
#' @importFrom DT renderDT
.drugSetEnrichmentAnalyserServer <- function(id, x, path=NULL) {
    x <- .convertToFunction(x)
    moduleServer(
        id,
        function(input, output, session) {
            getExtraDrugSets <- function() {
                extraDrugSets <- c("NCI60 2D", "NCI60 3D", "CMap 2D", "CMap 3D")
                extraDrugSets <- setNames(
                    extraDrugSets,
                    paste(extraDrugSets, "molecular descriptors"))
                return(extraDrugSets)
            }
            
            getSelectedObject <- reactive(x()[[input$object]])
            getSelectedSet <- reactive({
                sets <- .filterDatasetsByClass(x(), "drugSets")
                drugSet <- input$drugSet
                if (drugSet %in% names(sets)) {
                    res <- sets[[req(input$drugSet)]]
                } else if (drugSet %in% getExtraDrugSets()) {
                    drugSet <- strsplit(drugSet, " ")[[1]]
                    withProgress(message="Loading drug descriptors", {
                        res <- loadDrugSet(drugSet[[1]], drugSet[[2]],
                                           path=path)
                        incProgress(1)
                        return(res)
                    })
                } else {
                    res <- NULL
                }
                return(res)
            })
            
            # Update available datasets
            observe({
                .updateDatasetChoices(session, "object", x(),
                                      "referenceComparison")
            })
            
            # Update available drug sets
            observe({
                drugSets <- .filterDatasetsByClass(x(), "drugSets")
                drugSets <- names(drugSets)
                
                # Show all available drug sets from cTRAP
                if (is.null(drugSets)) {
                    choices <- list(
                        "Other available drug sets"=getExtraDrugSets())
                    selected <- list()
                } else {
                    choices <- list(
                        "Loaded drug sets"=list(drugSets),
                        "Other available drug sets"=getExtraDrugSets())
                    selected <- drugSets[[1]]
                }
                updateSelectizeInput(session, "drugSet", choices=choices,
                                     selected=selected)
            })
            
            getDSEAresult <- reactive({
                obj      <- req(getSelectedObject())
                sets     <- req(getSelectedSet())
                sort     <- input$sort
                statsKey <- input$statsKey
                setsKey  <- input$setsKey
                
                isValid <- function(e) !is.null(e) && e != ""
                if (is.null(obj) || !isValid(sort)) return(NULL)
                if (!isValid(statsKey) || !isValid(setsKey)) return(NULL)
                
                withProgress(message="Analysing drug set enrichment", {
                    res <- analyseDrugSetEnrichment(
                        sets, obj, col=sort,
                        keyColSets=setsKey, keyColStats=statsKey)
                    incProgress(1)
                    return(res)
                })
            })
            
            observeEvent(input$object, {
                obj <- getSelectedObject()
                numericCols <- names(obj)[vapply(obj, is.numeric, logical(1))]
                updateSelectizeInput(session, "sort", choices=numericCols)
            })
            
            # Update available keys to select for datasets
            getCompoundMatchKeys <- function(statsKey=NULL, setsKey=NULL) {
                obj <- req(getSelectedObject())
                sets <- req(getSelectedSet())
                
                statsInfo <- prepareStatsCompoundInfo(obj)$statsInfo
                setsInfo  <- prepareSetsCompoundInfo(sets)$setsCompoundInfo
                
                probableKey <- findIntersectingCompounds(
                    statsInfo, setsInfo, keys1=statsKey, keys2=setsKey)
                setsKey     <- probableKey$key2
                statsKey    <- probableKey$key1
                
                keyList      <- getCompoundIntersectingKeyList()
                statsOptions <- intersect(names(statsInfo), keyList)
                setsOptions  <- intersect(names(setsInfo), keyList)
                return(list(probableKey=probableKey,
                            setsKey=setsKey, statsKey=statsKey,
                            setsOptions=setsOptions, statsOptions=statsOptions))
            }
            
            # Update interface for selecting columns for compound matching
            output$matchCompounds <- renderUI({
                ns <- session$ns
                ui <- tagList(
                    fluidRow(
                        column(6, selectizeInput(ns("statsKey"), choices=NULL,
                                                 "Dataset key")),
                        column(6, selectizeInput(ns("setsKey"), choices=NULL,
                                                 "Drug set key"))),
                    tags$span(class="help-block", style="margin-top: -10px;",
                              textOutput(ns("msg"))),
                    actionButton(ns("analyse"), "Visualise",
                                 class="btn-primary"))
                
                res <- suppressMessages( getCompoundMatchKeys() )
                updateSelectizeInput(session, "statsKey", selected=res$statsKey,
                                     choices=res$statsOptions)
                updateSelectizeInput(session, "setsKey", selected=res$setsKey,
                                     choices=res$setsOptions)
                return(ui)
            })
            
            # Update number of intersecting compounds based on selected keys
            observe({
                res <- getCompoundMatchKeys(req(input$statsKey),
                                            req(input$setsKey))
                
                num <- length(res$probableKey[[3]])
                msg <- "cross-matches with selected columns"
                output$msg <- renderText(paste(num, msg))
            })
            
            observeEvent(input$analyse, {
                dsea <- getDSEAresult()
                updateSelectizeInput(session, "element", choices=dsea[[1]])
                output$table <- renderDT({
                    hiddenCols <- "leadingEdge"
                    hiddenCols <- match(hiddenCols, colnames(dsea))
                    columnDefs <- list(list(visible=FALSE,
                                            targets=hiddenCols - 1))
                    .prepareDT(dsea, columnDefs=columnDefs)
                })
                
                output$plot <- renderPlot({
                    obj  <- req(getSelectedObject())
                    sets <- req(getSelectedSet())
                    
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
                    suppressMessages(
                        plotDrugSetEnrichment(sets, obj, col=sort,
                                              selectedSets=element,
                                              keyColStats=statsKey,
                                              keyColSets=setsKey)[[1]])
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
#' @inheritParams loadENCODEsamples
#'
#' @return Differential expression data
#' @family visual interface functions
#' @export
launchDiffExprLoader <- function(cellLine=NULL, gene=NULL,
                                 file="ENCODEmetadata.rds", path=".") {
    metadata <- downloadENCODEknockdownMetadata(file=file)
    id       <- "diffExpr"
    ui       <- .prepareNavPage(.diffExprENCODEloaderUI(id))
    server   <- function(input, output, session) {
        .diffExprENCODEloaderServer(id, metadata, cellLine, gene, path=path)
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
        .dataPlotterUI(dataId),
        .targetingDrugsVSsimilarPerturbationsPlotterUI(comparePlotId),
        .datasetComparisonUI(compareId),
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
        .drugSetEnrichmentAnalyserUI(dseaId),
        .dataPlotterUI(dataId),
        .metadataViewerUI(metadataId))
    uiList <- Filter(length, uiList)
    ui     <- do.call(.prepareNavPage, uiList)
    
    server <- function(input, output, session) {
        .drugSetEnrichmentAnalyserServer(dseaId,
                                         c(elems, list("Custom drug set"=sets)))
        .dataPlotterServer(dataId, elems)
        .metadataViewerServer(metadataId, elems)
    }
    app <- runApp(shinyApp(ui, server))
    return(app)
}
