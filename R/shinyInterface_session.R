.setFileSizeLimit <- function(limitMiB) {
    options(shiny.maxRequestSize = limitMiB * 1024^2)
    message("cTRAP: file upload size limit set to ", limitMiB, " MiB")
}

.createToken <- function(len=10) {
    pool  <- list(LETTERS, letters, 0:9)
    size  <- sapply(pool, length)
    prob  <- rep(1/size, size)
    pool  <- unlist(pool)
    
    rand  <- sample(pool, len, replace=TRUE, prob=prob)
    token <- paste(rand, collapse="")
    return(token)
}

.prepareSessionModal <- function(title=NULL, createSession=TRUE,
                                 footer=modalButton("Dismiss"), ...) {
    newSessionUI <- tagList(
        h2("New session", style="margin-top: 0px;"),
        actionButton("createSession", "Create new session",
                     width="100%", icon=icon("plus"), class="btn-info"),
        tags$hr())
    
    loadTokenUI <- tagList(
        textInput("token", "Insert token of a previous session:"),
        actionButton("loadToken", "Load session with token",
                     width="100%", icon=icon("history"), class="btn-info"))
    loadDataUI <- tagList(
        fileInput("sessionFile", width="100%", multiple=TRUE, accept=".rds",
                  "Insert RDS file of a previous session:"),
        actionButton("loadData", "Load session with files", width="100%",
                     icon=icon("history"), class="btn-info"))
    
    modalDialog(
        title=title,
        if (createSession) newSessionUI,
        h2("Load session", style="margin-top: 0px;"),
        tabsetPanel(
            type="pills",
            tabPanel("Session token", loadTokenUI),
            tabPanel("Session data", loadDataUI)),
        size="s", footer=footer, ...)
}

.modifySessionUI <- function(ui, expire) {
    # Modify session
    session <- ui[[3]][[1]][[3]][[1]][[3]][[1]][[3]][[2]][[3]][[1]][[4]]
    ui[[3]][[1]][[3]][[1]][[3]][[1]][[3]][[2]][[3]][[1]][[4]] <- NULL
    
    # Add session buttons
    copyTokenButton <- actionLink(
        "copyToken", onclick="copyToken()", tagList(
            "Copy session token to clipboard",
            if (!is.null(expire)) {
                helpText(style="margin: 0px; padding: 3px 0px;",
                         paste("Sessions expire in", expire, "days"))  
            }))
    session[[3]][[2]][[3]] <- tagList(
        tags$li(role="presentation", copyTokenButton),
        tags$li(role="presentation", class="divider"),
        tags$li(role="presentation",
                downloadLink("downloadSession", "Download session data")),
        tags$li(role="presentation",
                actionLink("loadSessionModal", "Load another session")))
    session[[3]][[2]] <- tagAppendAttributes(session[[3]][[2]],
                                             class="pull-right")
    
    # Place session in the right side of the navigation bar
    ui[[3]][[1]][[3]][[1]][[3]][[1]][[3]][[3]] <- tags$ul(
        class="nav navbar-nav pull-right", session)
    return(ui)
}

globalUI <- function(elems, idList, expire) {
    elemClasses       <- sapply(lapply(elems, class), "[[", 1)
    hasSimilarPerts   <- "similarPerturbations" %in% elemClasses
    hasTargetingDrugs <- "targetingDrugs" %in% elemClasses
    showTwoKindPlot   <- hasSimilarPerts && hasTargetingDrugs
    
    uiList <- tagList(
        id="tab",
        navbarMenu("Load", icon=icon("table"),
                   "Differential gene data expression",
                   .diffExprENCODEloaderUI(idList$encode),
                   "----",
                   "CMap data",
                   .cmapDataLoaderUI(idList$cmap, shinyproxy=TRUE)),
        navbarMenu("Analyse", icon=icon("cogs"),
                   .rankSimilarPerturbationsUI(idList$rankPerts, elems, elems),
                   .drugSetEnrichmentAnalyserUI(idList$drugSet, elems, elems)),
        navbarMenu("Visualise", icon=icon("chart-bar"),
                   .dataPlotterUI(idList$data, elems),
                   # .targetingDrugsVSsimilarPerturbationsPlotterUI(
                   #     idList$comparePlot, elems, elemClasses),
                   .datasetComparisonUI(idList$compare, elems),
                   .metadataViewerUI(idList$metadata, elems)),
        navbarMenu(span("Session",
                        span(class="badge", textOutput("token", inline=TRUE))),
                   icon=icon("compass")))
    ui <- do.call(.prepareNavPage, uiList)
    ui <- .modifySessionUI(ui, expire=expire)
    
    # Add JS and CSS in header
    header <- tags$head(
        includeScript(system.file("shiny", "www", "cTRAP.js", package="cTRAP")),
        includeCSS(system.file("shiny", "www", "cTRAP.css", package="cTRAP")))
    ui <- tagList(header, ui)
    return(ui)
}

#' Complete visual interface with support for sessions
#'
#' @param ... Objects
#' @param expire Character: days until a session expires (message purposes only)
#' @param fileSizeLimitMB Numeric: file size limit in MiB
#'
#' @importFrom shiny tagList showModal modalButton modalDialog removeModal
#' reactiveValues tagAppendAttributes showNotification
#'
#' @return Launches result viewer and plotter (returns \code{NULL})
#' @family visual interface functions
#' @export
cTRAP <- function(..., expire=14, fileSizeLimitMiB=50) {
    .setFileSizeLimit(fileSizeLimitMiB)
    elems <- .prepareEllipsis(...)
    
    idList             <- list()
    idList$encode      <- "encodeDataLoader"
    idList$cmap        <- "cmapDataLoader"
    idList$compare     <- "datasetComparison"
    idList$comparePlot <- "comparePlotter"
    idList$data        <- "dataPlotter"
    idList$metadata    <- "metadataViewer"
    idList$rankPerts   <- "rankPerts"
    idList$drugSet     <- "drugSetAnalyser"
    ui <- globalUI(elems, idList, expire)
    
    sharedData <- reactiveValues()
    sharedData$elems <- elems
    
    server <- function(input, output, session) {
        .diffExprENCODEloaderServer(idList$encode)
        .cmapDataLoaderServer(idList$cmap)
        
        observe({
            elems <- sharedData$elems
            
            # analyse
            .rankSimilarPerturbationsServer(idList$rankPerts, elems, elems)
            .drugSetEnrichmentAnalyserServer(idList$drugSet, elems, elems)
            
            # visualise
            .dataPlotterServer(idList$data, elems)
            .targetingDrugsVSsimilarPerturbationsPlotterServer(
                idList$comparePlot, elems)
            .datasetComparisonServer(idList$compare, elems)
            .metadataViewerServer(idList$metadata, elems)
        })
        
        # When no token is found (i.e. new session), show welcome screen
        observe({
            if (!is.null(sharedData$token)) return(NULL)
            showModal(.prepareSessionModal("Welcome to cTRAP!", footer=NULL))
        })
        
        # Create new session
        observeEvent(input$createSession, {
            removeModal()
            sharedData$token <- .createToken()
        })
        
        # Update token display
        output$token <- renderText({
            token <- sharedData$token
            if (is.null(token)) token <- "?"
            return(token)
        })
        
        # Load session based on a token
        observeEvent(input$loadToken, {
            removeModal()
            sharedData$elems <- list(
                cmapKD=cmapKD, cmapCompounds=cmapCompounds, cmapPerts=cmapPerts,
                compareCompounds=compareCompounds, compareKD=compareKD)
        })
        
        # Load session based on a RDS file
        observeEvent(input$loadData, {
            file <- input$sessionFile
            
            if (is.null(file)) {
                showNotification("File input cannot be empty", type="error")
            }
            req(file)
            
            elems <- tryCatch({ readRDS(file$datapath) }, error=return)
            if (is(data, "error")) {
                showNotification("Error loading data", type="error")
            } else {
                sharedData$elems <- elems
                removeModal()
            }
        })
        
        observeEvent(input$loadSessionModal, {
            modal <- .prepareSessionModal(createSession=FALSE, easyClose=TRUE)
            showModal(modal)
        })
        
        # Notify when copying token
        observeEvent(input$copyToken, {
            msg <- tagList("Token", span(class="badge", sharedData$token),
                           "copied to your clipboard!")
            showNotification(msg, duration=3, closeButton=FALSE, type="message")
        })
        
        # Download objects in current session in a single RDS file
        output$downloadSession <- downloadHandler(
            filename=function() paste0("cTRAP-", sharedData$token, ".rds"),
            content=function(file) saveRDS(sharedData$elems, file))
    }
    app <- runApp(shinyApp(ui, server))
    return(app)
}
