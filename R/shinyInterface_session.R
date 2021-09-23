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
    if (dir.exists(token)) token <- .createToken(len)
    return(token)
}

.addToList <- function(x, data, name=NULL) {
    if (is.null(name)) name <- attr(data, "name")
    if (is.null(name) || name == "") name <- "Dataset"
    
    uniqName <- make.unique(c(names(x), name))
    name <- uniqName[[length(uniqName)]]
    x[[name]] <- data
    return(x)
}

.saveSession <- function(data, token) {
    if (is.null(token) || is.null(data)) return(NULL)
    if (!dir.exists(token)) dir.create(token)
    sessionRDS <- file.path(token, "session.rds")
    saveRDS(data, sessionRDS)
    message("     Session saved to ", sessionRDS)
}

#' @importFrom shiny textInput
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
                  "Upload RDS file of a previous session:"),
        actionButton("loadData", "Load session from RDS file", width="100%",
                     icon=icon("history"), class="btn-info"))
    pills <- tabsetPanel(
        type="pills",
        tabPanel("Session token", loadTokenUI),
        tabPanel("Session data", loadDataUI))
    pills[[3]][[1]] <- tagAppendAttributes(pills[[3]][[1]],
                                           class="nav-justified")
    modalDialog(
        title=title, size="s", footer=footer,
        if (createSession) newSessionUI,
        h2("Load session", style="margin-top: 0px;"), pills, ...)
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

#' @importFrom shiny tagAppendChildren
.addLoadingStatus <- function(ui) {
    loading <- conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                icon("circle-notch", "fa-spin"))
    loading$name <- "a"
    loading <- tags$li(loading)
    ui[[3]][[1]][[3]][[1]][[3]][[1]][[3]][[3]][[3]][[2]] <- 
        ui[[3]][[1]][[3]][[1]][[3]][[1]][[3]][[3]][[3]][[1]]
    ui[[3]][[1]][[3]][[1]][[3]][[1]][[3]][[3]][[3]][[1]] <- loading
    return(ui)
}

#' @importFrom shiny navbarMenu icon span textOutput tags includeScript
#' includeCSS
globalUI <- function(elems, idList, expire) {
    elemClasses       <- sapply(lapply(elems, class), "[[", 1)
    hasSimilarPerts   <- "similarPerturbations" %in% elemClasses
    hasTargetingDrugs <- "targetingDrugs" %in% elemClasses
    showTwoKindPlot   <- hasSimilarPerts && hasTargetingDrugs
    
    uiList <- tagList(
        id="tab",
        navbarMenu("Load", icon=icon("table"),
                   "Differential gene expression data",
                   .diffExprLoadUI(idList$diffExpr),
                   .diffExprENCODEloaderUI(idList$encode),
                   "----",
                   .cmapDataLoaderUI(idList$cmap, shinyproxy=TRUE)),
        navbarMenu("Analyse", icon=icon("cogs"),
                   .rankSimilarPerturbationsUI(idList$rankPerts, elems, elems),
                   .drugSetEnrichmentAnalyserUI(idList$drugSet, elems, elems)),
        navbarMenu("Visualise", icon=icon("chart-bar"),
                   .dataPlotterUI(idList$data, elems),
                   # .targetingDrugsVSsimilarPerturbationsPlotterUI(
                   #     idList$comparePlot, elems, elemClasses),
                   .datasetComparisonUI(idList$compare, elems),
                   .metadataViewerUI(idList$metadata)),
        navbarMenu(span("Session", span(class="badge",
                                        textOutput("token", inline=TRUE))),
                   icon=icon("compass"), menuName="session"))
    ui <- do.call(.prepareNavPage, uiList)
    ui <- .modifySessionUI(ui, expire=expire)
    ui <- .addLoadingStatus(ui)
    
    # Add JS and CSS in header
    header <- tags$head(
        includeScript(system.file("shiny", "www", "cTRAP.js", package="cTRAP")),
        includeCSS(system.file("shiny", "www", "cTRAP.css", package="cTRAP")))
    ui <- tagList(header, ui)
    return(ui)
}

.sessionManagementServer <- function(input, output, session, sharedData) {
    # Show welcome screen when no token is set (e.g. new cTRAP sessions)
    observe({
        if (!is.null(sharedData$token)) return(NULL)
        showModal(.prepareSessionModal("Welcome to cTRAP!", footer=NULL))
    })
    
    # Create new session
    observeEvent(input$createSession, {
        sharedData$token <- .createToken()
        removeModal()
    })
    
    # Update token badge
    output$token <- renderText({
        token <- sharedData$token
        if (is.null(token)) token <- "?"
        return(token)
    })
    
    # Load session based on a token
    observeEvent(input$loadToken, {
        token <- isolate(input$token)
        if (dir.exists(token)) {
            sharedData$elems <- readRDS(file.path(token, "session.rds"))
            sharedData$token <- token
            removeModal()
        } else {
            msg <- tagList("Token", span(class="badge", token),
                           "does not exist")
            showNotification(type="error", msg)
        }
    })
    
    # Load session based on a RDS file
    observeEvent(input$loadData, {
        file <- input$sessionFile
        
        if (is.null(file)) {
            showNotification("File input cannot be empty", type="error")
        }
        req(file)
        
        data <- tryCatch(readRDS(file$datapath), error=return)
        if (is(data, "error")) {
            showNotification(paste("Error loading data:", data),
                             type="error")
        } else {
            sharedData$elems <- data
            sharedData$token <- token <- .createToken()
            removeModal()
            .saveSession(data, token)
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

.newDataNotification <- function(what, total, ..., auto=FALSE) {
    plural <- ifelse(total == 1, "", "s")
    totalTxt <- sprintf("Total: %s dataset%s", total, plural)
    message(sprintf("  -> %s (%s)",
                    paste(paste(what, collapse=" + "), "loaded"),
                    tolower(totalTxt)))
    
    len <- length(what)
    auto <- ifelse(auto, "automatically ", "")
    head <- "New %sloaded dataset:"
    if (len != 1) head <- paste(length(what), "new %sloaded datasets:")
    head <- sprintf(head, auto)
    
    what <- do.call(tags$ul, lapply(what, tags$li))
    showNotification(tagList(tags$b(head), what, totalTxt),
                     type="message", ...)
}

.loadDataFromRDSinDirServer <- function(input, output, session, sharedData) {
    checkNewRDSfiles <- function(path=".") {
        if (is.null(path)) return(NULL)
        res <- list.files(path, ".rds", ignore.case=TRUE)
        res <- res[res != "session.rds"]
        if (length(res) == 0) res <- NULL
        res <- file.path(path, res)
        return(res)
    }
    
    getNewRDSfiles <- reactivePoll(
        5000, session,
        checkFunc=function() checkNewRDSfiles(sharedData$token),
        valueFunc=function() checkNewRDSfiles(sharedData$token))
    
    observe({
        req(getNewRDSfiles())
        
        elems <- isolate(sharedData$elems)
        token <- isolate(sharedData$token)
        
        added <- 0
        for (i in getNewRDSfiles()) {
            message("Adding data from ", i, "...")
            obj <- try(readRDS(i), silent=TRUE)
            if (is(obj, "try-error")) {
                warning(obj)
                return(NULL)
            }
            elems <- .addToList(elems, obj)
            unlink(i)
            added <- added + 1
        }
        if (added == 0) return(NULL)
        sharedData$elems <- elems
        .newDataNotification(tail(names(elems), added), length(elems),
                             duration=NULL, auto=TRUE)
        .saveSession(elems, token)
    })
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
    idList$diffExpr    <- "diffExprLoader"
    idList$encode      <- "encodeDataLoader"
    idList$cmap        <- "cmapDataLoader"
    idList$compare     <- "datasetComparison"
    idList$comparePlot <- "comparePlotter"
    idList$data        <- "dataPlotter"
    idList$metadata    <- "metadataViewer"
    idList$rankPerts   <- "rankPerts"
    idList$drugSet     <- "drugSetAnalyser"
    ui <- globalUI(elems, idList, expire)
    
    server <- function(input, output, session) {
        sharedData       <- reactiveValues()
        sharedData$elems <- elems
        elems <- reactive(sharedData$elems)
        
        updateSharedData <- function(x) {
            observe({
                obj <- x()
                elems <- .addToList(isolate(sharedData$elems), obj)
                sharedData$elems <- elems
                
                .newDataNotification(tail(names(elems), 1), length(elems))
                .saveSession(elems, isolate(sharedData$token))
            })
        }
        
        # load data
        diffExpr <- .diffExprLoadServer(idList$diffExpr, elems)
        updateSharedData(diffExpr)
        
        encodeDiffExpr <- .diffExprENCODEloaderServer(
            idList$encode, shinyproxy=TRUE)
        updateSharedData(encodeDiffExpr)
        
        cmapData <- .cmapDataLoaderServer(
            idList$cmap, shinyproxy=TRUE, tab=reactive(session$input$tab))
        updateSharedData(cmapData)
        
        # analyse
        ranking <- .rankSimilarPerturbationsServer(idList$rankPerts, elems,
                                                   elems, shinyproxy=TRUE)
        updateSharedData(ranking)
        # .drugSetEnrichmentAnalyserServer(idList$drugSet, elems, elems)
        
        # visualise
        .dataPlotterServer(idList$data, elems)
        .targetingDrugsVSsimilarPerturbationsPlotterServer(
            idList$comparePlot, elems)
        .datasetComparisonServer(idList$compare, elems)
        .metadataViewerServer(idList$metadata, elems)
        
        .sessionManagementServer(input, output, session, sharedData)
        .loadDataFromRDSinDirServer(input, output, session, sharedData)
    }
    app <- runApp(shinyApp(ui, server))
    return(app)
}
