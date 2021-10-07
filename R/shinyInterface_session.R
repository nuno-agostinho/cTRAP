# Set size limit for user-uploaded files
.setFileSizeLimit <- function(limitMiB) {
    options(shiny.maxRequestSize = limitMiB * 1024^2)
    message("cTRAP: file upload size limit set to ", limitMiB, " MiB")
}

# Generate random string of given length
.genRandomString <- function(len=10) {
    pool <- list(LETTERS, letters, 0:9)
    size <- sapply(pool, length)
    prob <- rep(1/size, size)
    pool <- unlist(pool)
    
    rand <- sample(pool, len, replace=TRUE, prob=prob)
    str  <- paste(rand, collapse="")
    return(str)
}

# Create unique token
# Avoids creating a token that matches the name of a local folder
.createToken <- function(len=10, path=".") {
    repeat {
        token <- .genRandomString(len)
        # Token is available if no existing folder is named after the token
        isTokenAvailable <- !dir.exists(file.path(path, token))
        if (isTokenAvailable) break
    }
    return(token)
}

# Add elements to a named list (ensures unique names for each element)
.addToList <- function(x, data, name=NULL) {
    if (is.null(name)) name <- attr(data, "name")
    if (is.null(name) || name == "") name <- "Dataset"
    
    uniqName <- make.unique(c(names(x), name))
    name <- uniqName[[length(uniqName)]]
    x[[name]] <- data
    return(x)
}

# Save session data in token-named directory
.saveSession <- function(data, token) {
    if (is.null(token) || is.null(data)) return(NULL)
    if (!dir.exists(token)) dir.create(token)
    sessionRDS <- file.path(token, "session.rds")
    saveRDS(data, sessionRDS)
    message("     Session saved to ", sessionRDS)
}

#' @importFrom shiny textInput tags tabsetPanel fileInput
.prepareSessionModal <- function(title=NULL, createSession=TRUE,
                                 footer=modalButton("Dismiss"), ...) {
    newSessionUI <- tagList(
        tags$h2("New session", style="margin-top: 0px;"),
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
    pills <- tagAppendAttributes(pills, class="nav-justified",
                                 .cssSelector=".nav")
    modalDialog(
        title=title, size="s", footer=footer,
        if (createSession) newSessionUI,
        tags$h2("Load session", style="margin-top: 0px;"), pills, ...)
}

#' Find an item in list of lists and return its coordinates
#' @keywords internal
.traceInList <- function(ll, item) {
    if (is.list(ll)) {
        for (elem in seq(ll)) {
            res <- .traceInList(ll[[elem]], item)
            if (!is.null(res)) return(c(elem, res))
        }
    } else if (is.character(ll)) {
        if (any(grepl(item, ll, fixed=TRUE))) return(numeric(0))
    }
}

# Add context menu to session button in navigation bar
#' @importFrom purrr pluck pluck<-
#' @importFrom rlang !!!
#' @importFrom shiny actionLink downloadLink
.modifySessionUI <- function(ui, expire) {
    # Modify session
    pos     <- .traceInList(ui, "session")
    pos     <- head(pos, -4)
    session <- pluck(ui, !!!pos)
    pluck(ui, !!!pos) <- NULL
    
    # Add session buttons
    expireTxt <- NULL
    if (!is.null(expire)) {
        expireTxt <- helpText(style="margin: 0px; padding: 3px 0px;",
                              paste("Sessions expire in", expire, "days"))
    }
    
    copyTokenButton <- actionLink("copyToken", onclick="copyToken()", tagList(
        "Copy session token to clipboard", expireTxt))
    pluck(session, 3, 2, 3) <- tagList(
        tags$li(role="presentation", copyTokenButton),
        tags$li(role="presentation", class="divider"),
        tags$li(role="presentation",
                downloadLink("downloadSession", "Download session data")),
        tags$li(role="presentation",
                actionLink("loadSessionModal", "Load another session")))
    pluck(session, 3, 2) <- tagAppendAttributes(
        pluck(session, 3, 2), class="pull-right")
    
    # Place session in the right side of the navigation bar
    pos <- head(pos, -3)
    pluck(ui, !!!pos)[[3]] <- tags$ul(
        class="nav navbar-nav pull-right", session)
    return(ui)
}

# Add loading status in navigation bar
#' @importFrom shiny tagAppendChildren
.addLoadingStatus <- function(ui) {
    loading <- conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                icon("circle-notch", "fa-spin"))
    loading$name <- "a"
    loading <- tags$li(loading)
    
    pos <- .traceInList(ui, "session")
    pluck(ui, !!!head(pos, -5)) <- tagList(
        loading, pluck(ui, !!!head(pos, -5), 1))
    return(ui)
}

#' @importFrom shiny navbarMenu icon span textOutput tags includeScript
#' includeCSS
globalUI <- function(elems, idList, expire) {
    elemClasses       <- sapply(lapply(elems, class), "[[", 1)
    hasSimilarPerts   <- "similarPerturbations" %in% elemClasses
    hasTargetingDrugs <- "targetingDrugs" %in% elemClasses
    showTwoKindPlot   <- hasSimilarPerts && hasTargetingDrugs
    
    ui <- .prepareNavPage(
        id="tab",
        # a non-dropdown tab needs to be selected (bug)
        # https://github.com/rstudio/shiny/issues/3519
        selected="Help",
        navbarMenu("Load", icon=icon("table"),
                   "Differential gene expression data",
                   .diffExprLoadUI(idList$diffExpr),
                   .diffExprENCODEloaderUI(idList$encode),
                   "----",
                   .cmapDataLoaderUI(idList$cmap, globalUI=TRUE)),
        navbarMenu("Analyse", icon=icon("cogs"),
                   .rankSimilarPerturbationsUI(idList$rankPerts, elems, elems),
                   .drugSetEnrichmentAnalyserUI(idList$drugSet, elems, elems)),
        navbarMenu("Visualise", icon=icon("chart-bar"),
                   .dataPlotterUI(idList$data, elems),
                   # .targetingDrugsVSsimilarPerturbationsPlotterUI(
                   #     idList$comparePlot),
                   .datasetComparisonUI(idList$compare, elems),
                   .metadataViewerUI(idList$metadata)),
        tabPanel("Help", icon=icon("question-circle"), "Hello!"),
        navbarMenu(span("Session", span(class="badge",
                                        textOutput("token", inline=TRUE))),
                   icon=icon("compass"), menuName="session"))
    ui <- .modifySessionUI(ui, expire=expire)
    ui <- .addLoadingStatus(ui)
    
    # Add JS and CSS in header
    header <- tags$head(
        includeScript(system.file("shiny", "www", "cTRAP.js", package="cTRAP")),
        includeCSS(system.file("shiny", "www", "cTRAP.css", package="cTRAP")))
    ui <- tagList(header, ui)
    return(ui)
}

#' @importFrom shiny downloadHandler renderText req
.sessionManagementServer <- function(input, output, session, appData) {
    # Show welcome screen when no token is set (e.g. new cTRAP sessions)
    observe({
        if (!is.null(appData$token)) return(NULL)
        showModal(.prepareSessionModal("Welcome to cTRAP!", footer=NULL))
    })
    
    # Create new session
    observeEvent(input$createSession, {
        appData$elems <- NULL
        appData$token <- .createToken()
        removeModal()
    })
    
    # Update token badge
    output$token <- renderText({
        token <- appData$token
        if (is.null(token)) token <- "?"
        return(token)
    })
    
    # Load session based on a token
    observeEvent(input$loadToken, {
        token <- isolate(input$token)
        if (dir.exists(token)) {
            appData$elems <- readRDS(file.path(token, "session.rds"))
            appData$token <- token
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
            appData$elems <- data
            appData$token <- token <- .createToken()
            removeModal()
            .saveSession(data, token)
        }
    })
    
    observeEvent(input$loadSessionModal, {
        modal <- .prepareSessionModal(createSession=TRUE, easyClose=TRUE)
        showModal(modal)
    })
    
    # Notify when copying token
    observeEvent(input$copyToken, {
        msg <- tagList("Token", span(class="badge", appData$token),
                       "copied to your clipboard!")
        showNotification(msg, duration=3, closeButton=FALSE, type="message")
    })
    
    # Download objects in current session in a single RDS file
    output$downloadSession <- downloadHandler(
        filename=function() paste0("cTRAP-", appData$token, ".rds"),
        content=function(file) saveRDS(appData$elems, file))
}

.newDataNotification <- function(names, total, ..., type="message",
                                 auto=FALSE) {
    plural <- ifelse(total == 1, "", "s")
    totalTxt <- sprintf("Total: %s dataset%s", total, plural)
    message(sprintf("  -> %s (%s)",
                    paste(paste(names, collapse=" + "), "loaded"),
                    tolower(totalTxt)))
    
    len <- length(names)
    auto <- ifelse(auto, "automatically ", "")
    head <- "New %sloaded dataset:"
    if (len != 1) head <- paste(length(names), "new %sloaded datasets:")
    head <- sprintf(head, auto)
    
    names <- do.call(tags$ul, lapply(names, tags$li))
    showNotification(tagList(tags$b(head), names, totalTxt), type=type, ...)
}

# Continually check in the background to load new RDS files
#' @importFrom shiny reactivePoll
.loadDataFromLocalRdsServer <- function(input, output, session, appData) {
    checkNewRDSfiles <- function(path=".") {
        if (is.null(path)) return(NULL)
        res <- list.files(path, "\\.rds$", ignore.case=TRUE)
        res <- res[res != "session.rds"]
        if (length(res) == 0) res <- NULL
        res <- file.path(path, res)
        return(res)
    }
    
    getNewRDSfiles <- reactivePoll(
        5000, session,
        checkFunc=function() checkNewRDSfiles(appData$token),
        valueFunc=function() checkNewRDSfiles(appData$token))
    
    observe({
        req(getNewRDSfiles())
        
        elems <- isolate(appData$elems)
        token <- isolate(appData$token)
        
        added <- character(0)
        for (i in getNewRDSfiles()) {
            message("Adding data from ", i, "...")
            obj <- try(readRDS(i), silent=TRUE)
            if (is(obj, "try-error")) {
                warning(obj)
                return(NULL)
            }
            
            # Check if file was expected and replace it accordingly
            expected  <- .filterDatasetsByClass(elems, "expected")
            fileMatch <- match(i, sapply(expected, "[[", "outputFile"))
            if (!is.na(fileMatch)) {
                id <- names(expected)[[fileMatch]]
                message(sprintf("Replacing expected dataset '%s'...", id))
                attr(obj, "formInput") <- attr(elems[[id]], "formInput")
                elems[[id]] <- obj
                added <- c(added, id)
            } else {
                elems <- .addToList(elems, obj)
                added <- c(added, tail(names(elems), 1))
            }
            unlink(i)
        }
        if (length(added) == 0) return(NULL)
        appData$elems <- elems
        .newDataNotification(added, length(elems), duration=NULL, auto=TRUE)
        .saveSession(elems, token)
    })
}

# Update data shared across the app
updateAppData <- function(appData, x) {
    observe({
        obj <- x()
        elems <- .addToList(isolate(appData$elems), obj)
        appData$elems <- elems
        
        dataset <- tail(names(elems), 1)
        if (is(obj, "expected")) {
            showNotification(
                sprintf("'%s' is being calculated and will be loaded when",
                        "ready", dataset),
                type="warning")
        } else {
            .newDataNotification(dataset, length(elems), type="default")
        }
        .saveSession(elems, isolate(appData$token))
    })
}

#' Complete visual interface with support for sessions
#' 
#' Optimised to run in ShinyProxy with Celery/Flower backend with argument
#' \code{shinyproxy = TRUE}.
#'
#' @param ... Objects
#' @param commonPath Character: path where to store data common to all sessions
#' @param expire Character: days until a session expires (message purposes only)
#' @param fileSizeLimitMiB Numeric: file size limit in MiB
#' @param flowerURL Character: Flower REST API's URL (\code{NULL} to avoid using
#' Celery/Flower backend)
#' @inheritParams shiny::runApp
#'
#' @importFrom shiny tagList showModal modalButton modalDialog removeModal
#' reactiveValues tagAppendAttributes showNotification
#'
#' @return Launches result viewer and plotter (returns \code{NULL})
#' @family visual interface functions
#' @export
cTRAP <- function(..., commonPath="data", expire=14, fileSizeLimitMiB=50,
                  flowerURL=NULL, port=getOption("shiny.port"),
                  host=getOption("shiny.host", "127.0.0.1")) {
    .setFileSizeLimit(fileSizeLimitMiB)
    elems <- .prepareEllipsis(...)
    
    # if in ShinyProxy, use Celery/Flower backend via floweRy
    if (!is.null(flowerURL)) {
        if (!require(floweRy)) remotes::install_github("nuno-agostinho/floweRy")
        options(flowerURL=flowerURL)
        flower <- TRUE
    } else {
        flower <- FALSE
    }
    
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
    
    # Get common data from specific folder
    data <- function(x, path=commonPath) file.path(path, x)
    
    server <- function(input, output, session) {
        appData       <- reactiveValues()
        appData$elems <- elems
        elems <- reactive(appData$elems)
        
        # load data
        diffExpr <- .diffExprLoadServer(idList$diffExpr, elems)
        updateAppData(appData, diffExpr)

        encodeDiffExpr <- .diffExprENCODEloaderServer(
            idList$encode, globalUI=TRUE, path=reactive(appData$token),
            metadata=downloadENCODEknockdownMetadata(
                file=data("ENCODEmetadata.rds")))
        updateAppData(appData, encodeDiffExpr)

        cmapData <- .cmapDataLoaderServer(
            idList$cmap, globalUI=TRUE, tab=reactive(session$input$tab),
            metadata=data("cmapMetadata.txt"),
            zscores=data("cmapZscores.gctx"),
            geneInfo=data("cmapGeneInfo.txt"),
            compoundInfo=data("cmapCompoundInfo.txt"))
        updateAppData(appData, cmapData)

        # analyse
        ranking <- .rankSimilarPerturbationsServer(
            idList$rankPerts, elems, elems, globalUI=TRUE, flower=flower,
            token=reactive(appData$token))
        updateAppData(appData, ranking)

        # .drugSetEnrichmentAnalyserServer(idList$drugSet, elems, elems)

        # visualise
        .dataPlotterServer(idList$data, elems)
        .targetingDrugsVSsimilarPerturbationsPlotterServer(
            idList$comparePlot, elems)
        .datasetComparisonServer(idList$compare, elems)
        .metadataViewerServer(idList$metadata, elems)

        .sessionManagementServer(input, output, session, appData)
        .loadDataFromLocalRdsServer(input, output, session, appData)
    }
    app <- runApp(shinyApp(ui, server), port=port, host=host)
    return(app)
}
