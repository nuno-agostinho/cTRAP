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
    
    # Unique names for each data element
    uniqName <- make.unique(c(names(x), name))
    name <- uniqName[[length(uniqName)]]
    
    x[[name]] <- data
    return(x)
}

# Save session data in token-named directory
#' @importFrom qs qsave
.saveSession <- function(data, token) {
    if (is.null(token) || is.null(data)) return(NULL)
    if (!dir.exists(token)) dir.create(token)
    sessionQS <- file.path(token, "session.qs")
    qsave(data, sessionQS)
    message("     Session saved to ", sessionQS)
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
                              paste("Session expires in", expire, "days"))
    }
    
    copyTokenButton <- actionLink("copyToken", onclick="copyToken()", tagList(
        "Copy session token to clipboard", expireTxt))
    pluck(session, 3, 2, 3) <- tagList(
        tags$li(role="presentation", copyTokenButton),
        tags$li(role="presentation", class="divider"),
        tags$li(role="presentation",
                downloadLink("downloadSession", "Download session data")),
        tags$li(role="presentation",
                actionLink("loadSessionModal", "Load another session...")))
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
        .metadataViewerUI(idList$metadata, icon=icon("layer-group")),
        navbarMenu("Load", icon=icon("table"),
                   "Differential gene expression data",
                   .diffExprLoadUI(idList$diffExpr),
                   .diffExprENCODEloaderUI(idList$encode),
                   "----",
                   .cmapDataLoaderUI(idList$cmap, globalUI=TRUE)),
        navbarMenu("Analyse", icon=icon("cogs"),
                   .rankSimilarPerturbationsUI(idList$rankPerts),
                   .predictTargetingDrugsUI(idList$predictDrugs)),
        navbarMenu("Visualise", icon=icon("chart-bar"),
                   .dataPlotterUI(idList$data),
                   .datasetComparisonUI(idList$compare),
                   .targetingDrugsVSsimilarPerturbationsPlotterUI(
                       idList$comparePlot),
                   .drugSetEnrichmentAnalyserUI(idList$drugSet)),
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

# Set app data and add tags
.setAppData <- function(appData, elems) {
    appData$elems <- .addDatasetTags(elems)
    return(appData)
}

#' @importFrom shiny downloadHandler renderText req
#' @importFrom qs qread
#' @importFrom utils packageVersion
.sessionManagementServer <- function(input, output, session, appData) {
    # Show welcome screen when no token is set (e.g. new cTRAP sessions)
    observe({
        if (!is.null(appData$token)) return(NULL)
        title <- sprintf("Welcome to cTRAP %s!", packageVersion("cTRAP"))
        showModal(.prepareSessionModal(title, footer=NULL))
    })
    
    # Create new session
    observeEvent(input$createSession, {
        .setAppData(appData, NULL)
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
            file <- file.path(token, "session")
            rds  <- paste0(file, ".rds")
            qs   <- paste0(file, ".qs")
            if (file.exists(qs)) {
                .setAppData(appData, qread(qs))
            } else if (file.exists(rds)) {
                .setAppData(appData, readRDS(rds))
            }
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
        
        data <- tryCatch(readRDS(file$datapath), error=function(e) e)
        if (is(data, "error")) {
            showNotification(paste("Error loading data:", data),
                             type="error")
        } else {
            .setAppData(appData, data)
            appData$token <- token <- .createToken()
            removeModal()
            .saveSession(data, token)
        }
    })
    
    observeEvent(input$loadSessionModal, {
        title <- sprintf("Welcome to cTRAP %s!", packageVersion("cTRAP"))
        modal <- .prepareSessionModal(title, createSession=TRUE, easyClose=TRUE)
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

.newDataNotification <- function(names, total, expected=NULL, ...,
                                 type="message", auto=FALSE) {
    plural <- ifelse(total == 1, "", "s")
    totalTxt <- sprintf("Total: %s dataset%s", total, plural)
    if (!is.null(expected) && expected > 0) {
        totalTxt <- paste(totalTxt, sprintf("(%s running)", expected))
    }
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

# Continually check if the output files from Celery tasks are ready to be loaded
#' @importFrom shiny reactivePoll
.loadCeleryOutputServer <- function(input, output, session, appData) {
    getExpectedAppDataTasks <- function(elems) {
        if (is.null(elems) || length(elems) == 0) return(NULL)
        expected <- .filterDatasetsByClass(elems, "expected")
        return(expected)
    }
    
    checkExpectedCeleryTasks <- function(elems) {
        expectedAppTasks <- getExpectedAppDataTasks(elems)
        if (is.null(expectedAppTasks)) return(NULL)
         
        expectedTaskID <- sapply(expectedAppTasks, "[[", "task-id")
        if (length(expectedTaskID) == 0) return(NULL)
        
        tasks <- taskList()
        tasks <- tasks[tasks$uuid %in% expectedTaskID, ]
        if (is.null(tasks)) return(NULL)
        return(tasks)
    }
    
    getExpectedCeleryTasks <- reactivePoll(
        5000, session,
        checkFunc=function() checkExpectedCeleryTasks(appData$elems),
        valueFunc=function() checkExpectedCeleryTasks(appData$elems))
    
    observe({
        tasks    <- req(getExpectedCeleryTasks())
        
        elems    <- isolate(appData$elems)
        token    <- isolate(appData$token)
        
        added <- character(0)
        updatedState <- FALSE
        for (id in names( getExpectedAppDataTasks(elems) )) {
            outputFile <- elems[[id]]$outputFile
            if (file.exists(outputFile)) {
                # Read output RDS file
                message(sprintf("Updating '%s' with data from %s...",
                                id, outputFile))
                obj <- try(readRDS(outputFile), silent=TRUE)
                if (is(obj, "try-error")) {
                    warning(obj)
                    return(NULL)
                }
                
                # Replace data accordingly
                attr(obj, "formInput") <- attr(elems[[id]], "formInput")
                elems[[id]] <- obj
                added <- c(added, id)
                
                # Remove output file
                unlink(outputFile)
            } else {
                # Skip if task state is not found
                if (!"state" %in% colnames(tasks)) next
                
                # Update state of tasks if needed
                taskID   <- elems[[id]][["task-id"]]
                matched  <- tasks$uuid == taskID
                if (!any(matched)) {
                    message("cTRAP task not found in Celery/Flower history...")
                    newState <- "Not found"
                } else {
                    newState <- tasks[matched, "state"]
                }
                newState <- tolower(newState)
                oldState <- tolower(elems[[id]]$state)
                
                if (newState != oldState) {
                    message(sprintf("Updating %s from '%s' to '%s'...",
                                    id, oldState, newState))
                    elems[[id]]$state <- capitalize(newState)
                    updatedState <- TRUE
                }
            }
        }
        newDatasets <- length(added) > 0
        if (!newDatasets && !updatedState) {
            return(NULL)
        } else if (newDatasets) {
            expected <- length(.filterDatasetsByClass(elems, "expected"))
            .newDataNotification(added, length(elems) - expected, duration=30,
                                 auto=TRUE)
        }
        .setAppData(appData, elems)
        .saveSession(elems, token)
    })
}

# Update data shared across the app
updateAppData <- function(appData, x) {
    observe({
        obj <- x()
        elems <- .addToList(isolate(appData$elems), obj)
        .setAppData(appData, elems)
        token <- isolate(appData$token)
        
        dataset <- tail(names(elems), 1)
        if (is(obj, "expected")) {
            msg <- tagList(
                sprintf("'%s' is being calculated", dataset),
                tags$br(), tags$br(), tags$b("You may close the browser"),
                "and use the session token", span(class="badge", token),
                "to load your data later")
            showNotification(msg, type="warning", duration=10)
        } else {
            expected <- length(.filterDatasetsByClass(elems, "expected"))
            .newDataNotification(dataset, length(elems) - expected,
                                 type="default")
        }
        .saveSession(elems, token)
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
#' @importFrom utils getFromNamespace assignInNamespace
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
        # if (!requireNamespace("floweRy")) {
        #     remotes::install_github("nuno-agostinho/floweRy")
        # }
        options(floweRy.url=flowerURL)
        flower <- TRUE
    } else {
        flower <- FALSE
    }
    
    # Avoid large JSON response from DT: github.com/rstudio/DT/issues/504
    dt_mod <- getFromNamespace("dataTablesFilter", "DT")
    dt_rows_all_line <- grep("DT_rows_all = iAll", body(dt_mod))

    if (length(dt_rows_all_line) == 1) {
        mod <- gsub("DT_rows_all = iAll", "DT_rows_all = iCurrent", fixed=TRUE,
                    body(dt_mod)[dt_rows_all_line])
        body(dt_mod)[[dt_rows_all_line]] <- parse(text=mod)[[1]]
        assignInNamespace("dataTablesFilter", dt_mod, "DT")
    }

    idList              <- list()
    idList$diffExpr     <- "diffExprLoader"
    idList$encode       <- "encodeDataLoader"
    idList$cmap         <- "cmapDataLoader"
    idList$compare      <- "datasetComparison"
    idList$comparePlot  <- "comparePlotter"
    idList$data         <- "dataPlotter"
    idList$metadata     <- "metadataViewer"
    idList$rankPerts    <- "rankPerts"
    idList$predictDrugs <- "predictDrugs"
    idList$drugSet      <- "drugSetAnalyser"
    ui <- globalUI(elems, idList, expire)
    
    # Get common data from specific folder
    loadCommonData <- function(x, path=commonPath) file.path(path, x)
    
    server <- function(input, output, session) {
        appData       <- reactiveValues()
        .setAppData(appData, elems)
        elems <- reactive(appData$elems)
        
        # load data
        diffExpr <- .diffExprLoadServer(idList$diffExpr, elems)
        updateAppData(appData, diffExpr)

        encodeDiffExpr <- .diffExprENCODEloaderServer(
            idList$encode, globalUI=TRUE, path=reactive(appData$token),
            metadata=downloadENCODEknockdownMetadata(
                file=loadCommonData("ENCODEmetadata.rds")))
        updateAppData(appData, encodeDiffExpr)

        cmapData <- .cmapDataLoaderServer(
            idList$cmap, globalUI=TRUE, tab=reactive(session$input$tab),
            metadata=loadCommonData("cmapMetadata.txt"),
            zscores=loadCommonData("cmapZscores.gctx"),
            geneInfo=loadCommonData("cmapGeneInfo.txt"),
            compoundInfo=loadCommonData("cmapCompoundInfo.txt"))
        updateAppData(appData, cmapData)

        # analyse
        ranking <- .rankSimilarPerturbationsServer(
            idList$rankPerts, elems, globalUI=TRUE, flower=flower,
            token=reactive(appData$token))
        updateAppData(appData, ranking)
        
        predicted <- .predictTargetingDrugsServer(
            idList$predictDrugs, elems, globalUI=TRUE, flower=flower,
            path=commonPath, token=reactive(appData$token))
        updateAppData(appData, predicted)

        .drugSetEnrichmentAnalyserServer(idList$drugSet, elems, path=commonPath)

        # visualise
        .dataPlotterServer(idList$data, elems)
        .targetingDrugsVSsimilarPerturbationsPlotterServer(
            idList$comparePlot, elems)
        .datasetComparisonServer(idList$compare, elems)
        .metadataViewerServer(idList$metadata, elems)

        .sessionManagementServer(input, output, session, appData)
        if (flower) .loadCeleryOutputServer(input, output, session, appData)
    }
    app <- runApp(shinyApp(ui, server), port=port, host=host)
    return(app)
}
