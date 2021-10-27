# Temprarily include floweRy functions until approved in CRAN (?)
# Source code: https://github.com/nuno-agostinho/floweRy

getFlowerURL <- function() { getOption("floweRy.url", "http://localhost:5555") }

#' @importFrom httr status_code stop_for_status
simplifyAPIerror <- function(res, errors) {
    msg <- errors[[as.character(status_code(res))]]
    if (!is.null(msg)) {
        stop(msg)
    } else {
        stop_for_status(res)
    }
}

#' @importFrom httr GET content
getInfo <- function(type, ..., url=getFlowerURL(), errors=NULL) {
    res <- GET(url=file.path(url, "api", type), query=list(...))
    simplifyAPIerror(res, errors)
    return(content(res))
}

taskList <- function(limit=NULL, offset=NULL,
                     sort_by=c("name", "state", "received", "started"),
                     workername=NULL, taskname=NULL, state=NULL,
                     received_start=NULL, received_end=NULL, table=TRUE,
                     url=getFlowerURL()) {
    sort_by <- match.arg(sort_by)
    res <- getInfo(type="tasks", limit=limit, offset=offset, sort_by=sort_by,
                   workername=workername, taskname=taskname, state=state,
                   received_start=received_start, received_end=received_end,
                   url=url)
    
    timestamps <- c("received", "started", "succeeded", "timestamp", "revoked")
    convertTime <- function(x) as.POSIXct(x, origin="1970-01-01")
    
    for (task in names(res)) {
        # Replace NULL with NA
        nulls <- vapply(res[[task]], is.null, logical(1))
        res[[task]][nulls] <- NA
        
        # Replace Unix timestamp with formatted time
        res[[task]][timestamps] <- lapply(res[[task]][timestamps], convertTime)
    }
    
    # Create table
    if (table) {
        cols <- unique(unlist(lapply(res, names)))
        df   <- data.frame(matrix(ncol=length(cols), nrow=length(res),
                                  dimnames=list(names(res), cols)))
        for (col in cols) {
            values <- sapply(res, "[[", col)
            if (col %in% timestamps) values <- convertTime(values)
            df[[col]] <- values
        }
        res <- df
    }
    return(res)
}

#' @importFrom httr POST content
runTask <- function(type="apply", task=NULL,
                    args=NULL, kwargs=NULL, options=NULL, url=getFlowerURL(),
                    errors=NULL) {
    path <- file.path("api", "task", type, task)
    
    convert2list <- function (x) if (!is.null(x) && !is.list(x)) list(x) else x
    args    <- convert2list(args)
    kwargs  <- convert2list(kwargs)
    options <- convert2list(options)
    
    body <- list(args=args, kwargs=kwargs, options=options)
    res  <- POST(url=file.path(url, path), body=body, encode="json")
    simplifyAPIerror(res, errors)
    return(content(res))
}

taskAsyncApply <- function(task, ..., kwargs=NULL, options=NULL,
                           url=getFlowerURL()) {
    args <- list(...)
    if (length(args) == 0) args <- NULL
    
    errors <- list("404"="unknown task")
    runTask(type="async-apply", task=task, args=prepareArgs(...), kwargs=kwargs,
            options=options, url=url, errors=errors)
}
