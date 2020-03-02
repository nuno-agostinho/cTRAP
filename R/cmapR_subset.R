# Functions to parse GCTX files based on cmapR (R package that is not available
# on CRAN nor Bioconductor): https://github.com/cmap/cmapR

# io.R -------------------------------------------------------------------------

#' An S4 class to represent a GCT object
#'
#' @slot mat a numeric matrix
#' @slot rid a character vector of row ids
#' @slot cid a character vector of column ids
#' @slot rdesc a \code{data.frame} of row descriptors
#' @slot rdesc a \code{data.frame} of column descriptors
#' @slot src a character indicating the source (usually file path) of the data
#'
#' @description The GCT class serves to represent annotated
#'   matrices. The \code{mat} slot contains said data and the
#'   \code{rdesc} and \code{cdesc} slots contain data frames with
#'   annotations about the rows and columns, respectively
#'
#' @seealso \url{http://clue.io/help} for more information on the GCT format
#'
#' @inherit readGctxMeta source
#' @keywords internal
setClass("GCT", representation(
    mat = "matrix", rid = "character", cid = "character", rdesc = "data.frame",
    cdesc = "data.frame", version = "character", src = "character"))

#' Adjust the data types for columns of a meta data frame
#'
#' @description GCT(X) parsing initially returns data frames
#'   of row and column descriptors where all columns are of
#'   type character. This is inconvenient for analysis, so
#'   the goal of this function is to try and guess the
#'   appropriate data type for each column.
#'
#' @param meta a data.frame
#'
#' @details This is a low-level helper function
#'   which most users will not need to access directly
#'
#' @return meta the same data frame with (potentially) adjusted column types
#' @keywords internal
#'
#' @family GCTX parsing functions
#' @inherit readGctxMeta source
fix.datatypes <- function(meta) {
    for (field.name in names(meta)) {
        field <- meta[[field.name]]
        # Check if string contains numeric values
        field.as.numeric <- suppressWarnings(as.numeric(field))
        if (!any(is.na(field.as.numeric))) {
            field <- field.as.numeric
        }
        if (is.numeric(field)) {
            # Check if floats are actually integers
            field.as.integer <- suppressWarnings(as.integer(field))
            if (!any(is.na(field.as.integer))) {
                diffs <- field - field.as.integer
                if (all(diffs == 0)) field <- field.as.integer
            }
        }
        meta[[field.name]] <- field
    }
    return(meta)
}

#' Parse row or column metadata from GCTX files
#'
#' @param gctx_path the path to the GCTX file
#' @param dimension which metadata to read (row or column)
#' @param ids a character vector of a subset of row/column ids
#'   for which to read the metadata
#' @param set_annot_rownames a boolean indicating whether to set the
#'   \code{rownames} attribute of the returned \code{data.frame} to the
#'   corresponding row/column ids.
#'
#' @importFrom rhdf5 h5read
#'
#' @return a \code{data.frame} of metadata
#' @keywords internal
#'
#' @source \url{https://github.com/cmap/cmapR}
#' @family GCTX parsing functions
readGctxMeta <- function(gctx_path, dimension="row", ids=NULL,
                           set_annot_rownames=TRUE) {
    if (!file.exists(gctx_path)) stop(paste(gctx_path, "not found"))
    if (dimension=="column") dimension <- "col"
    if (!(dimension %in% c("row", "col")))
        stop("dimension can be either row or col")

    if (dimension == "row")
        name <- "0/META/ROW"
    else
        name <- "0/META/COL"

    raw_annots <- h5read(gctx_path, name=name) # Returns a list
    fields <- names(raw_annots)
    # Empty data frame of the correct dimensions
    annots <-  data.frame(matrix(nrow=length(raw_annots[[fields[1]]]),
                                 ncol=length(fields)))
    names(annots) <-  fields
    # Loop through each field and fill the annots data.frame
    for (i in seq(length(fields))) {
        field <- fields[i]
        # Remove any trailing spaces and cast as vector
        annots[,i] <- as.vector(gsub("\\s*$", "", raw_annots[[field]],
                                     perl=TRUE))
    }
    annots <- fix.datatypes(annots)
    # Subset to the provided set of identifiers, if any
    if (is.null(ids))
        ids <- as.character(annots$id)
    else
        ids <- ids

    annots <- subsetToIds(annots, ids)
    annots$id <- as.character(annots$id)
    if (set_annot_rownames) rownames(annots) <- annots$id
    return(annots)
}

#' Read GCTX row or column ids
#'
#' @param gctx_path path to the GCTX file
#' @param dimension which ids to read (row or column)
#'
#' @return a character vector of row or column ids from the provided file
#' @keywords internal
#'
#' @inherit readGctxMeta source
#' @family GCTX parsing functions
readGctxIds <- function(gctx_path, dimension="row") {
    if (!file.exists(gctx_path)) {
        stop(paste(gctx_path, "not found"))
    }
    if (dimension=="column") dimension <- "col"
    if (!(dimension %in% c("row", "col"))) {
        stop("dimension can be either row or col")
    }
    if (dimension == "row") {
        name <- "0/META/ROW/id"
    } else {
        name <- "0/META/COL/id"
    }
    ids <- gsub("\\s*$", "", h5read(gctx_path, name=name), perl=TRUE)
    ids <- as.character(ids)
    return(ids)
}

#' Return a subset of requested GCTX row/column ids
#' out of the universe of all ids
#'
#' @details This is a low-level helper function
#'   which most users will not need to access directly
#'
#' @param ids vector of requested ids. If \code{NULL}, no
#'   subsetting is performed
#' @param all_ids vector of universe of ids
#' @param type flag indicating the type of ids being processed
#'
#' @return a list with the following elements
#'  \code{ids}: a character vector of the processed ids
#'  \code{idx}: an integer list of their corresponding indices in \code{all_ids}
#'
#' @inherit readGctxMeta source
#' @family GCTX parsing functions
#' @keywords internal
processIds <- function(ids, all_ids, type="rid") {
    if (!is.null(ids)) {
        if (is.numeric(ids)) {
            # Check if numeric
            idx <- ids
            is_invalid_idx <- (idx > length(all_ids)) | (idx <= 0)
            invalid_idx <- idx[is_invalid_idx]
            if (all(is_invalid_idx)) {
                stop(paste("none of the requested", type,
                           "indices were found in the dataset"))
            }
            if (any(is_invalid_idx)) {
                # Requested indices are outside of the possible range
                warning(paste(
                    "the following ", type,
                    " were are outside possible range and will be ignored:\n",
                    paste(invalid_idx, collapse="\n"), sep=""))
            }
            idx <- idx[!is_invalid_idx]
        } else {
            # Assume it is a character
            idx <- match(ids, all_ids)
            if (all(is.na(idx))) {
                stop(paste("none of the requested", type,
                           "were found in the dataset"))
            }
            if (any(is.na(idx))) {
                ids_not_found <- ids[is.na(idx)]
                warning(paste("multiple", type,
                              "were not found and will be ignored"))
            }
            idx <- idx[!is.na(idx)]
        }
    } else {
        # NULL IDs: just return an index vector along all_ids
        idx <- seq_along(all_ids)
    }
    # Subset the character identifiers to the ones we want
    id_keep <- as.character(all_ids[idx])
    return(list(idx=idx, ids=id_keep))
}

# Define the initialization method for the GCT class
setMethod("initialize", signature = "GCT", definition = function(
    .Object, mat=NULL, rdesc=NULL, cdesc=NULL, src=NULL, rid=NULL, cid=NULL,
    set_annot_rownames=FALSE, matrix_only=FALSE, verbose=TRUE) {

    if (!is.null(src)) src <- as.character(src)
    time <- Sys.time()

    # Use both matrix and annotation (if supplied)
    if (!is.null(mat)) {
        .Object@mat <- mat
        .Object@rid <- rownames(mat)

        if (!is.null(cid))
            .Object@cid <- cid
        else
            .Object@cid <- colnames(mat)
    }
    if (!is.null(rdesc)) .Object@rdesc <- rdesc
    if (!is.null(cdesc)) {
        .Object@cdesc <- cdesc
    } else if (!is.null(src)) {
        if (! (grepl(".gct$", src) || grepl(".gctx$", src) ))
            stop("either a .gct or .gctx file must be given")
        if (grepl(".gct$", src)) {
            if ( ! is.null(rid) || !is.null(cid) ) {
                msg <- paste(
                    "Ignoring rid and cid values (which are only relevant for",
                    ".gctx files)")
                warning(msg)
            }
            .Object@src = src
            # Get the GCT version from the first line
            .Object@version = scan(src, what = "", nlines = 1, sep = "\t",
                                   quiet = TRUE)[1]
            # Get matrix dimensions from the second line
            dimensions = scan(src, what = double(0), nlines = 1, skip = 1,
                              sep = "\t", quiet = TRUE)
            nrmat = dimensions[1]
            ncmat = dimensions[2]
            if (length(dimensions)==4) {
                if (verbose) message("Parsing as GCT v1.3...")
                nrhd <- dimensions[3]
                nchd <- dimensions[4]
            } else {
                if (verbose) message("Parsing as GCT v1.2...")
                nrhd <- 0
                nchd <- 0
            }
            if (verbose) {
                message(paste(src, nrmat, "rows,", ncmat, "cols,", nrhd,
                              "row descriptors,", nchd, "col descriptors"))
            }
            # Read header
            header <- scan(src, what="", nlines=1, skip=2, sep="\t", quote=NULL,
                           quiet=TRUE)
            # Prepare row and column identifiers from header
            if ( nrhd > 0 ) {
                rhd <- header[2:(nrhd+1)]
                cid <- header[-(nrhd+1):-1]
                col_offset <- 1
            }
            else {
                if (any(grepl("description", header, ignore.case=TRUE))) {
                    # Check description column in v1.2 files
                    col_offset <- 2
                } else {
                    col_offset <- col_offset <- 1
                }
                rhd = NULL
                cid = header[(1+col_offset):length(header)]
            }
            # Read next set of headers (column annotations), shape into a matrix
            if ( nchd > 0 ) {
                header = scan(src, what = "", nlines = nchd, skip = 3,
                              sep = "\t", quote = NULL, quiet = TRUE)
                header = matrix(header, nrow = nchd,
                                ncol = ncmat + nrhd + 1, byrow = TRUE)
                # Extract column header and descriptions
                chd = header[,1]
                cdesc = header[,-(nrhd+1):-1]
                # Transpose in the case there is only one column annotation
                if ( nchd == 1 ) cdesc = t(cdesc)
            } else {
                chd = NULL
                cdesc = data.frame()
            }
            # Read data matrix and row descriptions, shape into a matrix
            mat = scan(src, what = "", nlines = nrmat,
                       skip = 3 + nchd, sep = "\t", quote = NULL, quiet = TRUE)
            mat = matrix(mat, nrow = nrmat, ncol = ncmat + nrhd + col_offset,
                         byrow = TRUE)
            # if (verbose) message(paste(dim(mat), collapse="\t"))
            # Extract the row identifiers and descriptions, and data matrix
            rid = mat[,1]
            if ( nrhd > 0 ) {
                # Use as.matrix for when there is only one row annotation
                rdesc = as.matrix(mat[,2:(nrhd + 1)])
                mat = matrix(as.numeric(mat[,-(nrhd + 1):-1]),
                             nrow = nrmat, ncol = ncmat)
            }
            else {
                rdesc = data.frame()
                mat = matrix(as.numeric(mat[, (1+col_offset):ncol(mat)]),
                             nrow = nrmat, ncol = ncmat)
            }
            # Assign names to the data matrix and the row/column descriptions
            # if (verbose) message(paste(dim(mat), collapse="\t"))
            dimnames(mat) = list(rid, cid)
            if ( nrhd > 0 ) {
                dimnames(rdesc) = list(rid,rhd)
                rdesc = as.data.frame(rdesc, stringsAsFactors = FALSE)
            }
            if ( nchd > 0 ) {
                cdesc = t(cdesc)
                dimnames(cdesc) = list(cid,chd)
                cdesc = as.data.frame(cdesc, stringsAsFactors = FALSE)
            }
            # Assign to the GCT slots
            .Object@mat = mat
            .Object@rid = rownames(mat)
            .Object@cid = colnames(mat)
            if (!matrix_only) {
                # Return annotations as well as matrix
                .Object@rdesc = fix.datatypes(rdesc)
                .Object@cdesc = fix.datatypes(cdesc)
                # Add column identifiers to rdesc and cdesc
                .Object@rdesc$id <- rownames(.Object@rdesc)
                .Object@cdesc$id <- rownames(.Object@cdesc)
            }
        } else {
            if (verbose) message(sprintf("Reading %s...", src))
            .Object@src = src
            # Get all row and column identifiers
            all_rid <- readGctxIds(src, dimension="row")
            all_cid <- readGctxIds(src, dimension="col")
            # Only read rows/columns specified by rid/cid (if available)
            processed_rids <- processIds(rid, all_rid, type="rid")
            processed_cids <- processIds(cid, all_cid, type="cid")
            # Read data matrix
            .Object@mat <- h5read(src, name="0/DATA/0/matrix",
                                  index=list(processed_rids$idx,
                                             processed_cids$idx))
            # Set the row and column identifiers, casting as characters
            .Object@rid <- processed_rids$ids
            .Object@cid <- processed_cids$ids
            rownames(.Object@mat) <- processed_rids$ids
            colnames(.Object@mat) <- processed_cids$ids
            # Get metadata
            if (!matrix_only) {
                .Object@rdesc <- readGctxMeta(
                    src, dimension="row", ids=processed_rids$ids,
                    set_annot_rownames=set_annot_rownames)
                .Object@cdesc <- readGctxMeta(
                    src, dimension="col", ids=processed_cids$ids,
                    set_annot_rownames=set_annot_rownames)
            } else {
                .Object@rdesc <- data.frame(id=.Object@rid,
                                            stringsAsFactors = FALSE)
                .Object@cdesc <- data.frame(id=.Object@cid,
                                            stringsAsFactors = FALSE)
            }
            closeOpenHandles()
            diff <- Sys.time() - time
            if (verbose) {
                msg <- sprintf("Successfully read data from %s in %s", src,
                               format(diff, digits=3))
                message(msg)
            }
        }
    }
    return(.Object)
})

#' Close open handles
#'
#' @importFrom utils packageVersion
#'
#' @return Closes all open identifiers
#' @keywords internal
closeOpenHandles <- function() {
    if(packageVersion('rhdf5') < "2.23.0")
        rhdf5::H5close()
    else
        rhdf5::h5closeAll()
}

# utils.R ----------------------------------------------------------------------

#' Check whether \code{test_names} are columns in the \code{\link{data.frame}}
#' @param test_names a vector of column names to test
#' @param df the \code{\link{data.frame}} to test against
#' @param throw_error boolean indicating whether to throw an error if
#'   any \code{test_names} are not found in \code{df}
#' @return boolean indicating whether or not all \code{test_names} are
#'   columns of \code{df}
#'
#' @inherit readGctxMeta source
#' @keywords internal
checkColnames <- function(test_names, df, throw_error=TRUE) {
    # Check if test_names are valid; throw error if specified
    diffs <- setdiff(test_names, names(df))
    if (length(diffs) > 0) {
        if (throw_error) {
            stop(paste("the following column names are not found in",
                       deparse(substitute(df)), ":",
                       paste(diffs, collapse=" "), "\n"))
        } else {
            return(FALSE)
        }
    } else {
        return(TRUE)
    }
}

#' Do a robust \code{\link{data.frame}} subset to a set of ids
#' @param df \code{\link{data.frame}} to subset
#' @param ids the ids to subset to
#' @return a subset version of \code{df}
#'
#' @inherit readGctxMeta source
#' @keywords internal
subsetToIds <- function(df, ids) {
    checkColnames("id", df)
    newdf <- data.frame(df[match(ids, df$id), ])
    names(newdf) <- names(df)
    return(newdf)
}
