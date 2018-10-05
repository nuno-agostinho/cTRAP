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
#' @seealso \code{\link{parse.gctx}}, \code{\link{write.gctx}},
#' \code{\link{read.gctx.meta}}, \code{\link{read.gctx.ids}}
#' @seealso \link{http://clue.io/help} for more information on the GCT format
#'
#' @source https://github.com/cmap/cmapR
setClass("GCT", representation(
    mat = "matrix", rid = "character", cid = "character", rdesc = "data.frame",
    cdesc = "data.frame", version = "character", src = "character"))

# set up methods for checking GCT validity
setValidity("GCT", function(object) {
    # check whether dimensions of various
    # slots are in sync
    nrows <- nrow(object@mat)
    ncols <- ncol(object@mat)
    if (nrows != length(object@rid)) {
        return("rid must be the same length as number of matrix rows")
    }
    if (ncols != length(object@cid)) {
        return("cid must be the same length as number of matrix columns")
    }
    if (length(object@cid) > length(unique(object@cid))) {
        return("cid must be unique")
    }
    if (length(object@rid) > length(unique(object@rid))) {
        return("rid must be unique")
    }
    if (nrow(object@cdesc) != ncols & nrow(object@cdesc) != 0) {
        return(paste("cdesc must either have 0 rows or the same number of rows",
                     "as matrix has columns"))
    }
    if (nrow(object@rdesc) != nrows & nrow(object@rdesc) != 0) {
        return(paste("rdesc must either have 0 rows or the same number of rows",
                     "as matrix has rows"))
    } else {
        return(TRUE)
    }
})

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
#' @return meta the same data frame with (potentially) adjusted
#'   column types.
#'
#' @examples
#' # meta data table with all character types
#' str(cdesc_char)
#' fixed <- cmapR:::fix.datatypes(cdesc_char)
#' # note how some column classes have changed
#' str(fixed)
#'
#' @family GCTX parsing functions
#' @keywords internal
#'
#' @source https://github.com/cmap/cmapR
fix.datatypes <- function(meta) {
    for (field.name in names(meta)) {
        # get the field values
        field <- meta[[field.name]]
        # check if it's numeric. data may come in as a string
        # but actually contains numeric values. if so, as.numeric
        # will not result in a vector of NA values
        field.as.numeric <- suppressWarnings(as.numeric(field))
        if (!any(is.na(field.as.numeric))) {
            field <- field.as.numeric
        }
        if (is.numeric(field)) {
            # check if it's an integer. data may be floats but
            # if we coerce to an integer and the difference from
            # original values is zero, that means data are actually
            # integers. integer conversion will return NA if there
            # are any issues.
            field.as.integer <- suppressWarnings(as.integer(field))
            if (!any(is.na(field.as.integer))) {
                # integer conversion was fine, lets see if the
                # values are altered
                diffs <- field - field.as.integer
                if (all(diffs == 0)) {
                    # converting to integer didn't change value,
                    # set field to integer values
                    field <- field.as.integer
                }
            }
        }
        # insert back into the annotations
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
#'   \code{rownames} addtribute of the returned \code{data.frame} to
#'   the corresponding row/column ids.
#'
#' @importFrom rhdf5 h5read
#'
#' @return a \code{data.frame} of metadata
#'
#' @examples
#' gct_file <- system.file("extdata", "modzs_n272x978.gctx", package="cmapR")
#' # row meta
#' row_meta <- read.gctx.meta(gct_file)
#' str(row_meta)
#' # column meta
#' col_meta <- read.gctx.meta(gct_file, dimension="column")
#' str(col_meta)
#' # now for only the first 10 ids
#' col_meta_first10 <- read.gctx.meta(gct_file, dimension="column",
#'                                    ids=col_meta$id[1:10])
#' str(col_meta_first10)
#'
#' @source https://github.com/cmap/cmapR
#'
#' @family GCTX parsing functions
#' @export
read.gctx.meta <- function(gctx_path, dimension="row", ids=NULL,
                           set_annot_rownames=TRUE) {
    if (!file.exists(gctx_path)) stop(paste(gctx_path, "does not exist"))
    if (dimension=="column") dimension <- "col"
    if (!(dimension %in% c("row", "col")))
        stop("dimension can be either row or col")

    if (dimension == "row")
        name <- "0/META/ROW"
    else
        name <- "0/META/COL"

    raw_annots <- h5read(gctx_path, name=name) # returns a list
    fields <- names(raw_annots)
    # define an empty data frame of the correct dimensions
    annots <-  data.frame(matrix(nrow=length(raw_annots[[fields[1]]]),
                                 ncol=length(fields)))
    names(annots) <-  fields
    # loop through each field and fill the annots data.frame
    for (i in 1:length(fields)) {
        field <- fields[i]
        # remove any trailing spaces
        # and cast as vector
        annots[,i] <- as.vector(gsub("\\s*$", "", raw_annots[[field]], perl=T))
    }
    annots <- fix.datatypes(annots)
    # subset to the provided set of ids, if given
    if (is.null(ids))
        ids <- as.character(annots$id)
    else
        ids <- ids

    # make sure annots row ordering matches that of ids
    annots <- subset_to_ids(annots, ids)
    annots$id <- as.character(annots$id)
    # use the id field to set the rownames
    if (set_annot_rownames) rownames(annots) <- annots$id
    return(annots)
}

#' Read GCTX row or column ids
#'
#' @param gctx_path path to the GCTX file
#' @param dimension which ids to read (row or column)
#'
#' @return a character vector of row or column ids from the provided file
#'
#' @examples
#' gct_file <- system.file("extdata", "modzs_n272x978.gctx", package="cmapR")
#' # row ids
#' rid <- read.gctx.ids(gct_file)
#' head(rid)
#' # column ids
#' cid <- read.gctx.ids(gct_file, dimension="column")
#' head(cid)
#'
#' @source https://github.com/cmap/cmapR
#'
#' @family GCTX parsing functions
#' @export
read.gctx.ids <- function(gctx_path, dimension="row") {
    if (!file.exists(gctx_path)) {
        stop(paste(gctx_path, "does not exist"))
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
    # remove any spaces
    ids <- gsub("\\s*$", "", h5read(gctx_path, name=name), perl=T)
    # cast as character
    ids <- as.character(ids)
    return(ids)
}

#' Return a subset of requested GCTX row/colum ids
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
#' @examples
#' gct_file <- system.file("extdata", "modzs_n272x978.gctx", package="cmapR")
#' ids <- read.gctx.ids(gct_file)
#' processed_ids <- cmapR:::process_ids(ids[1:10], ids)
#' str(processed_ids)
#'
#' @source
#'
#' @family GCTX parsing functions
#' @keywords internal
process_ids <- function(ids, all_ids, type="rid") {
    if (!is.null(ids)) {
        if (is.numeric(ids)) {
            # is it numeric?
            idx <- ids
            is_invalid_idx <- (idx > length(all_ids)) | (idx <= 0)
            invalid_idx <- idx[is_invalid_idx]
            if (all(is_invalid_idx)) {
                stop(paste("none of the requested", type,
                           "indices were found in the dataset"))
            }
            if (any(is_invalid_idx)) {
                # requested indices are outside of the possible range
                warning(paste(
                    "the following ", type,
                    " were are outside possible range and will be ignored:\n",
                    paste(invalid_idx, collapse="\n"), sep=""))
            }
            idx <- idx[!is_invalid_idx]
        } else {
            # assume its a character
            idx <- match(ids, all_ids)
            if (all(is.na(idx))) {
                stop(paste("none of the requested", type,
                           "were found in the dataset"))
            }
            if (any(is.na(idx))) {
                ids_not_found <- ids[is.na(idx)]
                warning(paste("the following ", type,
                              " were not found and will be ignored:\n",
                              paste(ids_not_found, collapse="\n"), sep=""))
            }
            idx <- idx[!is.na(idx)]
        }
    } else {
        # ids were null, just return an index vector
        # allong all_ids
        idx <- seq_along(all_ids)
    }
    # subset the character ids to the ones we want
    id_keep <- as.character(all_ids[idx])
    return(list(idx=idx, ids=id_keep))
}

# define the initialization method for the GCT class
setMethod("initialize", signature = "GCT", definition = function(
    .Object, mat=NULL, rdesc=NULL, cdesc=NULL, src=NULL, rid=NULL, cid=NULL,
    set_annot_rownames=F, matrix_only=F) {
    # if we were supplied a matrix and annotations, use them
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
        # we were not given a matrix, were we given a src file?
        # check to make sure it's either .gct or .gctx
        if (! (grepl(".gct$", src) || grepl(".gctx$", src) ))
            stop("Either a .gct or .gctx file must be given")
        if (grepl(".gct$", src)) {
            if ( ! is.null(rid) || !is.null(cid) )
                warning(paste(
                    "rid and cid values may only be given for .gctx files, not",
                    ".gct files\nignoring"))
            # parse the .gct
            .Object@src = src
            # get the .gct version by reading first line
            .Object@version = scan(src, what = "", nlines = 1, sep = "\t",
                                   quiet = TRUE)[1]
            # get matrix dimensions by reading second line
            dimensions = scan(src, what = double(0), nlines = 1, skip = 1,
                              sep = "\t", quiet = TRUE)
            nrmat = dimensions[1]
            ncmat = dimensions[2]
            if (length(dimensions)==4) {
                # a #1.3 file
                message("parsing as GCT v1.3")
                nrhd <- dimensions[3]
                nchd <- dimensions[4]
            } else {
                # a #1.2 file
                message("parsing as GCT v1.2")
                nrhd <- 0
                nchd <- 0
            }
            message(paste(src, nrmat, "rows,", ncmat, "cols,", nrhd,
                          "row descriptors,", nchd, "col descriptors"))
            # read in header line
            header = scan(src, what = "", nlines = 1, skip = 2, sep = "\t",
                          quote = NULL, quiet = TRUE)
            # construct row header and column id's from the header line
            if ( nrhd > 0 ) {
                rhd <- header[2:(nrhd+1)]
                cid <- header[-(nrhd+1):-1]
                col_offset <- 1
            }
            else {
                if (any(grepl("description", header, ignore.case=T))) {
                    # check for presence of description column in v1.2 files
                    col_offset <- 2
                } else {
                    col_offset <- col_offset <- 1
                }
                rhd = NULL
                cid = header[(1+col_offset):length(header)]
            }
            # read in the next set of headers (column annotations) and shape
            # into a matrix
            if ( nchd > 0 ) {
                header = scan(src, what = "", nlines = nchd, skip = 3,
                              sep = "\t", quote = NULL, quiet = TRUE)
                header = matrix(header, nrow = nchd,
                                ncol = ncmat + nrhd + 1, byrow = TRUE)
                # extract the column header and column descriptions
                chd = header[,1]
                cdesc = header[,-(nrhd+1):-1]
                # need to transpose in the case where there's only one column
                # annotation
                if ( nchd == 1 ) cdesc = t(cdesc)
            } else {
                chd = NULL
                cdesc = data.frame()
            }
            # read in the data matrix and row descriptions, shape into a matrix
            mat = scan(src, what = "", nlines = nrmat,
                       skip = 3 + nchd, sep = "\t", quote = NULL, quiet = TRUE)
            mat = matrix(mat, nrow = nrmat, ncol = ncmat + nrhd + col_offset,
                         byrow = TRUE)
            # message(paste(dim(mat), collapse="\t"))
            # Extract the row id's row descriptions, and the data matrix
            rid = mat[,1]
            if ( nrhd > 0 ) {
                # need as.matrix for the case where there's only one row
                # annotation
                rdesc = as.matrix(mat[,2:(nrhd + 1)])
                mat = matrix(as.numeric(mat[,-(nrhd + 1):-1]),
                             nrow = nrmat, ncol = ncmat)
            }
            else {
                rdesc = data.frame()
                mat = matrix(as.numeric(mat[, (1+col_offset):ncol(mat)]),
                             nrow = nrmat, ncol = ncmat)
            }
            # assign names to the data matrix and the row and column
            # descriptions
            # message(paste(dim(mat), collapse="\t"))
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
            # assign to the GCT slots
            .Object@mat = mat
            .Object@rid = rownames(mat)
            .Object@cid = colnames(mat)
            if (!matrix_only) {
                # return annotations as well as matrix
                .Object@rdesc = fix.datatypes(rdesc)
                .Object@cdesc = fix.datatypes(cdesc)
                # add id columns to rdesc and cdesc
                .Object@rdesc$id <- rownames(.Object@rdesc)
                .Object@cdesc$id <- rownames(.Object@cdesc)
            }
        } else {
            # parse the .gctx
            message(paste("reading", src))
            .Object@src = src
            # if the rid's or column id's are .grp files, read them in
            if ( length(rid) == 1 && grepl(".grp$", rid) )
                rid <- parse.grp(rid)
            if ( length(cid) == 1 && grepl(".grp$", cid) )
                cid <- parse.grp(cid)
            # get all the row and column ids
            all_rid <- read.gctx.ids(src, dimension="row")
            all_cid <- read.gctx.ids(src, dimension="col")
            # if rid or cid specified, read only those rows/columns
            # if already numeric, use as is
            # else convert to numeric indices
            processed_rids <- process_ids(rid, all_rid, type="rid")
            processed_cids <- process_ids(cid, all_cid, type="cid")
            # read the data matrix
            .Object@mat <- h5read(src, name="0/DATA/0/matrix",
                                  index=list(processed_rids$idx,
                                             processed_cids$idx))
            # set the row and column ids, casting as characters
            .Object@rid <- processed_rids$ids
            .Object@cid <- processed_cids$ids
            rownames(.Object@mat) <- processed_rids$ids
            colnames(.Object@mat) <- processed_cids$ids
            # get the meta data
            if (!matrix_only) {
                .Object@rdesc <- read.gctx.meta(
                    src, dimension="row", ids=processed_rids$ids,
                    set_annot_rownames=set_annot_rownames)
                .Object@cdesc <- read.gctx.meta(
                    src, dimension="col", ids=processed_cids$ids,
                    set_annot_rownames=set_annot_rownames)
            } else {
                .Object@rdesc <- data.frame(id=.Object@rid,
                                            stringsAsFactors = F)
                .Object@cdesc <- data.frame(id=.Object@cid,
                                            stringsAsFactors = F)
            }
            # close any open handles and return the object
            closeOpenHandles()
            message("done")
        }
    }
    # finally, make sure object is valid before returning
    ok <- validObject(.Object)
    return(.Object)
})

#' Close open handles
#'
#' @importFrom utils packageVersion
#'
#' @return Closes all open identifiers
closeOpenHandles <- function() {
    if(packageVersion('rhdf5') < "2.23.0")
        rhdf5::H5close()
    else
        rhdf5::h5closeAll()
}

#' Parse a GCTX file into the workspace as a GCT object
#'
#' @param fname path to the GCTX file on disk
#' @param cid either a vector of character or integer
#'   column indices or a path to a grp file containing character
#'   column indices. Only these indicies will be parsed from the
#'   file.
#'
#' @importFrom methods new
#'
#' @details \code{parse.gctx} also supports parsing of plain text
#'   GCT files, so this function can be used as a general GCT parser.
#'
#' @examples
#' gct_file <- system.file("extdata", "modzs_n272x978.gctx", package="cmapR")
#' (ds <- parse.gctx(gct_file))
#'
#' # matrix only
#' (ds <- parse.gctx(gct_file, matrix_only=TRUE))
#'
#' # only the first 10 rows and columns
#' (ds <- parse.gctx(gct_file, rid=1:10, cid=1:10))
#'
#' @family GCTX parsing functions
#' @export
parse.gctx <- function(fname, cid=NULL)
    new("GCT", src = fname, cid = cid)

# utils.R ----------------------------------------------------------------------

#' Check whether \code{test_names} are columns in the \code{\link{data.frame}}
#' @param test_names a vector of column names to test
#' @param df the \code{\link{data.frame}} to test against
#' @param throw_error boolean indicating whether to throw an error if
#'   any \code{test_names} are not found in \code{df}
#' @return boolean indicating whether or not all \code{test_names} are
#'   columns of \code{df}
#'
#' @source https://github.com/cmap/cmapR
#'
#' @examples
#' check_colnames(c("pert_id", "pert_iname"), cdesc_char)            # TRUE
#' check_colnames(c("pert_id", "foobar"), cdesc_char, throw_error=FALSE) # FALSE
#' @export
check_colnames <- function(test_names, df, throw_error=T) {
    # check whether test_names are valid names in df
    # throw error if specified
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
#' @source https://github.com/cmap/cmapR
#'
#' @keywords internal
subset_to_ids <- function(df, ids) {
    # helper function to do a robust df subset
    check_colnames("id", df)
    newdf <- data.frame(df[match(ids, df$id), ])
    names(newdf) <- names(df)
    return(newdf)
}
