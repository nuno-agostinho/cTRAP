# Based on cmapR package at https://github.com/cmap/cmapR

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
#' @seealso \code{\link{parse.gctx}}, \code{\link{write.gctx}}, \code{\link{read.gctx.meta}}, \code{\link{read.gctx.ids}}
#' @seealso \link{http://clue.io/help} for more information on the GCT format
methods::setClass("GCT", methods::representation(
    mat = "matrix", rid = "character", cid = "character", rdesc = "data.frame",
    cdesc = "data.frame", version = "character", src = "character"))

# set up methods for checking GCT validity
methods::setValidity("GCT", function(object) {
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
        return("cdesc must either have 0 rows or the same number of rows as matrix has columns")
    }
    if (nrow(object@rdesc) != nrows & nrow(object@rdesc) != 0) {
        return("rdesc must either have 0 rows or the same number of rows as matrix has rows")
    } else {
        return(TRUE)
    }
})

suppressMessages({
    # set method for displaying a GCT object
    # just use the 'str' function to show its structure
    setMethod("show", methods::signature("GCT"), function(object) {
        utils::str(object)
    })

    # dim, nrow and ncol to display the # of rows and columns
    # for a GCT object's matrix
    setMethod("ncol", methods::signature("GCT"), function(x) ncol(x@mat))
    setMethod("nrow", methods::signature("GCT"), function(x) nrow(x@mat))
    setMethod("dim", methods::signature("GCT"), function(x) dim(x@mat))
    setMethod("range", methods::signature("GCT"), function(x, na.rm=F,
                                                           finite=F) {
        range(x@mat, na.rm=na.rm, finite=finite)
    })
    setMethod("max", methods::signature("GCT"), function(x, na.rm=F) {
        max(x@mat, na.rm=na.rm)
    })
    setMethod("min", methods::signature("GCT"), function(x, na.rm=F) {
        min(x@mat, na.rm=na.rm)
    })
    setMethod("diag", methods::signature("GCT"), function(x) diag(x@mat))
})


#### define some helper methods for parsing gctx files ###

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
#' col_meta_first10 <- read.gctx.meta(gct_file, dimension="column", ids=col_meta$id[1:10])
#' str(col_meta_first10)
#'
#' @family GCTX parsing functions
#' @export
read.gctx.meta <- function(gctx_path, dimension="row", ids=NULL, set_annot_rownames=T) {
    if (!file.exists(gctx_path)) {
        stop(paste(gctx_path, "does not exist"))
    }
    if (dimension=="column") dimension <- "col"
    if (!(dimension %in% c("row", "col"))) {
        stop("dimension can be either row or col")
    }
    if (dimension == "row") {
        name <- "0/META/ROW"
    } else {
        name <- "0/META/COL"
    }
    raw_annots <- rhdf5::h5read(gctx_path, name=name) # returns a list
    fields <- names(raw_annots)
    # define an empty data frame of the correct dimensions
    annots <-  data.frame(matrix(nrow=length(raw_annots[[fields[1]]]), ncol=length(fields)))
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
    if (is.null(ids)) {
        ids <- as.character(annots$id)
    } else {
        ids <- ids
    }
    # make sure annots row ordering matches that of ids
    annots <- subset_to_ids(annots, ids)
    annots$id <- as.character(annots$id)
    # use the id field to set the rownames
    if (set_annot_rownames) {
        rownames(annots) <- annots$id
    }
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
    ids <- gsub("\\s*$", "", rhdf5::h5read(gctx_path, name=name), perl=T)
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
                stop(paste("none of the requested", type, "indices were found in the dataset"))
            }
            if (any(is_invalid_idx)) {
                # requested indices are outside of the possible range
                warning(paste("the following ", type, " were are outside possible range and will be ignored:\n",
                              paste(invalid_idx, collapse="\n"), sep=""))
            }
            idx <- idx[!is_invalid_idx]
        } else {
            # assume its a character
            idx <- match(ids, all_ids)
            if (all(is.na(idx))) {
                stop(paste("none of the requested", type, "were found in the dataset"))
            }
            if (any(is.na(idx))) {
                ids_not_found <- ids[is.na(idx)]
                warning(paste("the following ", type, " were not found and will be ignored:\n",
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
methods::setMethod("initialize",
                   signature = "GCT",
                   definition = function(.Object, mat=NULL, rdesc=NULL, cdesc=NULL, src=NULL, rid=NULL, cid=NULL, set_annot_rownames=F,
                                         matrix_only=F) {
                       # if we were supplied a matrix and annotations, use them
                       if (!is.null(mat)) {
                           .Object@mat <- mat
                           # if given rid and cid, use those as well
                           if (!is.null(rid)) {
                               .Object@rid <- rid
                           } else {
                               .Object@rid <- rownames(mat)
                           }
                           if (!is.null(cid)) {
                               .Object@cid <- cid
                           } else {
                               .Object@cid <- colnames(mat)
                           }
                       }
                       if (!is.null(rdesc)) {
                           .Object@rdesc <- rdesc
                       }
                       if (!is.null(cdesc)) {
                           .Object@cdesc <- cdesc
                       } else if (!is.null(src)) {
                           # we were not given a matrix, were we given a src file?
                           # check to make sure it's either .gct or .gctx
                           if (! (grepl(".gct$", src) || grepl(".gctx$", src) ))
                               stop("Either a .gct or .gctx file must be given")
                           if (grepl(".gct$", src)) {
                               if ( ! is.null(rid) || !is.null(cid) )
                                   warning(paste("rid and cid values may only be given for .gctx files, not .gct files\n",
                                                 "ignoring"))
                               # parse the .gct
                               .Object@src = src
                               # get the .gct version by reading first line
                               .Object@version = scan(src, what = "", nlines = 1, sep = "\t", quiet = TRUE)[1]
                               # get matrix dimensions by reading second line
                               dimensions = scan(src, what = double(0), nlines = 1, skip = 1, sep = "\t", quiet = TRUE)
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
                               message(paste(src, nrmat, "rows,", ncmat, "cols,", nrhd, "row descriptors,", nchd, "col descriptors"))
                               # read in header line
                               header = scan(src, what = "", nlines = 1, skip = 2, sep = "\t", quote = NULL, quiet = TRUE)
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
                               # read in the next set of headers (column annotations) and shape into a matrix
                               if ( nchd > 0 ) {
                                   header = scan(src, what = "", nlines = nchd, skip = 3, sep = "\t",
                                                 quote = NULL, quiet = TRUE)
                                   header = matrix(header, nrow = nchd,
                                                   ncol = ncmat + nrhd + 1, byrow = TRUE)
                                   # extract the column header and column descriptions
                                   chd = header[,1]
                                   cdesc = header[,-(nrhd+1):-1]
                                   # need to transpose in the case where there's only one column annotation
                                   if ( nchd == 1 )
                                       cdesc = t(cdesc)
                               }
                               else {
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
                                   # need as.matrix for the case where there's only one row annotation
                                   rdesc = as.matrix(mat[,2:(nrhd + 1)])
                                   mat = matrix(as.numeric(mat[,-(nrhd + 1):-1]),
                                                nrow = nrmat, ncol = ncmat)
                               }
                               else {
                                   rdesc = data.frame()
                                   mat = matrix(as.numeric(mat[, (1+col_offset):ncol(mat)]), nrow = nrmat, ncol = ncmat)
                               }
                               # assign names to the data matrix and the row and column descriptions
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
                           }
                           else {
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
                               .Object@mat <- rhdf5::h5read(src, name="0/DATA/0/matrix",
                                                            index=list(processed_rids$idx, processed_cids$idx))
                               # set the row and column ids, casting as characters
                               .Object@rid <- processed_rids$ids
                               .Object@cid <- processed_cids$ids
                               rownames(.Object@mat) <- processed_rids$ids
                               colnames(.Object@mat) <- processed_cids$ids
                               # get the meta data
                               if (!matrix_only) {
                                   .Object@rdesc <- read.gctx.meta(src, dimension="row", ids=processed_rids$ids,
                                                                   set_annot_rownames=set_annot_rownames)
                                   .Object@cdesc <- read.gctx.meta(src, dimension="col", ids=processed_cids$ids,
                                                                   set_annot_rownames=set_annot_rownames)
                               }
                               else {
                                   .Object@rdesc <- data.frame(id=.Object@rid, stringsAsFactors = F)
                                   .Object@cdesc <- data.frame(id=.Object@cid, stringsAsFactors = F)
                               }
                               # close any open handles and return the object
                               if(utils::packageVersion('rhdf5') < "2.23.0") {
                                   rhdf5::H5close()
                               } else {
                                   rhdf5::h5closeAll()
                               }
                               message("done")
                           }
                       }
                       # finally, make sure object is valid before returning
                       ok <- methods::validObject(.Object)
                       return(.Object)
                   }
)


#' Parse a GCTX file into the workspace as a GCT object
#'
#' @param fname path to the GCTX file on disk
#' @param rid either a vector of character or integer
#'   row indices or a path to a grp file containing character
#'   row indices. Only these indicies will be parsed from the
#'   file.
#' @param cid either a vector of character or integer
#'   column indices or a path to a grp file containing character
#'   column indices. Only these indicies will be parsed from the
#'   file.
#' @param set_annot_rownames boolean indicating whether to set the
#'   rownames on the row/column metadata data.frames. Set this to
#'   false if the GCTX file has duplicate row/column ids.
#' @param matrix_only boolean indicating whether to parse only
#'   the matrix (ignoring row and column annotations)
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
parse.gctx <- function(fname, rid=NULL, cid=NULL, set_annot_rownames=F, matrix_only=F) {
    ds <- methods::new("GCT",
                       src = fname,
                       rid = rid,
                       cid = cid,
                       set_annot_rownames = set_annot_rownames,
                       matrix_only = matrix_only)
    return(ds)
}

#' Append matrix dimensions to filename
#'
#' @param ofile the file name
#' @param mat the matrix
#' @param extension the file extension
#'
#' @return a character string of the filename with
#'   matrix dimensions appended
#'
#' @details This is a helper function that most users
#'   will not use directly
#'
#' @examples
#' (filename <- cmapR:::append.dim("my.gctx.filename", matrix(nrow=10, ncol=15)))
#'
#' @keywords internal
#' @family GCTX parsing functions
append.dim <- function(ofile, mat, extension="gct") {
    nc <- ncol(mat)
    nr <- nrow(mat)
    filename <- basename(ofile)
    if (grepl("n[0-9]+x[0-9]+\\.gct", filename)) {
        # already has a dimensions token, ignore
        filename <- sub("_n[0-9]+x[0-9]+\\.gct.*", "", filename)
    }
    filename <- file.path(dirname(ofile),
                          sprintf('%s_n%dx%d.%s',filename,
                                  nc, nr, extension))
    return(filename)
}


#' Write a GCT object to disk in GCT format
#'
#' @param ds the GCT object
#' @param ofile the desired output filename
#' @param precision the numeric precision at which to
#'   save the matrix. See \code{details}.
#' @param appenddim boolean indicating whether to append
#'   matrix dimensions to filename
#' @param ver the GCT version to write. See \code{details}.
#'
#' @details Since GCT is text format, the higher \code{precision}
#'   you choose, the larger the file size.
#'   \code{ver} is assumed to be 3, aka GCT version 1.3, which supports
#'   embedded row and column metadata in the GCT file. Any other value
#'   passed to \code{ver} will result in a GCT version 1.2 file which
#'   contains only the matrix data and no annotations.
#'
#' @return NULL
#'
#' @examples
#' \dontrun{
#' write.gct(ds, "dataset", precision=2)
#' }
#' @family GCTX parsing functions
#' @export
write.gct <- function(ds, ofile, precision=4, appenddim=T, ver=3) {
    if (!class(ds)=="GCT") {
        stop("ds must be a GCT object")
    }
    # make sure it's valid
    ok <- methods::validObject(ds)
    # append the dimensions of the data set, if desired
    if (appenddim) ofile <- append.dim(ofile, ds@mat, extension="gct")

    precision = floor(precision)
    cat(sprintf('Saving file to %s\n',ofile))
    nr <- nrow(ds@mat)
    nc <- ncol(ds@mat)
    cat(sprintf('Dimensions of matrix: [%dx%d]\n',nr,nc))
    cat(sprintf('Setting precision to %d\n',precision))

    # open file and write
    if (ver==3) {
        # remove the 'id' columns
        ds@cdesc$id <- NULL
        ds@rdesc$id <- NULL
        # get the counts of meta data fields
        nrdesc = dim(ds@rdesc)[2]
        ncdesc = dim(ds@cdesc)[2]
        colkeys = colnames(ds@cdesc)
        # append header
        cat(sprintf('#1.%d\n%d\t%d\t%d\t%d', ver, nr, nc, nrdesc, ncdesc),
            file=ofile,sep='\n')
        # line 3: sample row desc keys and sample names
        cat(paste(c('id',colnames(ds@rdesc),ds@cid),collapse='\t'),
            file=ofile,sep='\n',append=T)
        # line 4 + ncdesc: sample desc
        filler = 'na'
        if (ncdesc > 0) {
            for (ii in 1:ncdesc) {
                if (is.numeric(ds@cdesc[,ii])) {
                    cat(paste(c(colkeys[ii],rep(filler,nrdesc),
                                round(ds@cdesc[,ii],precision)),
                              collapse='\t'),
                        file=ofile,sep='\n',append=T)
                } else {
                    cat(paste(c(colkeys[ii],rep(filler,nrdesc),
                                ds@cdesc[,ii]),
                              collapse='\t'),
                        file=ofile,sep='\n',append=T)
                }
            }
        }

        for (ii in 1:nr) {
            # print rows
            cat(paste(c(ds@rid[ii],
                        ds@rdesc[ii,],
                        round(ds@mat[ii,],precision)),collapse='\t'),
                sep='\n',file=ofile,append=T)
        }
    } else {
        # assume ver 1.2 and below, ignore descriptors
        # append header
        cat(sprintf('#1.%d\n%d\t%d', ver, nr, nc),
            file=ofile,sep='\n')
        # line 3: sample row desc keys and sample names
        cat(paste(c('id','Description',ds@cid),collapse='\t'),
            file=ofile,sep='\n',append=T)

        for (ii in 1:nr) {
            # print rows
            cat(paste(c(ds@rid[ii],
                        ds@rdesc[ii, 2],
                        round(ds@mat[ii,],precision)),collapse='\t'),
                sep='\n',file=ofile,append=T)
        }
    }

    cat(sprintf('Saved.\n'))
}


#' Write a GCT object to disk in GCTX format
#'
#' @param ds a GCT object
#' @param ofile the desired file path for writing
#' @param appenddim boolean indicating whether the
#'   resulting filename will have dimensions appended
#'   (e.g. my_file_n384x978.gctx)
#' @param compression_level integer between 1-9 indicating
#'   how much to compress data before writing. Higher values
#'   result in smaller files but slower read times.
#' @param matrix_only boolean indicating whether to write
#'   only the matrix data (and skip row, column annotations)
#' @param max_chunk_kb for chunking, the maximum number of KB
#'   a given chunk will occupy
#'
#' @examples
#' \dontrun{
#' # assume ds is a GCT object
#' write.gctx(ds, "my/desired/outpath/and/filename")
#' }
#' @family GCTX parsing functions
#' @export
write.gctx <- function(ds, ofile, appenddim=T, compression_level=0, matrix_only=F,
                       max_chunk_kb=1024) {
    if (!class(ds)=="GCT") {
        stop("ds must be a GCT object")
    }
    # make sure it's valid
    ok <- methods::validObject(ds)
    # add dimensions to filename if desired
    if (appenddim) ofile <- append.dim(ofile, ds@mat, extension="gctx")
    # check if the file already exists
    if (file.exists(ofile)) {
        message(paste(ofile, "exists, removing"))
        file.remove(ofile)
    }
    message(paste("writing", ofile))

    # start the file object
    rhdf5::h5createFile(ofile)

    # create all the necessary groups
    rhdf5::h5createGroup(ofile, "0")
    rhdf5::h5createGroup(ofile, "0/DATA")
    rhdf5::h5createGroup(ofile, "0/DATA/0")
    rhdf5::h5createGroup(ofile, "0/META")
    rhdf5::h5createGroup(ofile, "0/META/COL")
    rhdf5::h5createGroup(ofile, "0/META/ROW")

    # create and write matrix data, using chunking
    bits_per_element <- switch(storage.mode(ds@mat),
                               "double" = 64,
                               "integer" = 32)
    elem_per_kb <- max_chunk_kb * 8 / bits_per_element
    # assume matrix is of dimensions row_dim x col_dim
    row_dim <- nrow(ds)
    col_dim <- ncol(ds)
    row_chunk_size <- min(row_dim, 1000)
    # column chunk, such that row * col <= max_chunk_kb
    col_chunk_size <- min(((max_chunk_kb * elem_per_kb) %/% row_chunk_size), col_dim)
    chunking <- c(row_chunk_size, col_chunk_size)
    message(paste(c("chunk sizes:", chunking), collapse="\t"))
    rhdf5::h5createDataset(ofile, "0/DATA/0/matrix", dim(ds@mat), chunk=chunking, level=compression_level)
    rhdf5::h5write.default(ds@mat, ofile, "0/DATA/0/matrix")

    # write annotations
    rhdf5::h5write.default(as.character(ds@rid), ofile, "0/META/ROW/id")
    rhdf5::h5write.default(as.character(ds@cid), ofile, "0/META/COL/id")

    if (!matrix_only) {
        write.gctx.meta(ofile, ds@cdesc, dimension="column")
        write.gctx.meta(ofile, ds@rdesc, dimension="row")
    }

    # close any open handles
    if(utils::packageVersion('rhdf5') < "2.23.0") {
        rhdf5::H5close()
    } else {
        rhdf5::h5closeAll()
    }

    # add the version annotation and close
    fid <- rhdf5::H5Fopen(ofile)
    rhdf5::h5writeAttribute.character("GCTX1.0", fid, "version")
    rhdf5::H5Fclose(fid)

}

#' Write a \code{data.frame} of meta data to GCTX file
#'
#' @param ofile the desired file path for writing
#' @param df the \code{data.frame} of annotations
#' @param dimension the dimension to annotate
#'   (row or column)
#'
#' @examples
#' \dontrun{
#' # assume ds is a GCT object
#' cmapR:::write.gctx.meta("/my/file/path", cdesc_char, dimension="col")
#' }
#' @family GCTX parsing functions
#' @keywords internal
write.gctx.meta <- function(ofile, df, dimension="row") {
    path <- if ((dimension=="row")) "0/META/ROW/" else "0/META/COL/"
    # loop through all columns
    fields <- names(df)
    if (length(fields) > 0) {
        for (i in 1:length(fields)) {
            field <- fields[i]
            # if this is the id field, skip b/c that field is special
            # and is written as part of write.gctx
            if (field == "id") next
            v <- df[, i]
            # convert factors to character
            if(class(v) == "factor" || class(v) == "AsIs") {
                v <- as.character(v)
            }
            rhdf5::h5write.default(v, ofile, paste(path, field, sep=""))
        }
    }
}


###########################################
### functions for other CMap file types ###
###########################################

#' Read a GRP file and return a vector of its contents
#' @param fname the file path to be parsed
#' @return a vector of the contents of \code{fname}
#' @examples
#' grp_path <- system.file("extdata", "lm_epsilon_n978.grp", package="cmapR")
#' values <- parse.grp(grp_path)
#' str(values)
#' @family CMap parsing functions
#' @seealso \link{http://clue.io/help} for details on the GRP file format
#' @export
parse.grp <- function(fname) {
    grp <- scan(fname, what = "", quote = NULL, quiet = TRUE, sep="\n")
    return(grp)
}


#' Write a vector to a GRP file
#'
#' @param vals the vector of values to be written
#' @param fname the desired file name
#'
#' @examples
#' \dontrun{
#' write.grp(letters, "letter.grp")
#' }
#'
#' @family CMap parsing functions
#' @seealso \link{http://clue.io/help} for details on the GRP file format
#' @export
write.grp <- function(vals, fname) {
    if (is.list(vals)) vals <- unlist(vals)
    if (!is.vector(vals)) vals <- as.vector(vals)
    write(vals, fname, ncolumns=1)
}


#' Read a GMX file and return a list
#'
#' @param fname the file path to be parsed
#'
#' @return a list of the contents of \code{fname}. See details.
#'
#' @details \code{parse.gmx} returns a nested list object. The top
#'   level contains one list per column in \code{fname}. Each of
#'   these is itself a list with the following fields:
#'   - \code{head}: the name of the data (column in \code{fname})
#'   - \code{desc}: description of the corresponding data
#'   - \code{len}: the number of data items
#'   - \code{entry}: a vector of the data items
#'
#' @examples
#' gmx_path <- system.file("extdata", "lm_probes.gmx", package="cmapR")
#' gmx <- parse.gmx(gmx_path)
#' str(gmx)
#'
#' @family CMap parsing functions
#' @seealso \link{http://clue.io/help} for details on the GMX file format
#' @export
parse.gmx <- function(fname) {
    tmp <- utils::read.table(fname, sep = "\t",
                             header = TRUE, stringsAsFactors = FALSE)
    # preallocate a list for the gmx
    L <- list()
    # loop over the first row of the .gmx
    for ( n in names(tmp) ) {
        # get all the values; remove empties at the end
        values <- tmp[[n]][-1]
        remove.idx <- values == ""
        values <- values[!remove.idx]
        # put in a list
        L[[n]] <- list(head = n,
                       desc = tmp[[n]][1],
                       len = length(values),
                       entry = values)
    }
    return(L)
}


#' Read a GMT file and return a list
#' @param fname the file path to be parsed
#'
#' @return a list of the contents of \code{fname}. See details.
#'
#' @details \code{parse.gmt} returns a nested list object. The top
#'   level contains one list per row in \code{fname}. Each of
#'   these is itself a list with the following fields:
#'   - \code{head}: the name of the data (row in \code{fname})
#'   - \code{desc}: description of the corresponding data
#'   - \code{len}: the number of data items
#'   - \code{entry}: a vector of the data items
#'
#' @examples
#' gmt_path <- system.file("extdata", "query_up.gmt", package="cmapR")
#' gmt <- parse.gmt(gmt_path)
#' str(gmt)
#'
#' @family CMap parsing functions
#' @seealso \link{http://clue.io/help} for details on the GMT file format
#' @export
parse.gmt <- function(fname) {
    gmt.lines <- scan(fname, what = "", sep = "\n",
                      quote = NULL, quiet = TRUE)
    tmp <- lapply(gmt.lines, function(x) unlist(strsplit(x, "\t")))
    mk.gmt.entry <- function(x) {
        L <- list()
        L[["head"]] <- x[1]
        L[["desc"]] <- x[2]
        l.entry <- x[-c(1:2)]
        idx <- l.entry != ""
        L[["entry"]] <- l.entry[idx]
        L[["len"]] <- length(L[["entry"]])
        return(L)
    }
    L <- lapply(tmp, function(x) mk.gmt.entry(x))
    names(L) <- unlist(lapply(L, function(x) x$head))
    return(L)
}


#' Write a nested list to a GMT file
#'
#' @param lst the nested list to write. See \code{details}.
#' @param fname the desired file name
#'
#' @details \code{lst} needs to be a nested list where each
#'   sub-list is itself a list with the following fields:
#'   - \code{head}: the name of the data
#'   - \code{desc}: description of the corresponding data
#'   - \code{len}: the number of data items
#'   - \code{entry}: a vector of the data items
#'
#' @examples
#' \dontrun{
#' write.gmt(gene_set, "gene_set.gmt")
#' }
#'
#' @family CMap parsing functions
#' @seealso \link{http://clue.io/help} for details on the GMT file format
#' @export
write.gmt <- function(lst, fname) {
    # assumes that each element of the list will have the fields
    # head, desc, entry
    if (file.exists(fname)) {
        message(paste(fname, "exists, deleting..."))
        file.remove(fname)
    }
    for (i in 1:length(lst)) {
        el <- lst[[i]]
        ncolumns <- 2 + length(el$entry)
        write(c(el$head, el$desc, el$entry), file=fname, sep="\t", append=T, ncolumns=ncolumns)
    }
}


########################################
### Other Misc. utility functions ######
########################################


#' Write a \code{data.frame} to a tab-delimited text file
#'
#' @param tbl the \code{data.frame} to be written
#' @param ofile the desired file name
#' @param ... additional arguments passed on to \code{write.table}
#'
#' @details This method simply calls \code{write.table} with some
#'   preset arguments that generate a unquoated, tab-delimited file
#'   without row names.
#'
#' @examples
#' \dontrun{
#' write.tbl(cdesc_char, "col_meta.txt")
#' }
#'
#' @seealso \code{\link{write.table}}
#' @export
write.tbl <- function(tbl, ofile, ...) {
    utils::write.table(tbl, file = ofile, sep="\t", quote=F,
                       col.names=T, row.names=F, ...)
}

# utils.R ----------------------------------------------------------------------

#' Transform a GCT object in to a long form \code{\link{data.table}} (aka 'melt')
#'
#' @description Utilizes the \code{\link{data.table::melt}} function to transform the
#'   matrix into long form. Optionally can include the row and column
#'   annotations in the transformed \code{\link{data.table}}.
#'
#' @param g the GCT object
#' @param keep_rdesc boolean indicating whether to keep the row
#'   descriptors in the final result
#' @param keep_cdesc boolean indicating whether to keep the column
#'   descriptors in the final result
#' @param remove_symmetries boolean indicating whether to remove
#'   the lower triangle of the matrix (only applies if \code{g@mat} is symmetric)
#' @param suffixes the character suffixes to be applied if there are
#'   collisions between the names of the row and column descriptors
#' @param ... further arguments passed along to \code{data.table::merge}
#'
#' @return a \code{\link{data.table}} object with the row and column ids and the matrix
#'   values and (optinally) the row and column descriptors
#'
#' @examples
#' # simple melt, keeping both row and column meta
#' head(melt.gct(ds))
#'
#' # update row/colum suffixes to indicate rows are genes, columns experiments
#' head(melt.gct(ds, suffixes = c("_gene", "_experiment")))
#'
#' # ignore row/column meta
#' head(melt.gct(ds, keep_rdesc = FALSE, keep_cdesc = FALSE))
#'
#' @family GCT utilities
#' @export
setGeneric("melt.gct", function(g, suffixes=NULL, remove_symmetries=F,
                                keep_rdesc=T, keep_cdesc=T, ...) {
    standardGeneric("melt.gct")
})
setMethod("melt.gct", signature("GCT"),
          function(g, suffixes, remove_symmetries=F, keep_rdesc=T, keep_cdesc=T, ...) {
              # melt a gct object's matrix into a data.frame and merge row and column
              # annotations back in, using the provided suffixes
              # assumes rdesc and cdesc data.frames both have an 'id' field.
              # merges row and/or column annotations into the melted matrix as indicated by
              # keep_rdesc and keep_cdesc, respectively.
              # if remove_symmetries, will check whether matrix is symmetric
              # and return only values corresponding to the upper triangle
              # g@rdesc$id <- rownames(g@rdesc)
              # g@cdesc$id <- rownames(g@cdesc)
              # first, check if matrix is symmetric
              # if it is, use only the upper triangle
              message("melting GCT object...")
              mat <- g@mat
              if (remove_symmetries & isSymmetric(mat)) {
                  mat[upper.tri(mat, diag=F)] <- NA
              }
              mat <- data.table::data.table(mat)
              mat$rid <- g@rid
              d <- data.table::melt(mat, id.vars="rid")
              data.table::setattr(d, "names", c("id.x", "id.y", "value"))
              d$id.x <- as.character(d$id.x)
              d$id.y <- as.character(d$id.y)
              # standard data.frame subset here to comply with testthat
              d <- subset(d, !is.na(value))
              if (keep_rdesc && keep_cdesc) {
                  # merge back in both row and column descriptors
                  data.table::setattr(d, "names", c("id", "id.y", "value"))
                  d <- merge(d, data.table::data.table(g@rdesc), by="id", ...)
                  data.table::setnames(d, "id", "id.x")
                  data.table::setnames(d, "id.y", "id")
                  d <- merge(d, data.table::data.table(g@cdesc), by="id", ...)
                  data.table::setnames(d, "id", "id.y")
              } else if (keep_rdesc) {
                  # keep only row descriptors
                  rdesc <- data.table::data.table(g@rdesc)
                  data.table::setnames(rdesc, "id", "id.x")
                  d <- merge(d, rdesc, by="id.x", ...)
              } else if (keep_cdesc) {
                  # keep only column descriptors
                  cdesc <- data.table::data.table(g@cdesc)
                  data.table::setnames(cdesc, "id", "id.y")
                  d <- merge(d, cdesc, by="id.y", ...)
              }
              # use suffixes if provided
              if (!is.null(suffixes) & length(suffixes) == 2) {
                  newnames <- gsub("\\.x", suffixes[1], names(d))
                  newnames <- gsub("\\.y", suffixes[2], newnames)
                  data.table::setattr(d, "names", newnames)
              }
              message("done")
              return(d)
          })


#' Check if x is a whole number
#'
#' @param x number to test
#' @param tol the allowed tolerance
#' @return boolean indicating whether x is tol away from a whole number value
#' @examples
#' is.wholenumber(1)
#' is.wholenumber(0.5)
#' @export
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {
    return(abs(x - round(x)) < tol)
}

#' Check whether \code{test_names} are columns in the \code{\link{data.frame}} df
#' @param test_names a vector of column names to test
#' @param df the \code{\link{data.frame}} to test against
#' @param throw_error boolean indicating whether to throw an error if
#'   any \code{test_names} are not found in \code{df}
#' @return boolean indicating whether or not all \code{test_names} are
#'   columns of \code{df}
#' @examples
#' check_colnames(c("pert_id", "pert_iname"), cdesc_char)            # TRUE
#' check_colnames(c("pert_id", "foobar"), cdesc_char, throw_error=FALSE) # FALSE, suppress error
#' @export
check_colnames <- function(test_names, df, throw_error=T) {
    # check whether test_names are valid names in df
    # throw error if specified
    diffs <- setdiff(test_names, names(df))
    if (length(diffs) > 0) {
        if (throw_error) {
            stop(paste("the following column names are not found in", deparse(substitute(df)), ":",
                       paste(diffs, collapse=" "), "\n"))
        } else {
            return(F)
        }
    } else {
        return(T)
    }
}

#' Do a robust \code{\link{data.frame}} subset to a set of ids
#' @param df \code{\link{data.frame}} to subset
#' @param ids the ids to subset to
#' @return a subset version of \code{df}
#' @keywords internal
subset_to_ids <- function(df, ids) {
    # helper function to do a robust df subset
    check_colnames("id", df)
    newdf <- data.frame(df[match(ids, df$id), ])
    names(newdf) <- names(df)
    return(newdf)
}


#' Subset a gct object using the provided row and column ids
#'
#' @param g a gct object
#' @param rid a vector of character ids or integer indices for ROWS
#' @param cid a vector of character ids or integer indices for COLUMNS
#' @examples
#' # first 10 rows and columns by index
#' (a <- subset.gct(ds, rid=1:10, cid=1:10))
#'
#' # first 10 rows and columns using character ids
#' (b <- subset.gct(ds, rid=ds@rid[1:10], cid=ds@cid[1:10]))
#'
#' identical(a, b) # TRUE
#'
#' @family GCT utilities
#' @export
setGeneric("subset.gct", function(g, rid=NULL, cid=NULL) {
    standardGeneric("subset.gct")
})
setMethod("subset.gct", signature("GCT"),
          function(g, rid, cid) {
              # ids can either be a vector of character strings corresponding
              # to row / column ids in the gct object, or integer vectors
              # corresponding to row / column indices
              if (is.null(rid)) rid <- g@rid
              if (is.null(cid)) cid <- g@cid
              # see whether we were given characters or integers
              # and handle accordingly
              process_ids <- function(ids, ref_ids, param) {
                  # simple helper function to handle id/idx conversion
                  # for character or integer ids
                  if (is.character(ids)) {
                      idx <- which(ref_ids %in% ids)
                  } else if (all(is.wholenumber(ids))) {
                      idx <- ids
                      ids <- ref_ids[idx]
                  } else {
                      stop(paste(param, "must be character or ingeter"))
                  }
                  return(list(ids=ids, idx=idx))
              }
              processed_rid <- process_ids(rid, g@rid, "rid")
              processed_cid <- process_ids(cid, g@cid, "cid")
              rid <- processed_rid$ids
              ridx <- processed_rid$idx
              cid <- processed_cid$ids
              cidx <- processed_cid$idx
              sdrow <- setdiff(rid, g@rid)
              sdcol <- setdiff(cid, g@cid)
              if (length(sdrow) > 0) warning("the following rids were not found:\n", paste(sdrow, collapse="\n"))
              if (length(sdcol) > 0) warning("the following cids were not found:\n", paste(sdcol, collapse="\n"))
              newg <- g
              # make sure ordering is right
              rid <- g@rid[ridx]
              cid <- g@cid[cidx]
              newg@mat <- matrix(g@mat[ridx, cidx], nrow=length(rid), ncol=length(cid))
              colnames(newg@mat) <- cid
              rownames(newg@mat) <- rid
              # cdesc <- data.frame(g@cdesc)
              # rdesc <- data.frame(g@rdesc)
              # make sure annotations row ordering matches
              # matrix, rid, and cid
              newg@cdesc <- subset_to_ids(g@cdesc, cid)
              newg@rdesc <- subset_to_ids(g@rdesc, rid)
              newg@rid <- rid
              newg@cid <- cid
              if (any(dim(newg@mat) == 0)) {
                  warning("one or more returned dimension is length 0
                          check that at least some of the provided rid and/or
                          cid values have matches in the GCT object supplied")
              }
              return(newg)
              })

#' Merge two GCT objects together
#'
#' @param g1 the first GCT object
#' @param g2 the second GCT object
#' @param dimension the dimension on which to merge (row or column)
#' @param matrix_only boolean idicating whether to keep only the
#'   data matrices from \code{g1} and \code{g2} and ignore their
#'   row and column meta data
#' @examples
#' # take the first 10 and last 10 rows of an object
#' # and merge them back together
#' (a <- subset.gct(ds, rid=1:10))
#' (b <- subset.gct(ds, rid=969:978))
#' (merged <- merge.gct(a, b, dimension="row"))
#'
#' @family GCT utilities
#' @export
setGeneric("merge.gct", function(g1, g2, dimension="row", matrix_only=F) {
    standardGeneric("merge.gct")
})
setMethod("merge.gct", signature("GCT", "GCT"),
          function(g1, g2, dimension, matrix_only) {
              # given two gcts objects g1 and g2, merge them
              # on the specified dimension
              if (dimension == "column") dimension <- "col"
              if (dimension == "row") {
                  message("appending rows...")
                  newg <- g1
                  # we're just appending rows so don't need to do anything
                  # special with the rid or rdesc. just cat them
                  newg@rid <- c(g1@rid, g2@rid)
                  newg@rdesc <- data.frame(rbind(data.table::data.table(g1@rdesc), data.table::data.table(g2@rdesc), fill=T))
                  # need figure out the index for how to sort the columns of
                  # g2@mat so that they are in sync with g1@mat
                  idx <- match(g1@cid, g2@cid)
                  newg@mat <- rbind(g1@mat, g2@mat[, idx])
                  if (!matrix_only) {
                      # apply the same sort order to the rows of g2@cdesc so that
                      # it's in sync with the final merged matrix
                      # figure out which fields are common and keep from the first gct
                      cmn_names <- intersect(names(g1@cdesc), names(g2@cdesc))
                      newg@cdesc <- cbind(g1@cdesc, g2@cdesc[idx, !(names(g2@cdesc) %in% cmn_names)])
                  } else {
                      newg@cdesc <- data.frame()
                  }
              }
              else if (dimension == "col") {
                  message("appending columns...")
                  newg <- g1
                  # we're just appending columns so don't need to do anything
                  # special with cid or cdesc. just cat them
                  newg@cid <- c(g1@cid, g2@cid)
                  newg@cdesc <- data.frame(rbind(data.table::data.table(g1@cdesc), data.table::data.table(g2@cdesc), fill=T))
                  # need figure out the index for how to sort the rows of
                  # g2@mat so that they are in sync with g1@mat
                  idx <- match(g1@rid, g2@rid)
                  newg@mat <- cbind(g1@mat, g2@mat[idx, ])
                  if (!matrix_only) {
                      # apply the same sort order to the rows of g2@rdesc so that
                      # it's in sync with the final merged matrix
                      # figure out which fields are common and keep from the first gct
                      cmn_names <- intersect(names(g1@rdesc), names(g2@rdesc))
                      newg@rdesc <- cbind(g1@rdesc, g2@rdesc[idx, !(names(g2@rdesc) %in% cmn_names)])
                  } else {
                      newg@rdesc <- data.frame()
                  }
              } else {
                  stop("dimension must be either row or col")
              }
              return(newg)
          })


#' Merge two \code{\link{data.frame}}s, but where there are common fields
#' those in \code{x} are retained and those in \code{y} are dropped.
#'
#' @param x the \code{\link{data.frame}} whose columns take precedence
#' @param y another \code{\link{data.frame}}
#' @param by a vector of column names to merge on
#' @param allow.cartesian boolean indicating whether it's ok
#'   for repeated values in either table to merge with each other
#'   over and over again.
#' @param as_data_frame boolean indicating whether to ensure
#'   the returned object is a \code{\link{data.frame}} instead of a \code{\link{data.table}}.
#'   This ensures compatibility with GCT object conventions,
#'   that is, the \code{\link{rdesc}} and \code{\link{cdesc}} slots must be strictly
#'   \code{\link{data.frame}} objects.
#'
#' @return a \code{\link{data.frame}} or \code{\link{data.table}} object
#'
#' @examples
#' (x <- data.table(foo=letters[1:10], bar=1:10))
#' (y <- data.table(foo=letters[1:10], bar=11:20, baz=LETTERS[1:10]))
#' # the 'bar' column from y will be dropped on merge
#' cmapR:::merge_with_precedence(x, y, by="foo")
#'
#' @keywords internal
#' @seealso data.table::merge
merge_with_precedence <- function(x, y, by, allow.cartesian=T,
                                  as_data_frame = T) {
    trash <- check_colnames(by, x)
    trash <- check_colnames(by, y)
    # cast as data.tables
    x <- data.table::data.table(x)
    y <- data.table::data.table(y)
    # get rid of row names
    data.table::setattr(x, "rownames", NULL)
    data.table::setattr(y, "rownames", NULL)
    common_cols <- intersect(names(x), names(y))
    y_keepcols <- unique(c(by, setdiff(names(y), common_cols)))
    y <- y[, y_keepcols, with=F]
    # if not all ids match, issue a warning
    if (!all(x[[by]] %in% y[[by]])) {
        warning("not all rows of x had a match in y. some columns may contain NA")
    }
    # merge keeping all the values in x, making sure that the
    # resulting data.table is sorted in the same order as the
    # original object x
    merged <- merge(x, y, by=by, allow.cartesian=allow.cartesian, all.x=T)
    if (as_data_frame) {
        # cast back to a data.frame if requested
        merged <- data.frame(merged)
    }
    return(merged)
}


#' Add annotations to a GCT object
#'
#' @description Given a GCT object and either a \code{\link{data.frame}} or
#' a path to an annotation table, apply the annotations to the
#' gct using the given \code{keyfield}.
#'
#' @param g a GCT object
#' @param annot a \code{\link{data.frame}} or path to text table of annotations
#' @param dimension either 'row' or 'column' indicating which dimension
#'   of \code{g} to annotate
#' @param keyfield the character name of the column in \code{annot} that
#'   matches the row or column identifiers in \code{g}
#'
#' @return a GCT object with annotations applied to the specified
#'   dimension
#'
#' @examples
#' \dontrun{
#'  g <- parse.gctx('/path/to/gct/file')
#'  g <- annotate.gct(g, '/path/to/annot')
#' }
#'
#' @family GCT utilities
#' @export
setGeneric("annotate.gct", function(g, annot, dimension="row", keyfield="id") {
    standardGeneric("annotate.gct")
})
setMethod("annotate.gct", signature("GCT"), function(g, annot, dimension,
                                                     keyfield) {
    if (!(any(class(annot) == "data.frame"))) {
        # given a file path, try to read it in
        annot <- fread(annot)
    } else {
        # convert to data.table
        annot <- data.table::data.table(annot)
    }
    # convert the keyfield column to id for merging
    # assumes the gct object has an id field in its existing annotations
    if (!(keyfield %in% names(annot))) {
        stop(paste("column", keyfield, "not found in annotations"))
    }
    # rename the column to id so we can do the merge
    annot$id <- annot[[keyfield]]
    if (dimension == "column") dimension <- "col"
    if (dimension == "row") {
        orig_id <- g@rdesc$id
        merged <- merge_with_precedence(g@rdesc, annot, by="id",
                                        allow.cartesian=T, as_data_frame=T)
        idx <- match(orig_id, merged$id)
        merged <- merged[idx, ]
        g@rdesc <- merged
    } else if (dimension == "col") {
        orig_id <- g@cdesc$id
        merged <- merge_with_precedence(g@cdesc, annot, by="id",
                                        allow.cartesian=T, as_data_frame=T)
        idx <- match(orig_id, merged$id)
        merged <- merged[idx, ]
        g@cdesc <- merged
    } else {
        stop("dimension must be either row or column")
    }
    return(g)
})

#' Transpose a GCT object
#'
#' @param g the \code{GCT} object
#'
#' @return a modified verion of the input \code{GCT} object
#'   where the matrix has been transposed and the row and column
#'   ids and annotations have been swapped.
#'
#' @examples
#' transpose.gct(ds)
#'
#' @family GCT utilties
#' @export
setGeneric("transpose.gct", function(g) standardGeneric("transpose.gct"))
setMethod("transpose.gct", signature("GCT"), function(g) {
    # transpose matrix
    g@mat <- t(g@mat)
    # create new data
    rid.new <- g@cid
    cid.new <- g@rid
    rdesc.new <- g@cdesc
    cdesc.new <- g@rdesc
    # overwrite g
    g@rid <- rid.new
    g@cid <- cid.new
    g@rdesc <- rdesc.new
    g@cdesc <- cdesc.new
    return(g)
})
