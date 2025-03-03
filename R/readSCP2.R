library("SingleCellExperiment")

.readSingleCellExperiment2 <- function(assayData,
                                       quantCols = NULL,
                                       fnames = NULL,
                                       ecol = NULL, ...) {
    quantCols <- .checkWarnEcolCopy(quantCols, ecol)
    if (!is.vector(quantCols) || is.list(quantCols))
        stop("'quantCols' must be an atomics vector.")
    if (is.data.frame(assayData)) xx <- assayData
    else {
        args <- list(...)
        args$file <- assayData
        if ("rownames" %in% names(args)) {
            if (is.null(fnames)) fnames <- args$rownames
            args$rownames <- NULL
        }
        args$stringsAsFactors <- FALSE
        xx <- do.call(read.csv, args)
    }
    if (is.character(quantCols) || is.factor(quantCols)) {
        mis <- !quantCols %in% colnames(xx)
        if (any(mis))
            stop("Column identifiers ",
                 paste(quantCols[mis], collapse = ", "),
                 " not recognised among\n",
                 paste(colnames(xx), paste = ", "))
        quantCols <- which(colnames(xx) %in% quantCols)
    } else if (is.logical(quantCols)) {
        if (length(quantCols) != length(xx))
            stop("Length of 'quantCols' and 'assayData' do not match.")
        quantCols <- which(quantCols)
    }
    assay <- as.matrix(xx[, quantCols, drop = FALSE])
    fdata <- DataFrame(xx[, -quantCols, drop = FALSE])

    if (!missing(fnames)) {
        fnames <- fnames[1]
        if (is.numeric(fnames))
            fnames <- colnames(xx)[fnames]
        if (is.na(match(fnames, colnames(xx))))
            stop(fnames, " not found among\n",
                 paste(colnames(xx), paste = ", "))
        rownames(fdata) <- rownames(assay) <- fdata[, fnames]
    } else {
        rownames(fdata) <- rownames(assay) <- seq_len(nrow(assay))
    }
    SingleCellExperiment(assay, rowData = fdata)
}



##' @export
##'
##' @rdname readQFeatures
readSCP2 <- function(assayData,
                     colData = NULL,
                     quantCols = NULL,
                     runCol = NULL,
                     name = "quants",
                     removeEmptyCols = FALSE,
                     verbose = TRUE,
                     ecol = NULL,
                     ...) {
    if (verbose) message("Checking arguments.")
    assayData <- as.data.frame(assayData)
    if (!is.null(colData))
        colData <- data.frame(colData)
    quantCols <- .checkWarnEcolCopy(quantCols, ecol)
    quantCols <- .checkQuantColsCopy(assayData, colData, quantCols)
    runs <- .checkRunColCopy(assayData, colData, runCol)
    if (verbose) message("Loading data as a 'SummarizedExperiment' object.")
    se <- .readSingleCellExperiment2(assayData, quantCols, ...)
    rownames(se) <- make.unique(rownames(se))
    if (length(runs)) {
        if (verbose) message("Splitting data in runs.")
        el <- .splitSECopy(se, runs)
        el <- .createUniqueColnamesCopy(el, quantCols)
    } else {
        el <- structure(list(se), .Names = name[1])
    }
    if (removeEmptyCols) el <- .removeEmptyColumnsCopy(el)
    if (verbose) message("Formatting sample annotations (colData).")
    colData <- .formatColDataCopy(el, colData, runs, quantCols)
    if (verbose) message("Formatting data as a 'QFeatures' object.")
    QFeatures(experiments = el, colData = colData)
}


## ecol will be deprecated next release. This function warns if ecol
## is used (i.e. is not NULL), then sets quantCols with the value of
## ecol if quantCols wasn't used.
.checkWarnEcolCopy <- function(quantCols, ecol) {
    if (!is.null(ecol)) {
        if (!is.null(quantCols))
            stop("'quantCols' and 'ecols' can't be defined together. ",
                 "Use 'quantCols' only.")
        warning("'ecol' is deprecated, use 'quantCols' instead.")
        if (is.null(quantCols))
            quantCols <- ecol
    }
    quantCols
}

## This function will check the quantitation variable inputs. At the
## end, it will return a valid character quantCols, either as provided
## directly by the user, or generated from colData.
.checkQuantColsCopy <- function(assayData, colData, quantCols) {
    ## Fail early if both a missing
    if (is.null(colData) & is.null(quantCols))
        stop("Provide one of 'colData' or 'quantCols', both mustn't be NULL.")
    ## If we have a colData data.frame, it must contain a
    ## quantCols column and no quantCols should be provided.
    if (!is.null(colData)) {
        if (!"quantCols" %in% colnames(colData) &&
            length(quantCols) > 1)
            stop("'colData' must contain a column called 'quantCols'")
    }
    if (is.null(quantCols)) {
        ## if (is.null(colData))
        ##     stop("'quantCols' and 'colData' cannot both be NULL.")
        if (!"quantCols" %in% colnames(colData))
            stop("When 'quantCols' is NULL, 'colData' must ",
                 "contain a column called 'quantCols'.")
        quantCols <- unique(colData$quantCols)
    }
    if (is.numeric(quantCols) || is.logical(quantCols))
        quantCols <- colnames(assayData)[quantCols]
    mis <- quantCols[!quantCols %in% colnames(assayData)]
    if (length(mis))
        stop("Some column names in 'quantCols' are not found ",
             "in 'assayData': ", paste0(mis, collapse = ", "), ".")
    quantCols
}


## This function will check the batch/run variable inputs. At the end,
## it will return a vector of runs/batches or NULL (single-set
## case). Possible inputs combinations are:
##
## - `runCol` is NULL: single-set case
## - `runCol` only (i.e. `colAnnotion` is NULL), of length 1, refering
##   to a variable in `assayData`.
## - `runCol` and `colData`: in this case,
##   `colData$runCol` must exist.
.checkRunColCopy <- function(assayData, colData, runCol) {
    ## No runCol provided: single-set case
    if (is.null(runCol)) return(NULL)
    ## We have a runCol argument: multi-set case
    if (length(runCol) > 1)
        stop("'runCol' must contain the name of a single column ",
             "in 'assayData'.")
    if (!runCol %in% colnames(assayData))
        stop("'", runCol, "' (provided as 'runCol') not found ",
             "in 'assayData'.")
    runs <- assayData[[runCol]]
    if (!is.null(colData)) {
        ## We have a colData argument
        if (!"runCol" %in% colnames(colData))
            stop("When 'runCol' is not NULL, 'colData' must ",
                 "contain a column called 'runCol'.")
        mis <- !runs %in% colData$runCol
        if (any(mis)) {
            warning("Some runs are missing in 'colData': ",
                    paste0(unique(runs[mis]), collapse = ", "))
        }
    }
    assayData[[runCol]]
}

##' Split SummarizedExperiment into an ExperimentList
##'
##' The fonction creates an [ExperimentList] containing
##' [SummarizedExperiment] objects from a [SummarizedExperiment]
##' object (also works with [SingleCellExperiment] objects). `f` is
##' used to split `x`` along the rows (`f`` was a feature variable
##' name) or samples/columns (f was a phenotypic variable name). If f
##' is passed as a factor, its length will be matched to nrow(x) or
##' ncol(x) (in that order) to determine if x will be split along the
##' features (rows) or sample (columns). Hence, the length of f must
##' match exactly to either dimension.
##'
##' This function is not exported and was initially available as
##' scp::.splitSCE().
##'
##' @param x a single [SummarizedExperiment] object
##'
##' @param f a factor or a character of length 1. In the latter case,
##'     `f` will be matched to the row and column data variable names
##'     (in that order). If a match is found, the respective variable
##'     is extracted, converted to a factor if needed.
##' @noRd
.splitSECopy <- function(x, f) {
    ## Check that f is a factor
    if (length(f) == 1) {
        if (f %in% colnames(rowData(x))) {
            f <- rowData(x)[, f]
        }
        else if (f %in% colnames(colData(x))) {
            f <- colData(x)[, f]
        }
        else {
            stop("'", f, "' not found in rowData or colData")
        }
    }
    ## Check that the factor matches one of the dimensions
    if (!length(f) %in% dim(x))
        stop("length(f) not compatible with dim(x).")
    if (length(f) == nrow(x)) { ## Split along rows
        xl <- lapply(split(rownames(x), f = f), function(i) x[i, ])
    } else { ## Split along columns
        xl <- lapply(split(colnames(x), f = f), function(i) x[, i])
    }
    ## Convert list to an ExperimentList
    do.call(ExperimentList, xl)
}

.createUniqueColnamesCopy <- function(el, quantCols) {
    if (length(quantCols) == 1)  suffix <- ""
    else suffix <- paste0("_", quantCols)
    for (i in seq_along(el)) {
        colnames(el[[i]]) <- paste0(names(el)[[i]], suffix)
    }
    el
}

.removeEmptyColumnsCopy <- function(el) {
    for (i in seq_along(el)) {
        sel <- colSums(is.na(assay(el[[i]]))) != nrow(el[[i]])
        el[[i]] <- el[[i]][, sel]
    }
    el
}

## This function will create a colData from the different (possibly
## missing, i.e. NULL) arguments
.formatColDataCopy <- function(el, colData, runs, quantCols) {
    sampleNames <- unlist(lapply(el, colnames), use.names = FALSE)
    if (is.null(colData))
        return(DataFrame(row.names = sampleNames))
    if (!length(runs)) {
        rownames(colData) <- sampleNames
    } else {
        if (length(quantCols) == 1) {
            rownames(colData) <- colData$runCol
        } else {
            rownames(colData) <- paste0(colData$runCol, "_", colData$quantCols)
        }
    }
    colData <- colData[sampleNames, , drop = FALSE]
    rownames(colData) <- sampleNames ## clean NA in rownames
    colData
}


QFeatures <- function(..., assayLinks = NULL) {
    ans <- MultiAssayExperiment(...)
    if (isEmpty(ans)) assayLinks <- AssayLinks()
    else {
        if (is.null(assayLinks))
            assayLinks <- AssayLinks(names = names(ans))
    }
    new("QFeatures",
        ExperimentList = ans@ExperimentList,
        colData = ans@colData,
        sampleMap = ans@sampleMap,
        metadata = ans@metadata,
        assayLinks = assayLinks)
}
