library(scp)
library(tidyverse)
library(UpSetR)
library(muscat)
library(msqrob2)

sceBalanced <- readRDS("data/simulatedData/simulatedData_sampling_balanced_leduc.rds")

sceRes <- scpModelWorkflow(sceBalanced, ~ Mock)


daRes <- scpDifferentialAnalysis(
    sceRes,
    contrasts = list(c("Mock", "mock2", "mock1"))
)[[1]] %>%
    as.data.frame() %>%
    filter(padj <= 0.05)

resNames <- daRes$feature

trueNames <- row.names(rowData(sceRes)[rowData(sceRes)$is_DA_mock1_vs_mock2, ])

listUpSet <- list("model" = resNames, "real" = trueNames)
upset(fromList(listUpSet), order.by = "freq")

aggregateColumns <- function(object, fcol, fun, ...) {
    .makePeptideProteinVector <- function(m, collapse = ";") {
        stopifnot(inherits(m, "Matrix"))
        vec <- rep(NA_character_, nrow(m))
        for (i in seq_len(nrow(m)))
            vec[i] <- paste(names(which(m[i, ] != 0)), collapse = collapse)
        names(vec) <- rownames(m)
        vec
    }
    if (missing(fcol))
        stop("'fcol' is required.")
    m <- assay(object, 1)
    cd <- colData(object)
    if (!fcol %in% names(cd))
        stop("'fcol' not found in the assay's colData.")
    groupBy <- cd[[fcol]]

    ## Store class of assay i in case it is not a SummarizedExperiment
    ## so that the aggregated assay can be reverted to that class
    .class <- class(object)

    ## Message about NA values is quant/row data
    has_na <- character()
    if (anyNA(m))
        has_na <- c(has_na, "quantitative")
    if (anyNA(cd, recursive = TRUE))
        has_na <- c(has_na, "col")
    if (length(has_na)) {
        msg <- paste(paste("Your", paste(has_na, collapse = " and "),
                           " data contain missing values."),
                     "Please read the relevant section(s) in the",
                     "aggregateSamples manual page regarding the",
                     "effects of missing values on data aggregation.")
        message(paste(strwrap(msg), collapse = "\n"))
    }
    if (is.vector(groupBy) && !is.list(groupBy)) { ## atomic vectors
        aggregated_assay <- aggregate_columns_by_vector(m, groupBy, fun, ...)
        aggcount_assay <- aggregate_columns_by_vector(m, groupBy, rowCounts)
        aggregated_coldata <- QFeatures::reduceDataFrame(cd, cd[[fcol]],
                                                         simplify = TRUE,
                                                         drop = TRUE,
                                                         count = TRUE)<
            assays <- SimpleList(assay = aggregated_assay, aggcounts = aggcount_assay)
        coldata <- aggregated_coldata[colnames(aggregated_assay), , drop = FALSE]
    } else if (is(groupBy, "Matrix")) {
        aggregated_assay <- aggregate_columns_by_matrix(m, groupBy, fun, ...)
        ## Remove the adjacency matrix that should be dropped anyway
        cd[[fcol]] <- NULL
        ## Temp variable for unfolding and reducing - removed later
        rd[["._vec_"]] <- .makePeptideProteinVector(groupBy)
        rd <- unfoldDataFrame(rd, "._vec_")
        aggregated_coldata <- reduceDataFrame(rd, rd[["._vec_"]], drop = TRUE)
        aggregated_coldata[["._vec_"]] <- NULL
        ## Count the number of samples per conditions
        .n <- apply(groupBy != 0, 2, sum)
        aggregated_coldata[[".n"]] <- .n[rownames(aggregated_coldata)]

        assays <- SimpleList(assay = as.matrix(aggregated_assay)) ## to discuss
        coldata <- aggregated_coldata[rownames(aggregated_assay), , drop = FALSE]
    } else stop("'fcol' must refer to an atomic vector or a sparse matrix.")
    print(dim(coldata))
    print(dim(assays))
    se <- SummarizedExperiment(assays = assays,
                               colData = coldata,
                               rowData = rowData(object))

    ## If the input objects weren't SummarizedExperiments, then try to
    ## convert the merged assay into that class. If the conversion
    ## fails, keep the SummarizedExperiment, otherwise use the
    ## converted object (see issue #78).
    if (.class != "SummarizedExperiment")
        se <- tryCatch(as(se, .class),
                       error = function(e) se)
    return(se)
}

aggregate_columns_by_vector <- function (x, INDEX, FUN, ...)
{
    if (!(is.matrix(x) | inherits(x, "HDF5Matrix")))
        stop("'x' must be a matrix or an object that inherits from ",
             "'HDF5Matrix'.")
    if (!identical(length(INDEX), ncol(x)))
        stop("The length of 'INDEX' has to be identical to 'ncol(x).")
    FUN <- match.fun(FUN)
    res <- lapply(split(seq_len(ncol(x)), INDEX), FUN = function(i) FUN(x[, i
                                                                          , drop = FALSE], ...))
    nms <- names(res)
    res <- do.call(cbind, res)
    colnames(res) <- nms
    rownames(res) <- rownames(x)
    if (inherits(x, "HDF5Matrix"))
        res <- HDF5Array::writeHDF5Array(res, filepath = HDF5Array::path(x),
                                         with.dimnames = TRUE)
    res
}

aggregate_columns_by_matrix <- function (x, MAT, FUN, ...)
{
    if (!(is.matrix(x) | inherits(x, "HDF5Matrix")))
        stop("'x' must be a matrix or an object that inherits from ",
             "'HDF5Matrix'.")
    if (!is(MAT, "Matrix") && !is(MAT, "matrix"))
        stop("'MAT' must be a matrix.")
    if (!identical(ncol(MAT), ncol(x)))
        stop("nrow(MAT) must be identical to 'nrow(x).")
    res <- do.call(FUN, list(x, MAT, ...))
    if (!is.null(colnames(MAT)) && !identical(colnames(MAT),
                                              colnames(res)))
        stop("The colum names of 'MAT' have to be identical to the colnames of 'res'!")
    if (!is.null(colnames(x)) && !identical(colnames(x), colnames(res)))
        stop("The column names of 'x' have to be identical to the column names of 'res'!")
    rownames(res) <- rownames(x)
    if (inherits(x, "HDF5Matrix"))
        res <- HDF5Array::writeHDF5Array(res, filepath = HDF5Array::path(x),
                                         with.dimnames = TRUE)
    res
}

colData(sceBalanced)[["numMock"]] <- as.numeric(colData(sceBalanced)[["Mock"]])
aggregateColumns(object = sceBalanced, fcol = "numMock", fun = MatrixGenerics::rowMeans)

assay(test)
