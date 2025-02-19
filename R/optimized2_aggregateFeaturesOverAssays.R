library(profvis)
library(QFeatures)
library(scp)
library(scpdata)
library(MsCoreUtils)

source("R/generate_data.R")

.aggregateQFeaturesCopy <- function(object, fcol, fun, rowDataCols = NULL, ...) {
    ## Copied from the PSMatch package, given that it is not available
    ## on Bioconductor yet.
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
    if(is.null(rowDataCols)) rd <- rowData(object) else rd <- rowData(object)[, rowDataCols]
    if (!fcol %in% names(rd))
        stop("'fcol' not found in the assay's rowData.")
    groupBy <- rd[[fcol]]

    ## Store class of assay i in case it is not a SummarizedExperiment
    ## so that the aggregated assay can be reverted to that class
    .class <- class(object)

    ## Message about NA values is quant/row data
    has_na <- character()
    if (anyNA(m))
        has_na <- c(has_na, "quantitative")
    if (anyNA(rd, recursive = TRUE))
        has_na <- c(has_na, "row")
    if (length(has_na)) {
        msg <- paste(paste("Your", paste(has_na, collapse = " and "),
                           " data contain missing values."),
                     "Please read the relevant section(s) in the",
                     "aggregateFeatures manual page regarding the",
                     "effects of missing values on data aggregation.")
        message(paste(strwrap(msg), collapse = "\n"))
    }

    if (is.vector(groupBy) & !is.list(groupBy)) { ## atomic vectors
        aggregated_assay <- aggregate_by_vector(m, groupBy, fun, ...)
        aggcount_assay <- aggregate_by_vector(m, groupBy, colCounts)
        aggregated_rowdata <- QFeatures::reduceDataFrame(rd, rd[[fcol]],
                                                         simplify = TRUE,
                                                         drop = TRUE,
                                                         count = TRUE)
        assays <- SimpleList(assay = aggregated_assay, aggcounts = aggcount_assay)
        rowdata <- aggregated_rowdata[rownames(aggregated_assay), , drop = FALSE]
    } else if (is(groupBy, "Matrix")) {
        aggregated_assay <- aggregate_by_matrix(m, groupBy, fun, ...)
        ## Remove the adjacency matrix that should be dropped anyway
        rd[[fcol]] <- NULL
        ## Temp variable for unfolding and reducing - removed later
        rd[["._vec_"]] <- .makePeptideProteinVector(groupBy)
        rd <- unfoldDataFrame(rd, "._vec_")
        aggregated_rowdata <- reduceDataFrame(rd, rd[["._vec_"]], drop = TRUE)
        aggregated_rowdata[["._vec_"]] <- NULL
        ## Count the number of peptides per protein
        .n <- apply(groupBy != 0, 2, sum)
        aggregated_rowdata[[".n"]] <- .n[rownames(aggregated_rowdata)]

        assays <- SimpleList(assay = as.matrix(aggregated_assay)) ## to discuss
        rowdata <- aggregated_rowdata[rownames(aggregated_assay), , drop = FALSE]
    } else stop("'fcol' must refer to an atomic vector or a sparse matrix.")
    se <- SummarizedExperiment(assays = assays,
                               colData = colData(object),
                               rowData = rowdata)

    ## If the input objects weren't SummarizedExperiments, then try to
    ## convert the merged assay into that class. If the conversion
    ## fails, keep the SummarizedExperiment, otherwise use the
    ## converted object (see issue #78).
    if (.class != "SummarizedExperiment")
        se <- tryCatch(as(se, .class),
                       error = function(e) se)
    return(se)
}

invariant_col2 <- function(x) all(lengths(lapply(x, unique)) == 1)

invariant_cols2 <- function(x) {
    res <- rep(NA, length(x[[1]]))
    for (i in seq_along(res)) {
        res[i] <- invariant_col2(x[, i])
    }
    which(res)
}

aggregateFeaturesOverAssaysCopy <- function(object, i, fcol, name, fun, ...) {
    if (length(i) != length(name)) stop("'i' and 'name' must have same length")
    if (length(fcol) == 1) fcol <- rep(fcol, length(i))
    if (length(i) != length(fcol)) stop("'i' and 'fcol' must have same length")
    if (is.numeric(i)) i <- names(object)[i]

    ## Compute the aggregated assays
    el <- experiments(object)[i]
    rowDataColsKept <- colnames(rowData(el[[1]]))
    colKept <- lapply(X = seq_along(el), FUN = function(j) {
        rd <- rowData(el[[j]])
        rd <- split(rd, rd[[fcol[j]]])
        invariant_cols2(rd)
    })
    colKept <- unique(unlist(colKept))
    for (j in seq_along(el)) {
        suppressMessages(
            el[[j]] <- .aggregateQFeaturesCopy(el[[j]], fcol = fcol[j],
                                         fun = fun, rowDataCols = colKept, ...)
        )
        ## Print progress
        message("\rAggregated: ", j, "/", length(el), "\n")
    }
    names(el) <- name
    ## Get the AssayLinks for the aggregated assays
    alnks <- lapply(seq_along(i), function(j) {
        hits <- QFeatures:::.get_Hits(rdFrom = rowData(object[[i[j]]]),
                                      rdTo = rowData(el[[j]]),
                                      varFrom = fcol[[j]],
                                      varTo = fcol[[j]])
        AssayLink(name = name[j], from = i[j], fcol = fcol[j], hits = hits)
    })
    ## Append the aggregated assays and AssayLinks to the previous assays
    el <- c(object@ExperimentList, el)
    alnks <- append(object@assayLinks, AssayLinks(alnks))
    ## Update the sampleMapfrom the data
    smap <- MultiAssayExperiment:::.sampleMapFromData(colData(object), el)
    ## Create the new QFeatures object
    new("QFeatures",
        ExperimentList = el,
        colData = colData(object),
        sampleMap = smap,
        metadata = metadata(object),
        assayLinks = alnks)
}


leduc <- scpdata::leduc2022_pSCoPE()

PSM1000 <- generateTMTPSM(leduc, 1000)

opti <- function() {
    aggregateFeaturesOverAssaysCopy(PSM1000,
                                    i = seq_along(PSM1000),
                                    "modseq",
                                    name = paste0("peptide_", seq_along(PSM1000)),
                                    fun = MsCoreUtils::robustSummary)}


vanilla <- function() {
    aggregateFeaturesOverAssays(PSM1000,
                                 i = seq_along(PSM1000),
                                 "modseq",
                                 name = paste0("peptide_", seq_along(PSM1000)),
                                 fun = MsCoreUtils::robustSummary)
}
bench::mark({aggregateFeaturesOverAssaysCopy(PSM1000,
                                i = seq_along(PSM1000),
                                "modseq",
                                name = paste0("peptide_", seq_along(PSM1000)),
                                fun = MsCoreUtils::robustSummary)},
            {aggregateFeaturesOverAssays(PSM1000,
                                             i = seq_along(PSM1000),
                                             "modseq",
                                             name = paste0("peptide_", seq_along(PSM1000)),
                                             fun = MsCoreUtils::robustSummary)}, memory = FALSE, check = FALSE)


# Check of each col before the aggregation => then subset (unique(all the cols))


PSM1000Joined <- joinAssays(PSM1000, i = seq_along(PSM1000))
aggregateFeatures(PSM1000Joined, i = "joinedAssay", fcol = "modseq")
bench::mark(aggregateFeatures(PSM1000Joined, i = "joinedAssay", fcol = "modseq"),
            {aggregateFeaturesOverAssays(PSM1000,
                                         i = seq_along(PSM1000),
                                         "modseq",
                                         name = paste0("peptide_", seq_along(PSM1000)),
                                         fun = MsCoreUtils::robustSummary)}, memory = FALSE, check = FALSE)
