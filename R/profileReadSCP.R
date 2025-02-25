library(scp)
data("mqScpData")

quantCols <- grep("Reporter.intensity.\\d", colnames(mqScpData),
                   value = TRUE)

data("sampleAnnotation")

scp <- readSCP(assayData = mqScpData,
               colData = sampleAnnotation,
               runCol = "Raw.file",
               removeEmptyCols = TRUE)

qfeat <- readQFeatures(assayData = mqScpData,
                       colData = sampleAnnotation,
                       runCol = "Raw.file",
                       removeEmptyCols = TRUE)

results <- microbenchmark::microbenchmark(readSCP(assayData = mqScpData,
                                       colData = sampleAnnotation,
                                       runCol = "Raw.file",
                                       removeEmptyCols = TRUE),
                               readQFeatures(assayData = mqScpData,
                                       colData = sampleAnnotation,
                                       runCol = "Raw.file",
                                       removeEmptyCols = TRUE))
profvis::profvis(readSCP(assayData = mqScpData,
                         colData = sampleAnnotation,
                         runCol = "Raw.file",
                         removeEmptyCols = TRUE))

profvis::profvis({
    readSCP <- function(...) {
        ans <- readQFeatures(...)
        el <- ExperimentList(lapply(experiments(ans),
                                    as, "SingleCellExperiment"))
        experiments(ans) <- el
        ans
    }
    readSCP(assayData = mqScpData,
            colData = sampleAnnotation,
            runCol = "Raw.file",
            removeEmptyCols = TRUE)
})
profvis::profvis({
readSingleCellExperimentOpt <- function(assayData,
                                        quantCols = NULL,
                                        fnames = NULL,
                                        ecol = NULL, ...) {
    quantCols <- QFeatures:::.checkWarnEcol(quantCols, ecol)
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


readSCPopt <- function(assayData,
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
    quantCols <- QFeatures:::.checkWarnEcol(quantCols, ecol)
    quantCols <- QFeatures:::.checkQuantCols(assayData, colData, quantCols)
    runs <- QFeatures:::.checkRunCol(assayData, colData, runCol)
    if (verbose) message("Loading data as a 'SummarizedExperiment' object.")
    se <- readSingleCellExperimentOpt(assayData, quantCols, ...)
    rownames(se) <- make.unique(rownames(se))
    if (length(runs)) {
        if (verbose) message("Splitting data in runs.")
        el <- QFeatures:::.splitSE(se, runs)
        el <- QFeatures:::.createUniqueColnames(el, quantCols)
    } else {
        el <- structure(list(se), .Names = name[1])
    }
    if (removeEmptyCols) el <- QFeatures:::.removeEmptyColumns(el)
    if (verbose) message("Formatting sample annotations (colData).")
    colData <- QFeatures:::.formatColData(el, colData, runs, quantCols)
    if (verbose) message("Formatting data as a 'QFeatures' object.")
    QFeatures(experiments = el, colData = colData)
}


scp <- readSCPopt(assayData = mqScpData,
               colData = sampleAnnotation,
               runCol = "Raw.file",
               removeEmptyCols = TRUE)
})
results <- microbenchmark::microbenchmark(readSCP(assayData = mqScpData,
                                                  colData = sampleAnnotation,
                                                  runCol = "Raw.file",
                                                  removeEmptyCols = TRUE),
                                          readSCPopt(assayData = mqScpData,
                                                  colData = sampleAnnotation,
                                                  runCol = "Raw.file",
                                                  removeEmptyCols = TRUE),
                                          readQFeatures(assayData = mqScpData,
                                                        colData = sampleAnnotation,
                                                        runCol = "Raw.file",
                                                        removeEmptyCols = TRUE))


source("R/readSCP2.R")

profvis::profvis(readSCP2(assayData = mqScpData,
                          colData = sampleAnnotation,
                          runCol = "Raw.file",
                          removeEmptyCols = TRUE))

source("R/readQFeatures2.R")

profvis::profvis(readQFeatures2(assayData = mqScpData,
                          colData = sampleAnnotation,
                          runCol = "Raw.file",
                          removeEmptyCols = TRUE))

results <- microbenchmark::microbenchmark(readSCP(assayData = mqScpData,
                                                  colData = sampleAnnotation,
                                                  runCol = "Raw.file",
                                                  removeEmptyCols = TRUE),
                                          readSCP2(assayData = mqScpData,
                                                     colData = sampleAnnotation,
                                                     runCol = "Raw.file",
                                                     removeEmptyCols = TRUE),
                                          readQFeatures(assayData = mqScpData,
                                                        colData = sampleAnnotation,
                                                        runCol = "Raw.file",
                                                        removeEmptyCols = TRUE),
                                          readQFeatures2(assayData = mqScpData,
                                                         colData = sampleAnnotation,
                                                         runCol = "Raw.file",
                                                         removeEmptyCols = TRUE))
