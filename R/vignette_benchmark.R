#' leduc2022 vignette benchmarking
#'
#' @param nCellRange range of number of cells to benchmark
#'
#' @return results
#'
#' @importFrom scpdata leduc2022
#' @export
#'
leduc2022Benchmark <- function(nCellRange) {
    results <- list()
    baseLeduc <- leduc2022()
    for (nCell in nCellRange) {
        qfeatures <- leduc2022Generate(baseLeduc, nCell)
        result <- renderBenchmarking(qfeatures, nCell)
        results[[as.character(nCell)]] <- result
    }
    return(results)
}

#' Title
#'
#' @param qfeaturesObject qfeatures object
#' @param nCell number of cells
#'
#' @return time
#' @importFrom rmarkdown render
renderBenchmarking <- function(qfeaturesObject, nCell) {
    rmdPath <- system.file("rmd/leduc2022_benchmark.Rmd",
                       package = "benchmarkQFeatures")
    start <- Sys.time()
    rmarkdown::render(rmdPath,
                      output_file = paste0(nCell,"_leduc2022_benchmark.html"),
                      params = list(qfeatures = qfeaturesObject))
    end <- Sys.time()
    return(end - start)
}

#' Generate leduc2022 data
#'
#' @param nCell number of cells to generate
#'
#' @return qfeatures
#' @import scp
#'
leduc2022Generate <- function(base, nCell) {
    base <- base[, , -(135:138)]
    nRun <- nCell %/% 18
    nAssays <- length(base)
    if (nRun <= nAssays) {
        sampledAssays <- sample(seq_len(nAssays), nRun, replace = FALSE)
    } else {
        sampledAssays <- sample(seq_len(nAssays), nRun, replace = TRUE)
    }
    new_qfeatures <- base[, , 0]
    psmCounter <- 0
    for (i in seq_len(nRun)) {
        assay_idx <- sampledAssays[i]
        original_se <- getWithColData(base, assay_idx)
        newSampleNames <- paste0("run_", i, "_RI", seq_len(18))
        newAssay <- assay(original_se)
        newFeaturesNames <- paste0("PSM_",
                                   seq(from = psmCounter + 1,
                                       length.out = nrow(newAssay)))
        colnames(newAssay) <- newSampleNames
        rownames(newAssay) <- newFeaturesNames

        noise <- matrix(rnorm(n = length(newAssay), mean = 0, sd = 5),
                        nrow = nrow(newAssay),
                        ncol = ncol(newAssay))
        noisyNewAssay <- newAssay + noise

        newColData <- colData(original_se)
        newColData$Set <- paste0("run_", i)
        rownames(newColData) <- newSampleNames

        newRowData <- rowData(original_se)
        rownames(newRowData) <- newFeaturesNames

        psmCounter <- psmCounter + nrow(newRowData)
        noisy_se <- SummarizedExperiment(assays = list(noisyNewAssay),
                                         rowData = newRowData,
                                         colData = newColData)
        new_qfeatures <- addAssay(new_qfeatures, noisy_se, name = paste0("run_", i))
    }
    return(new_qfeatures)
}
