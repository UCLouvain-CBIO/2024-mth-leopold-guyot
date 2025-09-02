source(file.path("R", "minimumWorkflow.R"))
library(MultiAssayExperiment)
library(QFeatures)
library(peakRAM)
library(scp)
library(scpdata)
library(SingleCellExperiment)


leduc2022Generate <- function(base, nCell, SE) {
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
    newFeaturesNames <- paste0(
      "PSM_",
      seq(
        from = psmCounter + 1,
        length.out = nrow(newAssay)
      )
    )
    colnames(newAssay) <- newSampleNames
    rownames(newAssay) <- newFeaturesNames

    noise <- pmin(matrix(rnorm(n = length(newAssay), mean = 0, sd = 5),
                         nrow = nrow(newAssay),
                         ncol = ncol(newAssay)
    ), 0)
    noisyNewAssay <- newAssay + noise

    newColData <- colData(original_se)
    newColData$Set <- paste0("run_", i)
    newColData$Channel <- as.character(newColData$Channel)
    newColData$IsolationTimeStamp <- as.character(newColData$IsolationTimeStamp)
    rownames(newColData) <- newSampleNames

    newRowData <- rowData(original_se)
    rownames(newRowData) <- newFeaturesNames

    psmCounter <- psmCounter + nrow(newRowData)
    if (SE) {
      noisy_se <- SummarizedExperiment(
        assays = list(noisyNewAssay),
        rowData = newRowData,
        colData = newColData
      )
    } else {
      noisy_se <- SingleCellExperiment(
        assays = list(noisyNewAssay),
        rowData = newRowData,
        colData = newColData
      )
    }

    new_qfeatures <- addAssay(new_qfeatures, noisy_se, name = paste0("run_", i))
  }
  return(new_qfeatures)
}

benchWrapper <- function(replicate, nCell, subsetRowData = FALSE, SE) {
    results_list <- list()

    for (size in nCell) {
        leduc <- scpdata::leduc2022_pSCoPE()
        leduc <- leduc2022Generate(leduc, nCell = size, SE = SE)

        if (subsetRowData) {
            for (assay in seq_along(leduc)) {
                rowData(leduc[[assay]]) <- rowData(leduc[[assay]])[, c("dart_PEP",
                                                                       "PIF",
                                                                       "Proteins",
                                                                       "Leading.razor.protein",
                                                                       "dart_qval",
                                                                       "Potential.contaminant",
                                                                       "Reverse",
                                                                       "modseq")]
            }
        }

        results <- do.call(rbind, lapply(1:replicate, function(i) {
            res <- peakRAM(
                leduc <- filterRow(leduc),
                leduc <- filterCol(leduc),
                aggregatePSM(leduc)
            )
            res$replicate <- i
            res$subset_size <- size
            return(res)
        }))

        results_list[[as.character(size)]] <- results
    }

    return(do.call(rbind, results_list))
}

