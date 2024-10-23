library(peakRAM)
library(scpdata)
library(scp)
library(SingleCellExperiment)
library(bench)

source("R/vignette_leduc2022_script.R")

leduc2022Benchmark <- function(nCellRange, rep) {
    base <- scpdata::leduc2022()
    press(
        nCell = nCellRange,
        rep = rep,
        {
            set.seed(123)
            print0("Starting generation of data for ",
                   nCell,
                   " Cellules...")
            qfeatures <- leduc2022Generate(base, nCell)
            bench::mark(leduc2022script(qfeatures), memory = TRUE)
        }
    )
}

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

        noise <- pmin(matrix(rnorm(n = length(newAssay), mean = 0, sd = 5),
                        nrow = nrow(newAssay),
                        ncol = ncol(newAssay)), 0)
        noisyNewAssay <- newAssay + noise

        newColData <- colData(original_se)
        newColData$Set <- paste0("run_", i)
        newColData$Channel <- as.character(newColData$Channel)
        newColData$IsolationTimeStamp <- as.character(newColData$IsolationTimeStamp)
        rownames(newColData) <- newSampleNames

        newRowData <- rowData(original_se)
        rownames(newRowData) <- newFeaturesNames

        psmCounter <- psmCounter + nrow(newRowData)
        noisy_se <- SingleCellExperiment(assays = list(noisyNewAssay),
                                         rowData = newRowData,
                                         colData = newColData)
        new_qfeatures <- addAssay(new_qfeatures, noisy_se, name = paste0("run_", i))
    }
    return(new_qfeatures)
}

sd <- benchmarkme::get_sys_details()
cat("Machine: ", sd$sys_info$sysname, " (", sd$sys_info$release, ")\n",
    "R version: R.", sd$r_version$major, ".", sd$r_version$minor,
    " (svn: ", sd$r_version$`svn rev`, ")\n",
    "RAM: ", round(sd$ram / 1E9, 1), " GB\n",
    "CPU: ", sd$cpu$no_of_cores, " core(s) - ", sd$cpu$model_name, "\n",
    sep = "")
