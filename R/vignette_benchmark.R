library(peakRAM)
library(scpdata)
library(scp)
library(SingleCellExperiment)

source("R/vignette_leduc2022_script.R")

leduc2022Benchmark <- function(nCellRange,
    nreplicates,
    outputDir = "dataOutput/vignetteBenchmark") {
    if (file.exists(file.path(outputDir, "size_report.tsv"))) {
        unlink(file.path(outputDir, "size_report.tsv"))
        cat(file.path(outputDir, "size_report.tsv"), " has been deleted..\n")
    }
    if (file.exists(file.path(outputDir, "memoryOutput"))) {
        unlink(file.path(outputDir, "memoryOutput"), recursive = TRUE)
        cat(file.path(outputDir, "memoryOutput"), " has been deleted\n")
    }
    dir.create(file.path(outputDir, "memoryOutput"))
    write.table(
        data.frame(
            nCell = integer(),
            rep = integer(),
            size = integer()
        ),
        file = file.path(outputDir, "size_report.tsv"),
        append = FALSE,
        col.names = TRUE,
        row.names = FALSE
    )

    base <- scpdata::leduc2022()
    for (rep in seq(nreplicates)) {
        for (nCell in nCellRange) {
            set.seed(123)
            print(paste0(
                "Starting generation of data for ",
                nCell,
                " Cellules, replicate #",
                rep
            ))
            leduc <- leduc2022Generate(base, nCell)
            print(paste0(
                "Starting benchmarking for ",
                nCell,
                " Cellules, replicate #",
                rep
            ))
            res <- peakRAM(
                assaysNames <- names(leduc),
                # Remove contaminant, noisy and low-confidence spectra
                leduc <- leducFilterFeatures(leduc),
                # Sample to carrier filter
                leduc <- leducSCR(leduc),
                # Filter on summed single-cell signal
                leduc <- leducScSums(leduc),
                # Normalize to reference
                leduc <- leducNormToRef(leduc),
                # Aggregate PSM data to peptide data
                leduc <- leducAggPSM(leduc),
                # Consensus mapping of peptides to proteins
                leduc <- leducConsensus(leduc, assaysNames),
                # Cleaning missing data
                leduc <- leducMissingData(leduc, assaysNames),
                # Join assays
                leduc <- leducJoin(leduc, assaysNames),
                # Filter single-cells based on median CV
                leduc <- leducFilterCV(leduc, assaysNames),
                # Normalization peptides
                leduc <- leducNormPep(leduc),
                # Missing data filtering
                leduc <- leducFilteringNA(leduc),
                # Log-transformation
                leduc <- leducLogTransfo(leduc),
                # Aggregate peptide data to protein data
                leduc <- leducAggPep(leduc),
                # Normalization proteins
                leduc <- leducNormPro(leduc),
                # Imputation
                leduc <- leducImpute(leduc),
                # Batch correction
                leduc <- leducBatch(leduc),
                # Normalization batch corrected proteins
                leduc <- leducNormBatch(leduc),
                # PCA
                # leduc <- leducPCA(leduc)
                write.table(
                    data.frame(
                        nCell = as.integer(nCell),
                        rep = as.integer(rep),
                        size = object.size(leduc)
                    ),
                    file = file.path(outputDir, "size_report.tsv"),
                    append = TRUE,
                    col.names = FALSE,
                    row.names = FALSE
                )
            )
            write.csv(res, file = file.path(
                outputDir, "memoryOutput",
                paste0(paste(
                    nCell, rep,
                    sep = "_"
                ), ".csv")
            ))
            gc()
        }
    }
    sd <- benchmarkme::get_sys_details()
    if (file.exists(file.path(outputDir, "hardware_software.txt"))) {
        unlink(file.path(outputDir, "hardware_software.txt"))
    }
    sink(file.path(outputDir, "hardware_software.txt"))
    cat("Machine: ", sd$sys_info$sysname, " (", sd$sys_info$release, ")\n",
        "R version: R.", sd$r_version$major, ".", sd$r_version$minor,
        " (svn: ", sd$r_version$`svn rev`, ")\n",
        "RAM: ", round(sd$ram / 1E9, 1), " GB\n",
        "CPU: ", sd$cpu$no_of_cores, " core(s) - ", sd$cpu$model_name, "\n",
        sep = ""
    )
    cat(" ----------------------------------------------------------------- \n")
    print(sessionInfo())
    sink()
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
        noisy_se <- SingleCellExperiment(
            assays = list(noisyNewAssay),
            rowData = newRowData,
            colData = newColData
        )
        new_qfeatures <- addAssay(new_qfeatures, noisy_se, name = paste0("run_", i))
    }
    return(new_qfeatures)
}
