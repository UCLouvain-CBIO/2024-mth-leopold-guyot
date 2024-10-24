library(peakRAM)
library(scpdata)
library(scp)
library(SingleCellExperiment)
library(rmarkdown)

source("R/vignette_leduc2022_script.R")

leduc2022ResultGeneration <- function(sizeInputPath = "size_report.tsv",
                                      memoryUsageInputDir = "dataOutput/memoryOutput",
                                      outputFile = "leduc2022BenchmarkResults.html",
                                      outputDir = "reports"){
    rmarkdown::render("Rmd/leduc2022Results.rmd",
                      params = list(sizeInputPath = sizeInputPath,
                                    memoryUsageInputDir = memoryUsageInputDir),
                      output_file = outputFile,
                      output_dir = outputDir,
                      knit_root_dir = getwd())
}

leduc2022Benchmark <- function(nCellRange, rep) {
    base <- scpdata::leduc2022()
    results <- list()
    for (nCell in nCellRange){
        set.seed(123)
        print(paste0("Starting generation of data for ",
                     nCell,
                     " Cellules..."))
        qfeatures <- leduc2022Generate(base, nCell)
        print(paste0("Starting benchmarking for ",
                     nCell,
                     " Cellules..."))
        res <- peakRAM(leduc2022script(qfeatures))
        results[[as.character(nCell)]] <- res
        gc()
    }
    results
}

leduc2022BenchmarkDetails <- function(nCellRange,
                                      nreplicates,
                                      sizeOutputPath = "size_report.tsv",
                                      memoryUsageOutputDir = "dataOutput/memoryOutput") {
    if (file.exists(sizeOutputPath)) {
        unlink(sizeOutputPath)
        cat(sizeOutputPath, " has been deleted..\n")
    }
    if (file.exists(memoryUsageOutputDir)) {
        unlink(memoryUsageOutputDir,recursive = TRUE)
        cat(memoryUsageOutputDir, " has been deleted\n")
    }
    dir.create(memoryUsageOutputDir)
    write.table(data.frame(nCell = integer(),
                           rep = integer(),
                           size = integer()),
                file = sizeOutputPath,
                append = FALSE,
                col.names = TRUE,
                row.names = FALSE)

    base <- scpdata::leduc2022()
    for (rep in seq(nreplicates)){
        for (nCell in nCellRange){
            set.seed(123)
            print(paste0("Starting generation of data for ",
                         nCell,
                         " Cellules, replicate #",
                         rep))
            leduc <- leduc2022Generate(base, nCell)
            print(paste0("Starting benchmarking for ",
                         nCell,
                         " Cellules, replicate #",
                         rep))
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
                #leduc <- leducPCA(leduc)
                write.table(data.frame(nCell = as.integer(nCell),
                                     rep = as.integer(rep),
                                     size = object.size(leduc)),
                          file = sizeOutputPath,
                          append = TRUE,
                          col.names = FALSE,
                          row.names = FALSE)
            )
            write.csv(res, file = file.path(memoryUsageOutputDir,
                                            paste0(paste(
                                                nCell, rep, sep = "_"),".csv")))
            gc()
        }
    }
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
