profilingWrapper <- function(dataFolder, profOutPath, ...) {
    cat("test",profOutPath)
    if (!dir.exists(profOutPath)) {
        dir.create(profOutPath)
    }

    # Get the list of folders in the dataFolder
    folders <- list.dirs(dataFolder, full.names = TRUE, recursive = FALSE)

    for (folder in folders) {
        designPath <- file.path(folder, "design.csv")
        quantPath <- file.path(folder, "quant.csv")

        if (file.exists(designPath) && file.exists(quantPath)) {
            profFile <- file.path(profOutPath, paste0(basename(folder), "_Rprof.out"))

            Rprof(profFile) # Start profiling and direct output to the specified file
            profilingWorkflow(quantPath, designPath, runCol = "run", type = "TMT")
            Rprof(NULL)

            cat("Profiled folder:", folder, "\n")
        } else {
            cat("Missing design.csv or quant.csv in folder:", folder, "\n")
        }
    }
}

profilingWorkflow <- function(quantPath, designPath, runCol, dataDIA = NULL, type) {
    qfeatures <- stepDataImport(quantPath, designPath, runCol, dataDIA, type)
    assays <- 1:length(qfeatures)
    qfeatures <- stepPreProcessing(qfeatures, assays, type)
    qfeatures
}
