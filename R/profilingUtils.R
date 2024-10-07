profilingWrapper <- function(dataFolder, profOutPath, ...) {
    if (!dir.exists(profOutPath)) {
        dir.create(profOutPath)
    }

    # Get the list of folders in the dataFolder
    folders <- list.dirs(dataFolder, full.names = TRUE, recursive = FALSE)

    total_folders <- length(folders)
    pb <- txtProgressBar(min = 0, max = total_folders, style = 3)

    for (i in seq_along(folders)) {
        folder <- folders[i]
        designPath <- file.path(folder, "design.csv")
        quantPath <- file.path(folder, "quant.csv")
        cat("\nProfiling folder: ", folder, "\n")
        if (file.exists(designPath) && file.exists(quantPath)) {
            profFile <- file.path(profOutPath, paste0(basename(folder), "_Rprof.out"))

            Rprof(profFile) # Start profiling and direct output to the specified file
            profilingWorkflow(quantPath, designPath, runCol = "run", type = "TMT")
            Rprof(NULL)
        } else {
            cat("\nMissing design.csv or quant.csv in folder:", folder, "\n")
        }

        cat("\nFolder ", folder, "Profiled.\n")
        setTxtProgressBar(pb, i)
    }

    close(pb)
}

profilingWorkflow <- function(quantPath, designPath, runCol, dataDIA = NULL, type) {
    qfeatures <- stepDataImport(quantPath, designPath, runCol, dataDIA, type)
    assays <- 1:length(qfeatures)
    qfeatures <- stepPreProcessing(qfeatures, assays, type)
    qfeatures
}
