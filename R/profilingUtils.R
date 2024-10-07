library("progress")

profilingWrapper <- function(dataFolder, profOutPath, ...) {
    if (!dir.exists(profOutPath)) {
        dir.create(profOutPath)
    }

    # Get the list of folders in the dataFolder
    folders <- list.dirs(dataFolder, full.names = TRUE, recursive = FALSE)
    cat("Starting Workflow Profiling ...")
    pb <- progress_bar$new(
        format = paste0(":spin Workflow Profiling: (:percent)",
                        " [:bar] :current/:total | Elapsed: :elapsed"),
        total = length(folders),
        clear = FALSE,              # Keep the progress output in the console
        width = 80
    )
    pb$tick(0)
    for (i in seq_along(folders)) {
        folder <- folders[i]
        designPath <- file.path(folder, "design.csv")
        quantPath <- file.path(folder, "quant.csv")
        if (file.exists(designPath) && file.exists(quantPath)) {
            profFile <- file.path(profOutPath,
                                  paste0(basename(folder), "_Rprof.out"))
            # Start profiling and direct output to the specified file
            Rprof(profFile)
            profilingWorkflow(quantPath,
                              designPath,
                              runCol = "run",
                              type = "TMT")
            Rprof(NULL)
        } else {
            cat("\nMissing design.csv or quant.csv in folder:", folder, "\n")
        }
        pb$tick()
    }
    cat("Workflow Profiling Finished.\n")
}

profilingWorkflow <- function(quantPath,
                              designPath,
                              runCol,
                              dataDIA = NULL,
                              type) {
    qfeatures <- stepDataImport(quantPath, designPath, runCol, dataDIA, type)
    assays <- 1:length(qfeatures)
    qfeatures <- stepPreProcessing(qfeatures, assays, type)
    qfeatures
}
