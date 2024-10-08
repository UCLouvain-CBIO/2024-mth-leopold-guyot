library("progress", verbose = FALSE)

profilingWrapper <- function(dataFolder, profOutPath, ...) {
    if (!dir.exists(profOutPath)) {
        dir.create(profOutPath)
    }

    # Get the list of folders in the dataFolder
    folders <- list.dirs(dataFolder, full.names = TRUE, recursive = FALSE)
    cat("Starting Workflow Profiling ...")
    pb <- progress_bar$new(
        format = paste0("Workflow Profiling: (:percent)",
                        " [:bar] :current/:total | Elapsed: :elapsed"),
        total = length(folders),
        clear = FALSE,
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

    qfeatures <- stepAggregate(qfeatures, assays)
    assays <- grep("peptides_", names(qfeatures))
}

parseArgs <- function(args) {
  parsed_args <- list()
  i <- 1
  while (i <= length(args)) {
    # Check if the argument starts with '--'
    if (startsWith(args[i], "--")) {
      key <- sub("--", "", args[i])
      if (i + 1 <= length(args) && !startsWith(args[i + 1], "--")) {
        parsed_args[[key]] <- args[i + 1]
        i <- i + 1
      } else {
        parsed_args[[key]] <- TRUE
      }
    }
    i <- i + 1
  }
  return(parsed_args)
}