#' Execute Data Generation and Profiling Workflow
#'
#' This function generates synthetic datasets for profiling based on specified parameters
#' and then performs profiling across these datasets. Temporary files generated during the process
#' are automatically cleaned up on exit.
#'
#' @param nFeaturesRange Integer vector. Specifies the range of the number of features to generate. Defaults to `c(100, 1000)`.
#' @param nRunsRange Integer vector. Specifies the range of the number of runs to generate. Defaults to `c(4, 8)`.
#' @param naRateRange Numeric vector. Specifies the range of the rate of missing values to introduce. Defaults to `c(0.4, 0.6)`.
#' @param nReplicates Integer. Number of replicate datasets to generate for each combination. Defaults to `1`.
#' @param profilingTemp Character. Path to a temporary folder for storing intermediate profiling data. Defaults to `"profilingTemp"`.
#' @param outputFolder Character. Path where profiling results will be saved. Defaults to `"dataOutput/profilingResults"`.
#'
#' @return No return value. This function performs data generation and profiling, saving results to `outputFolder`.
#'
#' @details The function first generates synthetic data for TMT profiling by calling `replicateData`.
#'          Label-free DIA and plexed DIA generation functions can be added to extend functionality.
#'          Next, `profilingWrapper` is called to execute profiling on the generated data.
#'          Temporary files are cleaned up upon completion or if an error occurs.
#'
#' @examples
#' # Run profiling execution with default parameters:
#' # profilingExecution()
#'
#' # Customize parameters:
#' # profilingExecution(nFeaturesRange = c(200, 500), nRunsRange = c(6, 12), naRateRange = c(0.3, 0.5))
#'
#' @export
profilingExecution <- function(nFeaturesRange = c(500, 1000),
                               nRunsRange = c(4, 8),
                               naRateRange = c(0.4, 0.6),
                               nReplicates = 1,
                               profilingTemp = "profilingTemp",
                               outputFolder = "dataOutput/profilingResults") {
  # Ensure temp files are cleaned on exit
  on.exit({
    Rprof(NULL)
    unlink(profilingTemp, recursive = TRUE)
    cat("Temporary files cleaned up.\n")
  })

  ### Data Generation ###
  ### TMT ###

  replicateData(
    nFeaturesRange = nFeaturesRange,
    nRunsRange = nRunsRange,
    naRateRange = naRateRange,
    nReplicates = nReplicates,
    type = "TMT",
    folder = profilingTemp
  )

  ### label free DIA-NN ###

  ### plexed DIA-NN ###

  #### Profiling ####

  if (dir.exists(outputFolder)) {
    unlink(outputFolder, recursive = TRUE)
  }
  profilingWrapper(
    dataFolder = profilingTemp,
    outputFolder = outputFolder
  )
}

#' Profile a Workflow Across Multiple Datasets
#'
#' This function profiles a workflow over multiple datasets by iterating through
#' folders containing quantitative and design CSV files. Profiling results are saved
#' for each dataset into separate output files.
#'
#' @param dataFolder Character. Path to the folder containing subfolders with data files.
#'        Each subfolder should contain `quant.csv` and `design.csv`.
#' @param outputFolder Character. Path where profiling output files will be saved.
#' @param ... Additional arguments passed to the `profilingWorkflow` function.
#'
#' @return No return value. This function performs profiling and saves the output files.
#'
#' @details This function iterates through each subfolder in `dataFolder`, checks for the presence of
#'          `design.csv` and `quant.csv`, and then performs profiling of the specified workflow.
#'          Profiling data is saved to files in `outputFolder`, named after each dataset folder.
#'
#' @examples
#' # To profile workflows across multiple datasets:
#' # profilingWrapper(dataFolder = "data", outputFolder = "profile_results")
#'
#' @import progress
#' @importFrom utils Rprof
#'
#' @export
profilingWrapper <- function(dataFolder, outputFolder, ...) {
    if (!dir.exists(outputFolder)) {
        dir.create(outputFolder, recursive = TRUE)
    }

    # Get the list of folders in the dataFolder
    folders <- list.dirs(dataFolder, full.names = TRUE, recursive = FALSE)
    cat("Starting Workflow Profiling ...")
    pb <- progress_bar$new(
        format = paste0(
            "Workflow Profiling: (:percent)",
            " [:bar] :current/:total | Elapsed: :elapsed"
        ),
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
            profFile <- file.path(
                outputFolder,
                paste0(basename(folder), "_Rprof.out")
            )
            # Start profiling and direct output to the specified file
            Rprof(profFile)
            profilingWorkflow(quantPath,
                designPath,
                runCol = "run",
                type = "TMT"
            )
            Rprof(NULL)
        } else {
            cat("\nMissing design.csv or quant.csv in folder:", folder, "\n")
        }
        pb$tick()
    }
    cat("Workflow Profiling Finished.\n")
}

#' Perform Workflow Steps for Quantitative Data Profiling
#'
#' This function performs several steps of a profiling workflow, including data import,
#' preprocessing, and aggregation. It is intended for use within `profilingWrapper`.
#'
#' @param quantPath Character. Path to the quantitative data CSV file.
#' @param designPath Character. Path to the experimental design data CSV file.
#' @param runCol Character. Column name in the design data identifying experimental runs.
#' @param dataDIA Optional. Additional DIA data if applicable; defaults to `NULL`.
#' @param type Character. Type of the dataset, typically "TMT" or other supported types.
#'
#' @return The `qfeatures` object after processing steps are applied.
#'
#' @details This function processes quantitative data by performing import, preprocessing,
#'          and aggregation steps. These steps are modular and allow the function to handle
#'          various data processing tasks.
#'
profilingWorkflow <- function(quantPath,
    designPath,
    runCol,
    dataDIA = NULL,
    type) {
    qfeatures <- stepDataImport(quantPath, designPath, runCol, dataDIA, type)
    assays <- 1:length(qfeatures)

    qfeatures <- stepPreProcessing(qfeatures, assays, type)

    qfeatures <- stepAggregate(qfeatures, assays, "PSM")
    assays <- grep("peptides_", names(qfeatures))
}
