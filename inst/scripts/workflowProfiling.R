################################################################################
############################# WORKFLOW PROFILING ###############################
################################################################################

############################# Packages loading #################################

cat("Loading packages...\n")
suppressMessages({
    library("QFeatures", verbose = FALSE)
    library("scp", verbose = FALSE)
})
############################### Sourcing #######################################

cat("Sourcing R files...\n")
suppressMessages({
    source(file.path("R", "workflowStepOneDataImport.R"), verbose = FALSE)
    source(file.path("R", "workflowStepTwoPreProcessing.R"), verbose = FALSE)
    source(file.path("R", "workflowStepAggregation.R"), verbose = FALSE)

    source(file.path("R", "generateData.R"), verbose = FALSE)
    source(file.path("R", "profilingUtils.R"), verbose = FALSE)
})

################################# Parameters ###################################
cat("\nARGUMENTS:\n")
args <- commandArgs(trailingOnly = TRUE)

params <- if (length(args) == 0) {
    list()
} else {
    eval(parse(text = args[[2]]))
}


if (!is.null(params$naRateRange)) {
  naRateRange <- params$naRateRange
  cat("naRateRange: c(", paste(naRateRange, collapse = ", "), ")\n")
} else {
  cat("No naRateRange argument provided. Using c(0.4, 0.6) as default.\n")
  naRateRange <- c(0.4, 0.6)
}

if (!is.null(params$nFeaturesRange)) {
  nFeaturesRange <- params$nFeaturesRange
  cat("nFeaturesRange: c(", paste(nFeaturesRange, collapse = ", "), ")\n")
} else {
  cat("No nFeaturesRange argument provided. Using c(500, 1000) as default.\n")
  nFeaturesRange <- c(500, 1000)
}

if (!is.null(params$nRunsRange)) {
  nRunsRange <- params$nRunsRange
  cat("nRunsRange: c(", paste(nRunsRange, collapse = ", "), ")\n")
} else {
  cat("No nRunsRange argument provided. Using c(4, 8) as default.\n")
  nRunsRange <- c(4, 8)
}

if (!is.null(params$nReplicates)) {
  cat("nReplicates:", params$nReplicates, "\n")
  nReplicates <- params$nReplicates
} else {
  cat("No nReplicates argument provided. Using 1 as default.\n")
  nReplicates <- 1
}

################################## Main ########################################

main <- function() {
    # Ensure temp files are cleaned on exit
    on.exit({
        unlink("profilingTemp", recursive = TRUE)
        cat("Temporary files cleaned up.\n")
    })

    ############################# Data Generation ##############################

    ### TMT ###

    replicateData(
        nFeaturesRange = nFeaturesRange,
        nRunsRange = nRunsRange,
        naRateRange = naRateRange,
        nReplicates = nReplicates,
        type = "TMT",
        folder = "profilingTemp"
    )

    ### label free DIA-NN ###

    ### plexed DIA-NN ###

    ############################# Profiling ####################################
    profOutPath <- "dataOutput/profilingResults"

    if (dir.exists(profOutPath)) {
        unlink(profOutPath, recursive = TRUE)
    }
    profilingWrapper(
        dataFolder = "profilingTemp",
        profOutPath = profOutPath
    )
}

main()
