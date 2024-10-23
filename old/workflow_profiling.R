#!/usr/bin/env Rscript

################################################################################
############################# WORKFLOW PROFILING ###############################
################################################################################

############################# Package loading #################################

message("Loading packages...")
suppressMessages({
    library("benchmarkQFeatures", verbose = FALSE)
})

message("ARGUMENTS:")
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

if (!is.null(params$outputFolder)) {
  cat("outputFolder:", params$outputFolder, "\n")
  outputFolder <- params$outputFolder
} else {
  cat(paste("No outputFolder argument provided.",
            "Using dataOutput/profilingResults folder as default.\n"))
  outputFolder <- "dataOutput/profilingResults"
}

################################## Profiling ###################################

profilingExecution(
  nFeaturesRange = nFeaturesRange,
  nRunsRange = nRunsRange,
  naRateRange = naRateRange,
  nReplicates = nReplicates,
  profilingTemp = "profilingTemp",
  outputFolder = outputFolder
)

