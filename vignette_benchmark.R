#!/usr/bin/env Rscript

################################################################################
############################# WORKFLOW PROFILING ###############################
################################################################################

############################# Scripts loading #################################

source("R/vignette_leduc2022_script.R")
source("R/vignette_benchmark.R")


start <- Sys.time()
leduc2022BenchmarkDetails(c(1000, 2000, 4000, 8000, 16000), 3)
end <- Sys.time()

print(start - end)
