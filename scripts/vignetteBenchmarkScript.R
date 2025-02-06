#!/usr/bin/env Rscript


message("Sourcing ...")
suppressMessages({
    library("optparse")
    source(file.path("R", "vignette_benchmark.R"))
    source(file.path("R", "utils.R"))
})

option_list <- list(
    make_option(c("-c", "--cellRange"),
        type = "character", default = "c(1000, 2000, 4000)",
        help = 'Integer vector of the different size of dataset to use. Should be passed as a character of the form "c()"',
        metavar = "character"
    ),
    make_option(c("-r", "--nReplicates"),
        type = "integer", default = 3,
        help = "Number of replicate by dataset", metavar = "integer"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

cellRange <- eval(parse(text = opt$cellRange))
nReplicates <- opt$nReplicates

print("Starting benchmark...")
leduc2022Benchmark(cellRange, nReplicates)

print("Rendering Vignette...")
renderRmarkdown("leduc2022Results.rmd")
