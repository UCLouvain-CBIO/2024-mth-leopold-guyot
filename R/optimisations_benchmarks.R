
leducBenchWrapper <- function(git, sizes, replicate,
                              outputDir, subsetRowData = FALSE) {
    remotes::install_github(git, force = TRUE)
    source("R/optimisations_benchmarks_helpers.R")
    leduc2022Benchmark(sizes, replicate, outputDir, subsetRowData)
}

scpBase <- "UCLouvain-CBIO/scp@ed48ba714a4a354658c67345da1dcffd408cbe8a" # 1.15.2
qfeaturesBase <- "rformassspectrometry/QFeatures@1349182ca75f0cf209ba3ba564325fb6b31d4451" # 1.17.1
qfeaturesAgg <- "rformassspectrometry/QFeatures@b0f6cffdf5697f2871707e9fe74c30ca2b2a18e1" # 1.17.4
MultiAssayExperimentBase <- "waldronlab/MultiAssayExperiment@325c85ca26582027bcdd643df673a06ae649a36b" # 1.33.1
MultiAssayExperimentSubset <- "waldronlab/MultiAssayExperiment@7d6b9c4290259bfd2b59ce81f3ee8be3a447f618" # 1.33.2
git <- list("base" = c(MultiAssayExperimentBase, qfeaturesBase, scpBase),
         "agg" = c(MultiAssayExperimentBase, qfeaturesAgg, scpBase),
         "subset" = c(MultiAssayExperimentSubset, qfeaturesBase, scpBase),
         "all" = c(MultiAssayExperimentSubset, qfeaturesAgg, scpBase)
         )

callr::r(
    func = leducBenchWrapper,
    args = list(
        git = git[["base"]],
        sizes = 2000,
        rep = 1L,
        outputDir = paste0("dataOutput/optimisationsBench/subsetRowData"),
        subsetRowData = TRUE
    )
)

for (verName in names(git)) {
    ver <- git[[verName]]
    callr::r(
        func = leducBenchWrapper,
        args = list(
            git = ver,
            sizes = 2000,
            rep = 1L,
            outputDir = paste0("dataOutput/optimisationsBench/", verName)
        )
    )
}

callr::r(
    func = leducBenchWrapper,
    args = list(
        git = git[["all"]],
        sizes = 2000,
        rep = 1L,
        outputDir = paste0("dataOutput/optimisationsBench/allSubset"),
        subsetRowData = TRUE
    )
)
