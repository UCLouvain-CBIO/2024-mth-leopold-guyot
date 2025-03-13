source(file = "R/minimumWorkflow.R")
leducBenchWrapper <- function(git, replicate,
                              outputDir, subsetRowData = FALSE) {
    remotes::install_github(git, force = TRUE)
    library(MultiAssayExperiment)
    library(QFeatures)
    library(microbenchmark)
    source(file = "R/minimumWorkflow.R")
    leduc <- scpdata::leduc2022_pSCoPE()[,, 1:134]
    if (subsetRowData) {
        for (assay in seq_along(leduc)) {
        rowData(leduc[[assay]]) <- rowData(leduc[[assay]])[, c("dart_PEP",
                                                        "PIF",
                                                        "Proteins",
                                                        "Leading.razor.protein",
                                                        "dart_qval",
                                                        "Potential.contaminant",
                                                        "Reverse",
                                                        "modseq")]
        }
    }
    write.csv(microbenchmark::microbenchmark(minimalWorkflow(leduc), times = replicate),
        file = paste0(outputDir, "/", "microbenchres.csv"))
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
cat("Starting subset config benchmark\n")
callr::r(
    func = leducBenchWrapper,
    args = list(
        git = git[["base"]],
        rep = 3L,
        outputDir = paste0("dataOutput/optimisationsBench/subsetRowData"),
        subsetRowData = TRUE
    )
)

for (verName in names(git)) {
    ver <- git[[verName]]
    cat("Starting", verName, "config benchmark\n")
    callr::r(
        func = leducBenchWrapper,
        args = list(
            git = ver,
            rep = 3L,
            outputDir = paste0("dataOutput/optimisationsBench/", verName)
        )
    )
}

cat("Starting allSubset config benchmark\n")
callr::r(
    func = leducBenchWrapper,
    args = list(
        git = git[["all"]],
        rep = 3L,
        outputDir = paste0("dataOutput/optimisationsBench/allSubset"),
        subsetRowData = TRUE
    )
)
