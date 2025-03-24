library(VerR)
MultiAssayExperimentBase <- "waldronlab/MultiAssayExperiment@325c85ca26582027bcdd643df673a06ae649a36b"
subsetByColDataOpti <- "waldronlab/MultiAssayExperiment@a293ab1fd1d5ac6671a4fb7a6cc91d14b21d0758"
QFeaturesBase <- "rformassspectrometry/QFeatures@dc229a06c4feb7f61d1bc1ad7c9d626d58deb93a"
QFeaturesAgg <- "rformassspectrometry/QFeatures@9299c7bf86b9598a4db49e04f76e999cdfdd4952"
scpBase <- "UCLouvain-CBIO/scp@8818589510ffd69673b533f614fc3f3eb7f9598d"

envCreate("subsetByColData",
          packages = c(subsetByColDataOpti,
                       QFeaturesBase,
                       scpBase))

envCreate("base",
          packages = c(MultiAssayExperimentBase,
                       QFeaturesBase,
                       scpBase))

envCreate("aggregation",
          packages = c(MultiAssayExperimentBase,
                       QFeaturesAgg,
                       scpBase))

envCopyTo(file.path("R", "scriptVerR.R"),
          targetPath = "R")
envCopyTo(file.path("R", "minimumWorkflow.R"),
          targetPath = "R")

write.csv(runInEnv({
    source(file.path("R", "scriptVerR.R"))
    benchWrapper(replicate = 1, subsetRowData = FALSE)
}), file = "dataOutput/optimisationsBench/noSubsetRowData.csv")

write.csv(runInEnv({
    source(file.path("R", "scriptVerR.R"))
    benchWrapper(replicate = 1, subsetRowData = TRUE)
}), file = "dataOutput/optimisationsBench/subsetRowData.csv")
