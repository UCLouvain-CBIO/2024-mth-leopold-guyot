library(VerR)
library(dplyr)
MultiAssayExperimentBase <- "waldronlab/MultiAssayExperiment@325c85ca26582027bcdd643df673a06ae649a36b"
subsetByColDataOpti <- "waldronlab/MultiAssayExperiment@a293ab1fd1d5ac6671a4fb7a6cc91d14b21d0758"
QFeaturesBase <- "rformassspectrometry/QFeatures@dc229a06c4feb7f61d1bc1ad7c9d626d58deb93a"
QFeaturesAgg <- "leopoldguyot/QFeatures@333786236576d393d84415280ee563cd6d953469"
scpBase <- "UCLouvain-CBIO/scp@8818589510ffd69673b533f614fc3f3eb7f9598d"
scpOpti <- "leopoldguyot/scp@68de84b559f1e92b31d1190c1add678914012e87"

envCreate("base",
          packages = c(MultiAssayExperimentBase,
                       QFeaturesBase,
                       scpBase,
                       "bioc::scpdata",
                       "peakRAM"))

envCreate("subsetByColData",
          packages = c(subsetByColDataOpti,
                       QFeaturesBase,
                       scpBase,
                       "bioc::scpdata",
                       "peakRAM"))

envCreate("aggregation",
          packages = c(MultiAssayExperimentBase,
                       QFeaturesAgg,
                       scpBase,
                       "bioc::scpdata",
                       "peakRAM"))


envCreate("allOpti",
          packages = c(subsetByColDataOpti,
                       QFeaturesAgg,
                       scpBase,
                       "bioc::scpdata",
                       "peakRAM"))


envCopyTo(file.path("R", "scriptVerR.R"),
          targetPath = file.path("R", "scriptVerR.R")
)

envCopyTo(file.path("R", "minimumWorkflow.R"),
          targetPath = file.path("R", "minimumWorkflow.R"))

noSub <- runInEnv({
    source(file.path("R", "scriptVerR.R"))
    benchWrapper(replicate = 1, nCell = c(4000), subsetRowData = FALSE, SE = FALSE)
})
noSub <- bind_rows(noSub, .id = "version")

write.csv(noSub, file = "dataOutput/optimisationsBench/noSubsetRowDataSCE.csv")

sub <- runInEnv({
    source(file.path("R", "scriptVerR.R"))
    benchWrapper(replicate = 1, nCell = c(4000), subsetRowData = TRUE, SE = FALSE)
})

sub <- bind_rows(sub, .id = "version")

write.csv(sub, file = "dataOutput/optimisationsBench/subsetRowDataSCE.csv")

sub <- runInEnv({
    source(file.path("R", "scriptVerR.R"))
    benchWrapper(replicate = 1, nCell = c(4000), subsetRowData = FALSE, SE = TRUE)
})

sub <- bind_rows(sub, .id = "version")

write.csv(sub, file = "dataOutput/optimisationsBench/noSubsetRowDataSE.csv")

sub <- runInEnv({
    source(file.path("R", "scriptVerR.R"))
    benchWrapper(replicate = 1, nCell = c(4000), subsetRowData = TRUE, SE = TRUE)
})

sub <- bind_rows(sub, .id = "version")

write.csv(sub, file = "dataOutput/optimisationsBench/subsetRowDataSE.csv")

unlink(".envs", recursive = TRUE, force = TRUE)
