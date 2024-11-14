library(scp)
library(scpdata)
library(peakRAM)
source(file.path("R", "generate_data.R"))

# Variables
destPath <- file.path("dataOutput", "individualStepsBenchmark")
replicate <- 3

benchmarkFilterFeatures <- function(qfeatures) {
    filterFeatures(qfeatures, ~ filterBench > 1)
}

benchmarkFilterSamples <- function(variables) {
    subsetByColData(qfeatures, leduc$filterBench > 1)
}

benchmarkZeroisNA <- function(variables) {
    zeroIsNA(qfeatures)
}

stepsPSM <- list(
    "benchmarkFilterFeatures" = benchmarkFilterFeatures,
    "benchmarkFilterSamples" = benchmarkFilterSamples,
    "benchmarkZeroisNA" = benchmarkZeroisNA
    )

stepsPep <- c()

stepsPro <- c()

unlink(destPath, recursive = TRUE)
dir.create(destPath)

leduc <- scpdata::leduc2022_pSCoPE()

sizes <- c(500, 1000)

for (i in 1:replicate) {
    for (size in sizes) {
        qfeatures <- generateTMTPSM(leduc, size)
        for (stepPSM in names(stepsPSM)) {
            destFile <- file.path(destPath,
                paste0(size, "_", stepPSM, "_", i))
            write.csv(peakRAM(stepsPSM[[stepPSM]](qfeatures)),
                destFile)
        }
        qfeatures <- generateTMTPeptides(qfeatures)
        for (stepPep in stepsPep) {
            peakRAM(stepPep(qfeatures))
        }
        qfeatures  <- generateTMTProteins(qfeatures)
        for (stepPro in stepsPro) {
            peakRAM(stepPro(qfeatures))
        }
    }
}
