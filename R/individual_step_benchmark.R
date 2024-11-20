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

benchmarkFilterSamples <- function(qfeatures) {
    subsetByColData(qfeatures, qfeatures$filterBench > 1)
}

benchmarkZeroisNA <- function(qfeatures) {
    zeroIsNA(qfeatures, i = seq_along(qfeatures))
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
write.table(data.frame(nCell = integer(),
                       rep = integer(),
                       step = character(),
                       sizeBefore = integer(),
                       sizeAfter = integer()),
            file = file.path(destPath, "qfeatures_size_report.tsv"),
            append = FALSE,
            col.names = TRUE,
            row.names = FALSE)

leduc <- scpdata::leduc2022_pSCoPE()

sizes <- c(500, 1000)

for (i in 1:replicate) {
    for (size in sizes) {
        qfeatures <- generateTMTPSM(leduc, size)
        qfeaturesSizeBefore <- object.size(qfeatures)
        for (stepPSM in names(stepsPSM)) {
            destFile <- file.path(destPath,
                paste0(size, "_", stepPSM, "_", i, ".csv"))
            write.csv(peakRAM(qfeaturesAfter <- stepsPSM[[stepPSM]](qfeatures),
                              write.table(data.frame(nCell = as.integer(size),
                                                     rep = as.integer(i),
                                                     step = stepPSM,
                                                     sizeBefore = qfeaturesSizeBefore,
                                                     sizeAfter = object.size(qfeaturesAfter)),
                                          file = file.path(destPath, "qfeatures_size_report.tsv"),
                                          append = TRUE,
                                          col.names = FALSE,
                                          row.names = FALSE)
                              ),
                destFile)


        }
        # qfeatures <- generateTMTPeptides(qfeatures)
        # for (stepPep in stepsPep) {
        #     peakRAM(stepPep(qfeatures))
        # }
        # qfeatures  <- generateTMTProteins(qfeatures)
        # for (stepPro in stepsPro) {
        #     peakRAM(stepPro(qfeatures))
        # }

    }
}
