library(scp)
library(scpdata)
library(peakRAM)
library(pcaMethods)
source(file.path("R", "generate_data.R"))
source(file.path("R", "utils.R"))

# Variables
destPath <- file.path("dataOutput", "individualStepsBenchmark")
replicate <- 3
sizes <- c(500, 1000, 2000, 4000)


benchmarkFilterFeatures <- function(qfeatures) {
    filterFeatures(qfeatures, ~ filterBench > 1)
}

benchmarkFilterSamples <- function(qfeatures) {
    subsetByColData(qfeatures, qfeatures$filterBench > 1)
}

benchmarkZeroisNA <- function(qfeatures) {
    zeroIsNA(qfeatures, i = seq_along(qfeatures))
}

benchmarkAggPSM <- function(qfeatures) {
    aggregateFeaturesOverAssays(
        qfeatures,
        i = seq_along(qfeatures),
        "modseq",
        name = paste0("peptide_", seq_along(qfeatures)),
        fun = MsCoreUtils::robustSummary
    )
}

benchmarkJoinPSM <- function(qfeatures) {
    joinAssays(qfeatures,
        i = paste0("peptide_", 1:(length(qfeatures) / 2)),
        name = "peptides"
    )
}

benchmarkNormSampPep <- function(qfeatures, assayName = "peptides") {
    sweep(qfeatures,
        i = assayName,
        MARGIN = 2,
        FUN = "/",
        STATS = colMedians(assay(qfeatures[[assayName]]), na.rm = TRUE),
        name = paste0(assayName, "_norm")
    )
}

benchmarkNormFeatPep <- function(qfeatures, assayName = "peptides") {
    sweep(qfeatures,
        i = assayName,
        MARGIN = 1,
        FUN = "/",
        STATS = colMedians(assay(qfeatures[[assayName]]), na.rm = TRUE),
        name = paste0(assayName, "_norm")
    )
}

benchmarkLogPep <- function(qfeatures, assayName = "peptides") {
    logTransform(qfeatures,
        base = 2,
        i = assayName,
        name = paste0(assayName, "_log")
    )
}

benchmarkAggPep <- function(qfeatures, assayName = "peptides") {
    aggregateFeatures(qfeatures,
        i = assayName,
        fcol = "Leading.razor.protein.symbol",
        fun = colMedians, na.rm = TRUE,
        name = sub("peptides", "proteins", assayName)
    )
}

benchmarkImputPro <- function(qfeatures, assayName = "proteins") {
    impute(qfeatures,
        i = assayName,
        method = "knn",
        k = 3, rowmax = 1, colmax = 1,
        name = "imputed_proteins"
    )
}

benchmarkPCAPro <- function(qfeatures, assayName = "proteins") {
    nipals_res <- assay(qfeatures[[assayName]]) %>%
        as.data.frame() %>%
        mutate_all(~ ifelse(is.nan(.), NA, .)) %>%
        t() %>%
        pcaMethods::pca(method = "nipals", nPcs = 20)
    reducedDim(qfeatures[[assayName]], "NIPALS") <- pcaMethods::scores(nipals_res)
    qfeatures
}

stepsPSM <- list(
    "benchmarkFilterFeatures" = benchmarkFilterFeatures,
    "benchmarkFilterSamples" = benchmarkFilterSamples,
    "benchmarkZeroisNA" = benchmarkZeroisNA,
    "benchmarkAggPSM" = benchmarkAggPSM,
    "benchmarkJoinPSM" = benchmarkJoinPSM
)

stepsPep <- list(
    "benchmarkNormSampPep" = benchmarkNormSampPep,
    "benchmarkNormFeatPep" = benchmarkNormFeatPep,
    "benchmarkLogPep" = benchmarkLogPep,
    "benchmarkAggPep" = benchmarkAggPep
)

stepsPro <- list(
    "benchmarkImputePro" = benchmarkImputPro # ,
    # "benchmarkPCAPro" = benchmarkPCAPro # not useful
)

unlink(destPath, recursive = TRUE)
dir.create(destPath)
write.table(
    data.frame(
        nCell = integer(),
        rep = integer(),
        step = character(),
        sizeTotalBefore = numeric(),
        sizeTotalAfter = numeric(),
        sizeAssayBefore = numeric(),
        sizeAssayAfter = numeric(),
        sizeRowDataBefore = numeric(),
        sizeRowDataAfter = numeric(),
        sizeColDataBefore = numeric(),
        sizeColDataAfter = numeric()
    ),
    file = file.path(destPath, "qfeatures_size_report.tsv"),
    append = FALSE,
    col.names = TRUE,
    row.names = FALSE
)

leduc <- scpdata::leduc2022_pSCoPE()


for (i in 1:replicate) {
    for (size in sizes) {
        print(paste0("Starting dataset of size: ", size, " replicate no: ", i))
        qfeatures <- generateTMTPSM(leduc, size)
        qfeaturesTotalSizeBefore <- object.size(qfeatures)
        qfeaturesAssaySizeBefore <- getAssaySize(qfeatures)
        qfeaturesRowDataSizeBefore <- getRowDataSize(qfeatures)
        qfeaturesColDataSizeBefore <- getColDataSize(qfeatures)
        for (stepPSM in names(stepsPSM)) {
            destFile <- file.path(
                destPath,
                paste0(size, "_", stepPSM, "_", i, ".csv")
            )
            if (stepPSM == "benchmarkJoinPSM") {
                qfeatures_step <- benchmarkAggPSM(qfeatures)
            } else {
                qfeatures_step <- qfeatures
            }
            suppressMessages(
                write.csv(
                    peakRAM(
                        qfeaturesAfter <- stepsPSM[[stepPSM]](qfeatures_step),
                        write.table(
                            data.frame(
                                nCell = as.integer(size),
                                rep = as.integer(i),
                                step = stepPSM,
                                sizeTotalBefore = qfeaturesTotalSizeBefore,
                                sizeTotalAfter = object.size(qfeaturesAfter),
                                sizeAssayBefore = qfeaturesAssaySizeBefore,
                                sizeAssayAfter = getAssaySize(qfeaturesAfter),
                                sizeRowDataBefore = qfeaturesRowDataSizeBefore,
                                sizeRowDataAfter = getRowDataSize(qfeaturesAfter),
                                sizeColDataBefore = qfeaturesColDataSizeBefore,
                                sizeColDataAfter = getColDataSize(qfeaturesAfter)
                            ),
                            file = file.path(destPath, "qfeatures_size_report.tsv"),
                            append = TRUE,
                            col.names = FALSE,
                            row.names = FALSE
                        )
                    ),
                    destFile
                )
            )
        }
        qfeatures <- generateTMTPeptides(qfeatures)
        qfeaturesTotalSizeBefore <- object.size(qfeatures)
        qfeaturesAssaySizeBefore <- getAssaySize(qfeatures)
        qfeaturesRowDataSizeBefore <- getRowDataSize(qfeatures)
        qfeaturesColDataSizeBefore <- getColDataSize(qfeatures)
        for (stepPep in names(stepsPep)) {
            destFile <- file.path(
                destPath,
                paste0(size, "_", stepPep, "_", i, ".csv")
            )
            suppressMessages(
                write.csv(
                    peakRAM(
                        qfeaturesAfter <- stepsPep[[stepPep]](qfeatures),
                        write.table(
                            data.frame(
                                nCell = as.integer(size),
                                rep = as.integer(i),
                                step = stepPep,
                                sizeTotalBefore = qfeaturesTotalSizeBefore,
                                sizeTotalAfter = object.size(qfeaturesAfter),
                                sizeAssayBefore = qfeaturesAssaySizeBefore,
                                sizeAssayAfter = getAssaySize(qfeaturesAfter),
                                sizeRowDataBefore = qfeaturesRowDataSizeBefore,
                                sizeRowDataAfter = getRowDataSize(qfeaturesAfter),
                                sizeColDataBefore = qfeaturesColDataSizeBefore,
                                sizeColDataAfter = getColDataSize(qfeaturesAfter)
                            ),
                            file = file.path(destPath, "qfeatures_size_report.tsv"),
                            append = TRUE,
                            col.names = FALSE,
                            row.names = FALSE
                        )
                    ),
                    destFile
                )
            )
        }
        qfeatures <- generateTMTProteins(qfeatures)
        qfeaturesTotalSizeBefore <- object.size(qfeatures)
        qfeaturesAssaySizeBefore <- getAssaySize(qfeatures)
        qfeaturesRowDataSizeBefore <- getRowDataSize(qfeatures)
        qfeaturesColDataSizeBefore <- getColDataSize(qfeatures)
        for (stepPro in names(stepsPro)) {
            destFile <- file.path(
                destPath,
                paste0(size, "_", stepPro, "_", i, ".csv")
            )
            suppressMessages(
                write.csv(
                    peakRAM(
                        qfeaturesAfter <- stepsPro[[stepPro]](qfeatures),
                        write.table(
                            data.frame(
                                nCell = as.integer(size),
                                rep = as.integer(i),
                                step = stepPro,
                                sizeTotalBefore = qfeaturesTotalSizeBefore,
                                sizeTotalAfter = object.size(qfeaturesAfter),
                                sizeAssayBefore = qfeaturesAssaySizeBefore,
                                sizeAssayAfter = getAssaySize(qfeaturesAfter),
                                sizeRowDataBefore = qfeaturesRowDataSizeBefore,
                                sizeRowDataAfter = getRowDataSize(qfeaturesAfter),
                                sizeColDataBefore = qfeaturesColDataSizeBefore,
                                sizeColDataAfter = getColDataSize(qfeaturesAfter)
                            ),
                            file = file.path(destPath, "qfeatures_size_report.tsv"),
                            append = TRUE,
                            col.names = FALSE,
                            row.names = FALSE
                        )
                    ),
                    destFile
                )
            )
        }
    }
}
