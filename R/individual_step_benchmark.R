library(scp)
library(scpdata)
library(peakRAM)
source(file.path("R", "generate_data.R"))

# Variables
destPath <- file.path("dataOutput", "individualStepsBenchmark")
replicate <- 2

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
consensusMapping <- function(qfeatures) {
    ## Generate a list of DataFrames with the information to modify
    rbindRowData(qfeatures, i = grep("^pep", names(qfeatures))) %>%
        data.frame %>%
        group_by(modseq) %>%
        ## The majority vote happens here
        mutate(Leading.razor.protein.symbol =
                   names(sort(table(Leading.razor.protein),
                              decreasing = TRUE))[1]) %>%
        select(modseq, Leading.razor.protein.symbol) %>%
        filter(!duplicated(modseq, Leading.razor.protein.symbol)) ->
        ppMap
    consensus <- lapply(peptideAssays, function(i) {
        ind <- match(rowData(qfeatures[[i]])$modseq, ppMap$modseq)
        DataFrame(Leading.razor.protein.symbol =
                      ppMap$Leading.razor.protein.symbol[ind])
    })
    ## Name the list
    names(consensus) <- peptideAssays
    ## Modify the rowData
    rowData(qfeatures) <- consensus

    qfeatures
}
benchmarkJoinPSM <- function(qfeatures) {
    joinAssays(qfeatures,
               i = paste0("peptide_", 1:(length(qfeatures)/2)),
               name = "peptides")
}

stepsPSM <- list(
    "benchmarkFilterFeatures" = benchmarkFilterFeatures,
    "benchmarkFilterSamples" = benchmarkFilterSamples,
    "benchmarkZeroisNA" = benchmarkZeroisNA,
    "benchmarkAggPSM" = benchmarkAggPSM,
    "benchmarkJoinPSM" = benchmarkJoinPSM
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

sizes <- c(500, 550)

for (i in 1:replicate) {
    for (size in sizes) {
        qfeatures <- generateTMTPSM(leduc, size)
        qfeaturesSizeBefore <- object.size(qfeatures)
        for (stepPSM in names(stepsPSM)) {
            destFile <- file.path(destPath,
                paste0(size, "_", stepPSM, "_", i, ".csv"))
            if (stepPSM == "benchmarkJoinPSM") {
                qfeatures_step <- benchmarkAggPSM(qfeatures)} else {
                    qfeatures_step <- qfeatures
                }
            write.csv(peakRAM(qfeaturesAfter <- stepsPSM[[stepPSM]](qfeatures_step),
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
