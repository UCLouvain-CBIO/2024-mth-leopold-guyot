stepPreProcessing <- function(qfeatures, assays, type) {
    switch(type,
        "TMT" = stepPreProcessingTMT(qfeatures, assays),
        "LFDIA" = stepPreProcessingLFDIA(qfeatures, assays),
        "plexDIA" = stepPreProcessingPlexDIA(qfeatures, assays),
        stop("Invalid type")
    )
}


stepPreProcessingTMT <- function(qfeatures, assays) {
    qfeatures <- zeroIsNA(qfeatures, i = assays)

    qfeatures <- .filterPSMTMT(qfeatures, assays)
    qfeatures <- .filterAssayTMT(qfeatures)
    qfeatures <- .filterSCRTMT(qfeatures, assays)
    qfeatures <- .filterFDRTMT(qfeatures, assays)
}

.filterPSMTMT <- function(qfeatures, assays) {
    filterFeatures(
        qfeatures,
        ~ Reverse != "+" &
            Potential.contaminant != "+" &
            !is.na(PIF) & PIF > 0.8,
        i = assays
    )
}

.filterAssayTMT <- function(qfeatures) {
    keepAssay <- dims(qfeatures)[1, ] > 120
    qfeatures[, , keepAssay]
}

.filterSCRTMT <- function(qfeatures, assays) {
    qfeatures <- computeSCR(qfeatures,
        i = assays,
        colvar = "SampleType",
        carrierPattern = "Carrier",
        samplePattern = "Macrophage|Monocyte",
        sampleFUN = "mean",
        rowDataName = "MeanSCR"
    )
    qfeatures <- filterFeatures(qfeatures,
        ~ !is.na(MeanSCR) &
            MeanSCR < 0.1,
        i = assays
    )
}

.filterFDRTMT <- function(qfeatures, assays) {
    qfeatures <- pep2qvalue(qfeatures,
                     i = assays,
                     PEP = "PEP",
                     groupBy = "Leading.razor.protein",
                     rowDataName = "qvalueProteins")
    filterFeatures(qfeatures,
                          ~ qvalueProteins < 0.01)
}

.divideReferenceTMT <- function(qfeatures, assays) {
    qfeatures <- divideByReference(qfeatures,
                             i = assays,
                             colvar = "SampleType",
                             samplePattern = ".",
                             refPattern = "Reference")
}
