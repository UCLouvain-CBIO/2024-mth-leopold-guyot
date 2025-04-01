minimalWorkflow <- function(qfeatures) {
    qfeatures <- filterRow(qfeatures)
    qfeatures <- filterCol(qfeatures)
    qfeatures <- aggregatePSM(qfeatures)
}

filterRow <- function(qfeatures) {
    filterFeatures(qfeatures, ~ Potential.contaminant != "+" &
                       !grepl("CON", Proteins) &
                       Reverse != "+" &
                       !grepl("REV", Leading.razor.protein) &
                       (is.na(PIF) | PIF > 0.8) &
                       dart_qval < 0.001)
}

filterCol <- function(qfeatures) {
    qfeatures[, colData(qfeatures)$SampleType != "Unused", ]
}

aggregatePSM <- function(qfeatures) {
    peptideAssays <- paste0("peptides_", names(qfeatures))


    suppressMessages(aggregateFeatures(qfeatures,
                               i = names(qfeatures),
                               fcol = "modseq",
                               name = peptideAssays,
                               fun = colMeans
    ))
}

aggregatePSMscp <- function(qfeatures) {
    peptideAssays <- paste0("peptides_", names(qfeatures))


    aggregateFeaturesOverAssays(qfeatures,
                               i = names(qfeatures),
                               fcol = "modseq",
                               name = peptideAssays,
                               fun = colMeans
    )
}
