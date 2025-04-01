source(file.path("R", "minimumWorkflow.R"))
library(MultiAssayExperiment)
library(QFeatures)
library(peakRAM)
library(scp)
library(scpdata)

benchWrapper <- function(replicate, subsetRowData = FALSE) {
    leduc <- scpdata::leduc2022_pSCoPE()[,, 1:30]
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
    results <- do.call(rbind, lapply(1:replicate, function(i) {
        res <- peakRAM(
          leduc <- filterRow(leduc),
          leduc <- filterCol(leduc),
          aggregatePSM(leduc)
                       )
        res$replicate <- i
        return(res)
    }))

    return(results)
}
