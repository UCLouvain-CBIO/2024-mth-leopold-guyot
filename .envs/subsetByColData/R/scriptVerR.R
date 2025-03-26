source(file.path("R", "minimumWorkflow.R"))
library(MultiAssayExperiment)
library(QFeatures)
library(peakRAM)
library(scpdata)

benchWrapper <- function(replicate, subsetRowData = FALSE) {
    leduc <- scpdata::leduc2022_pSCoPE()[,, 1:125]
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
    peakRAM::peakRAM(leduc <- filterRow(leduc),
                     leduc <- filterCol(leduc),
                     aggregatePSM(leduc)
                     )
}
