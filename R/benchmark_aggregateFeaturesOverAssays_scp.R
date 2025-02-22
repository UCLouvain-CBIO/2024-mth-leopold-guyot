remotes::install_github("leopoldguyot/scp@aggregateFeaturesOverAssays_optimisation", force = TRUE)

library(QFeatures)
library(scp)
library(scpdata)
library(MsCoreUtils)
library(peakRAM)

source("R/generate_data.R")


leduc <- scpdata::leduc2022_pSCoPE()

aggPSM <- function(qfeatures, opti) {
    scp::aggregateFeaturesOverAssays(
        qfeatures,
        i = seq_along(qfeatures),
        "modseq",
        name = paste0("peptide_", seq_along(qfeatures)),
        fun = MsCoreUtils::robustSummary,
        uniformRowData = opti
    )
}

print(scp::aggregateFeaturesOverAssays)

leduc <- scpdata::leduc2022_pSCoPE()

replicate <- 3
sizes <- c(500, 1000, 2000, 4000)



for (size in sizes) {
    qfeat <- generateTMTPSM(leduc, size)
    for (rep in seq_len(replicate)) {
        suppressMessages(
            write.csv(
                peakRAM(
                    aggPSM(qfeat, FALSE),
                    aggPSM(qfeat, TRUE)
                ),
                file.path("dataOutput", "aggregateBenchSCP", paste0(size, "_", rep, ".csv"))
            )
        )
    }
}
