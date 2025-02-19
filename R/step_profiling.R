library(profvis)
library(QFeatures)
library(scp)
library(scpdata)

source("R/generate_data.R")

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

benchmarkImputPro <- function(qfeatures, assayName = "proteins") {
    impute(qfeatures,
           i = assayName,
           method = "knn",
           k = 3, rowmax = 1, colmax = 1,
           name = "imputed_proteins"
    )
}

leduc <- scpdata::leduc2022_pSCoPE()

PSM1000 <- generateTMTPSM(leduc, 1000)
PEP1000 <- generateTMTPeptides(PSM1000)
PRO1000 <- generateTMTProteins(PEP1000)

aggregateProf <- profvis(expr = {
    aggregateFeaturesOverAssays(
        PSM1000,
        i = seq_along(PSM1000),
        "modseq",
        name = paste0("peptide_", seq_along(PSM1000)),
        fun = MsCoreUtils::robustSummary
    )
})

aggPSM1000 <- aggregateFeaturesOverAssays(
    PSM1000,
    i = seq_along(PSM1000),
    "modseq",
    name = paste0("peptide_", seq_along(PSM1000)),
    fun = MsCoreUtils::robustSummary
)

joinProf <- profvis(expr = {
    joinAssays(aggPSM1000,
               i = paste0("peptide_", 1:(length(aggPSM1000) / 2)),
               name = "peptides"
    )
})

imputeProf <- profvis(expr = {
    impute(PRO1000,
           i = "proteins",
           method = "knn",
           k = 3, rowmax = 1, colmax = 1,
           name = "imputed_proteins"
    )
})


saveRDS(aggregateProf, file = "dataOutput/profRes/aggregateProf.rds")

saveRDS(joinProf, file = "dataOutput/profRes/joinProf.rds")

saveRDS(imputeProf, file = "dataOutput/profRes/imputeProf.rds")
