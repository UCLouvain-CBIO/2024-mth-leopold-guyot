library(scpdata)
library(scp)
library(SingleCellExperiment)
library(msqrob2)
library(iCOBRA)
library(UpSetR)
library(tidyverse)
library(MsCoreUtils)
library(BiocParallel)
library(ggplot2)

register(BPPARAM = MulticoreParam(workers = 4))

compute_performance <- function(modRes, rowdata) {
  df <- data.frame()
  for (name in names(modRes)) {
    curr <- data.frame(
      sid = rownames(modRes[[name]]),
      pval = modRes[[name]]$pval,
      adjPval = modRes[[name]]$adjPval,
      is_da = rowdata[rownames(modRes[[name]]), "TreatmentShiftedProtein"],
      method = rep(name, nrow(modRes[[name]]))
    )
    df <- rbind(df, curr)
  }
  rownames(df) <- NULL
  truth_df <- df %>%
    select(sid, is_da) %>%
    distinct() %>%
    column_to_rownames("sid")

  adjPval_df <- df %>%
    select(sid, method, adjPval) %>%
    pivot_wider(names_from = method, values_from = adjPval) %>%
    column_to_rownames("sid")

  cobdata <- COBRAData(
    padj = as.data.frame(adjPval_df),
    truth = as.data.frame(truth_df)
  )

  perf <- calculate_performance(cobdata,
                                binary_truth = "is_da",
                                aspects = "fdrtpr",
                                maxsplit = Inf,
                                thrs = seq(from = 0.001,  to = 1, by = 0.001)
  ) %>%
    fdrtpr() %>%
    mutate(thr = as.numeric(sub("thr", "", thr)))
  return(perf)
}


simulateCellPatientData <- function(expMetrics, rowdata,
                                    nPatient, nPopulation,  nCellPatPop,
                                    patientEffect, patientShift, patientSD,
                                    populationEffect, populationShift, populationSD,
                                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  mu <- expMetrics$meanIntensity
  proteinNames <- expMetrics$protNames
  residualSD <- expMetrics$sdIntensity
  missingProb <- expMetrics$pctMissing/100
  nProt <- length(mu)

  nCellTot <- nPatient * nPopulation * nCellPatPop

  # Initialize matrix
  simMatrix <- matrix(NA, nrow = nProt, ncol = nCellTot)
  rownames(simMatrix) <- proteinNames

  nPatientNames <- paste0("Patient", 1:nPatient)
  nPopulationNames <- paste0("CellType", 1:nPopulation)
  nCellPatPop <- 1:nCellPatPop

  dfCellId <- expand.grid(
    Patient = nPatientNames,
    CellType = nPopulationNames,
    Cell = nCellPatPop
  )
  cellIds <- with(dfCellId, paste(Patient, CellType, Cell, sep = "_"))
  colnames(simMatrix) <- cellIds

  # Number of proteins to shift
  nPatientShiftProt <- floor(patientEffect * nProt)
  nPopulationShiftProt <- floor(populationEffect * nProt)

  # Select proteins affected by patient and population effect
  patientShiftProteins <- sample(proteinNames, nPatientShiftProt)
  populationShiftProteins <- sample(proteinNames, nPopulationShiftProt)

  for (i in 1:nProt) {
    simMatrix[i, ] <- rnorm(nCellTot, mean = mu[i], sd = residualSD[i])
  }

  # Apply patient effect shifts
  for (prot in patientShiftProteins) {
    protIdx <- which(proteinNames == prot)

    for (pat in nPatientNames) {
      shiftVal <- rnorm(1,
                        mean = patientShift * mu[protIdx],
                        sd = patientSD * mu[protIdx])

      cols <- grep(paste0("^", pat, "_"), colnames(simMatrix))
      simMatrix[protIdx, cols] <- simMatrix[protIdx, cols] + shiftVal
    }
  }

  # Apply population (cell type) effect shifts
  for (prot in populationShiftProteins) {
    protIdx <- which(proteinNames == prot)

    for (pop in nPopulationNames) {
      shiftVal <- rnorm(1,
                        mean = populationShift * mu[protIdx],
                        sd = populationSD * mu[protIdx])

      cols <- grep(paste0("_", pop, "_"), colnames(simMatrix))
      simMatrix[protIdx, cols] <- simMatrix[protIdx, cols] + shiftVal
    }
  }

  # Add missingness at random
  for (p in seq_len(nProt)) {
    missing <- runif(nCellTot) < missingProb[p]
    simMatrix[p, missing] <- NA
  }
  simSCE <- SingleCellExperiment(assays = SimpleList(simMatrix),
                                 rowData = rowdata,
                                 colData = dfCellId)
  simSCE
}

addTreatmentEffect <- function(sce, expMetrics,
                               treatmentEffect,
                               treatmentShift,
                               treatmentSD,
                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  mu <- expMetrics$meanIntensity
  proteinNames <- expMetrics$protNames
  nProt <- length(mu)

  nTreatmentShiftProt <- floor(treatmentEffect * nProt)

  treatmentShiftProteins <- sample(proteinNames, nTreatmentShiftProt)

  # Pick half the patients (randomly)
  patientLevels <- unique(colData(sce)$Patient)
  nPatients <- length(patientLevels)
  nTreated <- floor(nPatients / 2)

  treatedPatients <- sample(patientLevels, nTreated)

  assayMatrix <- assay(sce)

  # Apply shifts
  for (prot in treatmentShiftProteins) {
    protIdx <- which(rownames(assayMatrix) == prot)

    for (pat in treatedPatients) {
      shiftVal <- rnorm(1,
                        mean = treatmentShift * mu[protIdx],
                        sd = treatmentSD * mu[protIdx])

      cols <- which(colData(sce)$Patient == pat)

      assayMatrix[protIdx, cols] <- assayMatrix[protIdx, cols] + shiftVal
    }
  }

  # Save treated/untreated status in colData
  colData(sce)$Treatment <- ifelse(colData(sce)$Patient %in% treatedPatients,
                                   "Treated", "Control")
  # Add logical column for shifted proteins
  rowData(sce)$TreatmentShiftedProtein <- rownames(sce) %in% treatmentShiftProteins

  assay(sce) <- assayMatrix

  return(sce)
}

customRobustSummary <- function(x, ...) {
  # For each feature (row), estimate robust mean
  res <- apply(x, 1, function(row) {
    if (all(is.na(row))) return(NA_real_)
    fit <- MASS::rlm(row ~ 1, ...)
    coef(fit)[1]
  })

  names(res) <- rownames(x)
  res
}


aggregate_se <- function(se, group_by_cols, fun = mean, robustSummary = FALSE) {
  if (!all(group_by_cols %in% colnames(colData(se)))) {
    stop("Some specified group_by_cols do not exist in colData.")
  }
  colData(se)$Group <- apply(colData(se)[, group_by_cols, drop = FALSE], 1, paste, collapse = "_")
  group_factor <- factor(colData(se)$Group)
  assay_matrix <- assay(se)

  # Aggregate assay matrix by group_factor levels order
  grouped_indices <- split(seq_along(group_factor), group_factor)
  levels_group <- levels(group_factor)  # factor levels in order
  aggregated_assay <- do.call(cbind, lapply(levels_group,
                                            function(g) {
                                              idx <- grouped_indices[[g]]
                                              if (robustSummary) {
                                                customRobustSummary(assay_matrix[, idx, drop = FALSE])
                                              } else {
                                                apply(assay_matrix[, idx, drop = FALSE], 1, fun, na.rm = TRUE)
                                              }
                                            }))
  new_colData <- colData(se) %>%
    as.data.frame() %>%
    group_by(Group) %>%
    summarise(across(all_of(group_by_cols), dplyr::first), .groups = "drop") %>%
    as.data.frame()
  # Reorder new_colData rows to match factor levels
  new_colData <- new_colData[match(levels_group, new_colData$Group), ]
  rownames(new_colData) <- new_colData$Group
  new_se <- SummarizedExperiment(assays = list(aggregated_assay),
                                 colData = new_colData,
                                 rowData = rowData(se))
  return(new_se)
}


benchmarkMethods <- function(expMetrics, rowdata,
                             nPatient, nPopulation,  nCellPatPop,
                             patientEffect, patientShift, patientSD,
                             populationEffect, populationShift, populationSD,
                             treatmentEffect, treatmentShift, treatmentSD,
                             seed = NULL) {

  simSCE <- simulateCellPatientData(expMetrics = protMetrics, rowdata = rowdata,
                          nPatient = nPatient, nPopulation = nPopulation,  nCellPatPop = nCellPatPop,
                          patientEffect = patientEffect, patientShift = patientShift, patientSD = patientSD,
                          populationEffect = populationEffect, populationShift = populationShift, populationSD = populationSD,
                          seed = seed)

  simSCE <- addTreatmentEffect(simSCE, expMetrics = protMetrics,
                               treatmentEffect = treatmentEffect, treatmentShift = treatmentShift, treatmentSD = treatmentSD, seed = seed)

  aggMean <- aggregate_se(simSCE, group_by_cols = c("CellType", "Patient", "Treatment"), fun = mean)
  aggMedian <- aggregate_se(simSCE, group_by_cols = c("CellType", "Patient", "Treatment"), fun = median)
  aggSum <- aggregate_se(simSCE, group_by_cols = c("CellType", "Patient", "Treatment"), fun = sum)
  aggRobust <- aggregate_se(simSCE, group_by_cols = c("CellType", "Patient", "Treatment"), fun = NULL, robustSummary = TRUE)

  scpModel <- scpModelWorkflow(simSCE, formula = ~ 1 + CellType + Treatment, verbose = FALSE)
  scpRes <- scpDifferentialAnalysis(
    scpModel,
    contrasts = list(c("Treatment", "Control", "Treated"))
  )[[1]]

  colnames(scpRes) <- c("feature", "logFC", "se", "df", "t", "pval", "adjPval")
  rownames(scpRes) <- scpRes$feature

  L <- makeContrast("TreatmentTreated=0", parameterNames = c("TreatmentTreated"))
  msqModel <- suppressMessages(suppressWarnings(msqrob(simSCE, formula = ~ CellType + Treatment + (1 | Patient))))
  msqRes <- rowData(hypothesisTest(object = msqModel, contrast = L))$TreatmentTreated

  aggMeanMod <- suppressWarnings(msqrob(aggMean, ~ 1 + Treatment + CellType))
  aggMeanRes <- rowData(hypothesisTest(object = aggMeanMod, contrast = L))$TreatmentTreated

  aggMedianMod <- suppressWarnings(msqrob(aggMedian, ~ 1 + Treatment + CellType))
  aggMedianRes <- rowData(hypothesisTest(object = aggMedianMod, contrast = L))$TreatmentTreated

  aggSumMod <- suppressWarnings(msqrob(aggSum, ~ 1 + Treatment + CellType))
  aggSumRes <- rowData(hypothesisTest(object = aggSumMod, contrast = L))$TreatmentTreated

  aggRobustMod <- suppressWarnings(msqrob(aggRobust, ~ 1 + Treatment + CellType))
  aggRobustRes <- rowData(hypothesisTest(object = aggRobustMod, contrast = L))$TreatmentTreated

  fdrtpr <- compute_performance(list("scp" = scpRes,
                                     "msqrob2" = msqRes,
                                     "pseudobulkMean" = aggMeanRes,
                                     "pseudobulkMedian" = aggMedianRes,
                                     "pseudobulkSum" = aggSumRes,
                                     "pseudobulkRobustSummary" = aggRobustRes
                                     ), rowdata = rowData(simSCE))

}

base <- scpdata::brunner2022()
sce <- getWithColData(base, "proteins")
sce <- sce[, colData(sce)$CellCycleStage == "UB"]
sce <- logTransform(sce) # log transform to obtain normal intensities distribution
sce <- sce[rowSums(!is.na(assay(sce))) > 1, ]
sce

protAssay <- as.data.frame(assay(sce))

# Compute metric
meanIntensity <- rowMeans(protAssay, na.rm = TRUE)
sdIntensity <- apply(protAssay, 1, sd, na.rm = TRUE)
pctMissing <- rowMeans(is.na(protAssay)) * 100
protMetrics <- data.frame(
  protNames = rownames(protAssay),
  meanIntensity = meanIntensity,
  sdIntensity = sdIntensity,
  pctMissing = pctMissing,
  stringsAsFactors = FALSE
)

benchRes <- list()

for (nCellPatPop in c(10, 25, 50, 100)) {
  for (treatmentShift in c(0.05, 0.1, 0.15, 0.2)) {
    cat("Starting simulation:", "nCell = ", nCellPatPop, ", shift = ", treatmentShift, "\n")
    benchRes[[paste0("nCell", nCellPatPop, "_", "shift", treatmentShift)]] <-
      benchmarkMethods(protMetrics, rowdata = rowData(sce),
                     nPatient = 16, nPopulation = 5,  nCellPatPop = nCellPatPop,
                     patientEffect = 1, patientShift = 0.2, patientSD = 0.05,
                     populationEffect = 0.33, populationShift = 0.4, populationSD = 0.1,
                     treatmentEffect = 0.33, treatmentShift = treatmentShift, treatmentSD = treatmentShift/4,
                     seed = 123)
  }
}

saveRDS(benchRes, file = "dataOutput/artificialSimPB/tprfdrRes.rds")
