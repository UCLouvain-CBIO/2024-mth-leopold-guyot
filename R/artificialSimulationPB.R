library(scpdata)
library(scp)
library(SingleCellExperiment)
library(msqrob2)
library(iCOBRA)
library(UpSetR)
library(tidyverse)
library(MsCoreUtils)


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
  
  assay(sce) <- assayMatrix
  
  return(sce)
}

customRobustSummary <- function(x, ...) {
  apply(x, 1, function(row) MASS::rlm(row ~ 1, ...)$coefficients[1])
}

aggregate_se <- function(se, group_by_cols, fun = mean) {
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
                                              apply(assay_matrix[, idx, drop = FALSE], 1, fun, na.rm = TRUE)
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

t <- simulateCellPatientData(protMetrics, rowdata = rowData(sce),
                             nPatient = 4, nPopulation = 4,  nCellPatPop = 10,
                             patientEffect = 1, patientShift = 0.3, patientSD = 0.1,
                             populationEffect = 0.3, populationShift = 0.3, populationSD = 0.1,
                             seed = 123)

treat <- addTreatmentEffect(t, expMetrics = protMetrics, treatmentEffect = 0.3, treatmentShift = 0.3, treatmentSD = 0.1
                   )
