source("R/artificialSimulationPB.R")

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

simSCE <- simulateCellPatientData(protMetrics, rowdata = rowData(sce),
                                  nPatient = 16, nPopulation = 5,  nCellPatPop = 25,
                                  patientEffect = 1, patientShift = 0.2, patientSD = 0.05,
                                  populationEffect = 0.33, populationShift = 0.4, populationSD = 0.1,
                                  seed = 123)


scpMod <- scpModelWorkflow(simSCE, formula = ~ 1 + Patient + CellType)
