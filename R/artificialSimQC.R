library(scp)
library(patchwork)

source("R/artificialSimulationPB.R")

base <- scpdata::brunner2022()
sce <- getWithColData(base, "proteins")
sce <- sce[, colData(sce)$CellCycleStage == "UB"]
sce <- logTransform(sce) # log transform to obtain normal intensities distribution
sce <- sce[rowSums(!is.na(assay(sce))) > 1, ]
sce

all_values <- as.vector(assay(sce))
df <- data.frame(Expression = all_values)
ggplot(df, aes(x = Expression)) +
    geom_density(fill = "lightblue", alpha = 0.6) +
    labs(title = "Density Plot of All Log-Intensity Values",
         x = "Expression",
         y = "Density") +
    theme_minimal()

ggsave("Figs/report/brunnerDensity.pdf")


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
                                  patientEffect = 1, patientShift = 0, patientSD = 0.1,
                                  populationEffect = 0.33, populationShift = 0, populationSD = 0.2,
                                  seed = 123)


scpMod <- scpModelWorkflow(simSCE, formula = ~ 1 + Patient + CellType)
colData(scpMod)$cell <- rownames(colData(scpMod))
pcs <- scpComponentAnalysis(scpMod, method = "APCA", residuals = FALSE, unmodelled = FALSE)

bySamplePCs <- scpAnnotateResults(
    pcs$bySample, colData(scpMod), by = "cell"
)

plotPatient <- scpComponentPlot(
    bySamplePCs[1],
    pointParams = list(
        aes(colour = Patient),
        alpha = 0.6
    )
)

ggsave("Figs/report/simApcaPatient.pdf")

plotCell <- scpComponentPlot(
    bySamplePCs[2],
    pointParams = list(
        aes(colour = CellType),
        alpha = 0.6
    )
)

ggsave("Figs/report/simApcaCell.pdf")


