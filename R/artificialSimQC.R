library(scp)
library(patchwork)
library(pcaMethods)
source("R/artificialSimulationPB.R")

base <- scpdata::brunner2022()
sce <- getWithColData(base, "proteins")
sce <- sce[, colData(sce)$CellCycleStage == "UB"]
sce <- logTransform(sce) # log transform to obtain normal intensities distribution
sce <- sce[rowSums(!is.na(assay(sce))) > 1, ]
sce

plot_data <- do.call(rbind, lapply(top_peptides, function(peptide) {
    values <- assay(sce)[peptide, ]
    values <- values[!is.na(values)]
    data.frame(
        Intensity = values,
        Peptide = peptide,
        Mean = mean(values),
        SD = sd(values)
    )
}))

normal_lines <- plot_data %>%
    group_by(Peptide) %>%
    summarise(
        Mean = unique(Mean),
        SD = unique(SD)
    ) %>%
    rowwise() %>%
    mutate(
        x = list(seq(Mean - 4*SD, Mean + 4*SD, length.out = 200)),
        y = list(dnorm(x, mean = Mean, sd = SD))
    ) %>%
    tidyr::unnest(cols = c(x, y))

ggplot(plot_data, aes(x = Intensity)) +
    geom_density(fill = "lightblue", alpha = 0.6) +
    geom_line(data = normal_lines, aes(x = x, y = y), color = "red", linetype = "dashed") +
    facet_wrap(~Peptide, ncol = 1, scales = "free_y") +
    labs(title = "Empirical vs. Fitted Normal Distribution (Top 3 Peptides)",
         x = "Log-Intensity",
         y = "Density") +
    theme_minimal()

ggsave("Figs/report/brunnerDensity.pdf", width = 5, height = 10)


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

pca <- pca(t(assay(simSCE)), method = "nipals")
df <- merge(scores(pca), colData(simSCE))

ggplot(df, aes(PC1, PC2, color = CellType)) +
    geom_point(alpha = 0.3) +
    xlab(paste("PC1", pca@R2[1] * 100, "% of the variance")) +
    ylab(paste("PC2", pca@R2[2] * 100, "% of the variance"))

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

ggsave("Figs/report/simApcaPatient.pdf", width = 10, height = 6)

plotCell <- scpComponentPlot(
    bySamplePCs[2],
    pointParams = list(
        aes(colour = CellType),
        alpha = 0.6
    )
)

ggsave("Figs/report/simApcaCell.pdf", width = 10, height = 6)


