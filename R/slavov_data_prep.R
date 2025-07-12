library(tidyverse)
library(QFeatures)
library(parallel)
library(lme4)
library(limma)

pep <- read.csv("/mnt/disk3/slavov_data/scope_pep-001.csv")
prot <- read.csv("/mnt/disk3/slavov_data/scope_prot_max.csv")
coldata <- read.csv("/mnt/disk3/slavov_data/scprot_metadata.csv")

#coldata
coldata$quantCols <- coldata$sc_id
coldata$runCols <- coldata$raw.file
rownames(coldata) <- coldata$sc_id
subpepColdata <- pep[, c("sc_id", "raw.file")] %>% 
  group_by(sc_id) %>% 
  summarize("raw.file" = unique(raw.file)) %>% 
  filter(sc_id %in% coldata$sc_id)
coldataJoined <- left_join(coldata, subpepColdata, by = "sc_id")
rownames(coldataJoined) <- rownames(coldata)

# peptides
subpep <- pep[, c("modseq", "sc_id", "quantitation")]
subpep <- subpep[subpep$sc_id %in% rownames(coldataJoined),]
widepep <- pivot_wider(subpep, names_from = "sc_id", values_from = "quantitation")

pepse <- readSummarizedExperiment(widepep, quantCols = 2:ncol(widepep), fnames = "modseq")
coldataJoined <- coldataJoined[coldataJoined$sc_id %in% colnames(pepse), ]
colData(pepse) <- as(coldataJoined, "DataFrame")
pepse <- logTransform(pepse)
pepse <- normalize(pepse, method = "center.median")
pepse <- pepse[, colData(pepse)$day =="d0"]
pepse <- pepse[, colData(pepse)$cell_type_lowerres %in% c("CD4T", "CD8T", "monocyte", "NK")]


# proteins

subprot <- prot[, c("protein", "sc_id", "prot_quant")]
subprot <- subprot[subprot$sc_id %in% rownames(coldataJoined),]
wideprot <- pivot_wider(subprot, names_from = "sc_id", values_from = "prot_quant")

protse <- readSummarizedExperiment(wideprot, quantCols = 2:ncol(wideprot), fnames = "protein")
coldataJoined <- coldataJoined[coldataJoined$sc_id %in% colnames(protse), ]
colData(protse) <- as(coldataJoined, "DataFrame")
protse <- logTransform(protse)

# Explore with scplainer
# 
# modscp <- scpModelWorkflow(pepse, formula = ~ patient_id + cell_type_lowerres)
# pca <- scpComponentAnalysis(modscp, method = "APCA", residuals = FALSE, unmodelled = FALSE)
# colData(pepse)$cell <- colData(pepse)$sc_id
# bySamplePCs <- scpAnnotateResults(
#   pca$bySample, colData(pepse), by = "cell"
# )
# 
# scpComponentPlot(
#   bySamplePCs,
#   pointParams = list(
#     aes(colour = patient_id, shape = cell_type_lowerres),
#     alpha = 0.6
#   )
# ) |>
#   wrap_plots(guides = "collect")
# 
# var <- scpVarianceAnalysis(modscp)
# scpVariancePlot(var)
# model peptides

longpep <- longForm(pepse, colvars = c("sc_id", "cell_type_lowerres", "patient_id", "raw.file", "day", "plate", "lc_batch", "leiden"))

splitpep <- split(longpep, longpep$rowname)
modelList <- mclapply(splitpep, FUN = function(pep) {
  tryCatch({
    mod <- lmer(formula = value ~ cell_type_lowerres + patient_id + (1|leiden), as.data.frame(pep))
    list("coefficients" = summary(mod)$coefficients, "sigma" = summary(mod)$sigma)  
  },
           error = function(e) NA)
  },
  mc.cores = 12L
)

saveRDS(modelList, "dataOutput/slavovModels/pepMod.RDS")
