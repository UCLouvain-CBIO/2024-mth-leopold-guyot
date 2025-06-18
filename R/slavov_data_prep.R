library(tidyverse)
library(QFeatures)
library(parallel)
library(lme4)

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
subpep <- pep[, c("modseq", "sc_id", "pep_quant")]
subpep <- subpep[subpep$sc_id %in% rownames(coldataJoined),]
widepep <- pivot_wider(subpep, names_from = "sc_id", values_from = "pep_quant")

pepse <- readSummarizedExperiment(widepep, quantCols = 2:ncol(widepep), fnames = "modseq")
coldataJoined <- coldataJoined[coldataJoined$sc_id %in% colnames(pepse), ]
colData(pepse) <- as(coldataJoined, "DataFrame")
pepse <- logTransform(pepse)

# proteins

subprot <- prot[, c("protein", "sc_id", "prot_quant")]
subprot <- subprot[subprot$sc_id %in% rownames(coldataJoined),]
wideprot <- pivot_wider(subprot, names_from = "sc_id", values_from = "prot_quant")

protse <- readSummarizedExperiment(wideprot, quantCols = 2:ncol(wideprot), fnames = "protein")
coldataJoined <- coldataJoined[coldataJoined$sc_id %in% colnames(protse), ]
colData(protse) <- as(coldataJoined, "DataFrame")
protse <- logTransform(protse)

# model peptides

longpep <- longForm(pepse, colvars = c("sc_id", "cell_type_lowerres", "patient_id", "raw.file", "day"))
longpepD0 <- longpep[longpep$day == "d0", ]
longpepD0 <- longpepD0[longpepD0$cell_type_lowerres %in% c("CD4T", "CD8T", "monocyte", "NK"),]

splitpep <- split(longpepD0, longpepD0$rowname)

modelList <- mclapply(splitpep, FUN = function(pep) {
  tryCatch({
    mod <- lmer(formula = value ~ cell_type_lowerres + patient_id + (1 | raw.file), as.data.frame(pep))
    list("coefficients" = summary(mod)$coefficients, "sigma" = summary(mod)$sigma)  
  },
           error = function(e) NA)
  },
  mc.cores = 12L
)

saveRDS(modelList, "dataOutput/slavovModels/pepMod.RDS")
