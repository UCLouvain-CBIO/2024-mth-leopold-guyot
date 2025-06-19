library(tidyverse)
library(QFeatures)
library(parallel)
library(lme4)
library(scp)
library(pcaMethods)

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

pca <- pca(t(assay(pepse)), "nipals")

df <- merge(scores(pca), colData(pepse), by = 0)
ggplot(df, aes(PC1, PC2, shape=cell_type_lowerres, color=lc_batch)) +
  geom_point() +
  xlab(paste("PC1", pca@R2[1] * 100, "% of the variance")) +
  ylab(paste("PC2", pca@R2[2] * 100, "% of the variance")) + 
  theme(legend.position = "none")

saveRDS(df, "dataOutput/slavovModels/pca.rds")  
