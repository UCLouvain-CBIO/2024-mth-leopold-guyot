library(tidyverse)
library(QFeatures)
library(parallel)
library(lme4)
library(scp)
library(pcaMethods)
library(patchwork)

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

pca <- pca(t(assay(pepse)), "nipals", )

df <- merge(scores(pca), colData(pepse), by = 0)
combinedPlot <- list()
for (var in c("label", "plate", "patient_id", "lc_batch", "day", "cell_type_lowerres", "raw.file")) {
  combinedPlot[[var]] <- ggplot(df, aes(PC1, PC2, color=.data[[var]])) +
    geom_point() +
    xlab(paste("PC1", pca@R2[1] * 100, "% of the variance")) +
    ylab(paste("PC2", pca@R2[2] * 100, "% of the variance")) + 
    theme(legend.position = "none") + ggtitle(var)
}

ggsave(wrap_plots(combinedPlot), "Figs/nipalsMultiPatient.pdf", width = 15, height = 10)

scpMod <- scpModelWorkflow(pepse, formula = ~ 1 + patient_id + cell_type_lowerres)

saveRDS(df, "dataOutput/slavovModels/pca.rds")  
