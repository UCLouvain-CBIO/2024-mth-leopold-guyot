library(tidyverse)
library(QFeatures)
library(parallel)
library(lme4)
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
colData(pepse) <- colData(pepse) %>%
  as.data.frame() %>%
  separate_wider_delim(cols = patient_id, delim = "_", names = c("patient", "dayp")) %>%
  as("DataFrame")
colnames(pepse) <- colData(pepse)$sc_id

# Processing
pepse <- pepse[, colData(pepse)$day %in% c("d0", "d2")]
pepse <- zeroIsNA(pepse)

pepse <- filterNA(pepse, pNA = 0.98)

pepse <- sweep(pepse,
      MARGIN = 2,
      FUN = "/",
      STATS = colMedians(assay(pepse), na.rm = TRUE))
pca <- pcaMethods::nipalsPca(t(assay(pepse)), maxSteps = 5000)
df <- merge(scores(pca), colData(pepse), by = 0)

ggplot(df, aes(x = V1, y = V2, color = patient)) + geom_point() + xlim(c(-500, 500)) + ylim(c(-500, 300))

pepseMod <- scpModelWorkflow(pepse, formula = ~ 1 + patient + day + lc_batch + cell_type_lowerres)

dayRes <- scpDifferentialAnalysis(
  pepseMod,
  contrasts = list(c("day", "d0", "d2"))
)[[1]] %>%
  as.data.frame() %>% 
  filter(padj <= 0.05)

pepse <- pepse[-dayRes$feature,]

