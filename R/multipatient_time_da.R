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
colData(pepse) <- colData(pepse) %>%
  as.data.frame() %>%
  separate_wider_delim(cols = patient_id, delim = "_", names = c("patient_id", "dayp")) %>%
  as("DataFrame")


# Processing

pepse <- zeroIsNA(pepse)

pepse <- filterNA(pepse, pNA = 0.98)

pepse <- sweep(pepse,
      MARGIN = 2,
      FUN = "/",
      STATS = colMedians(assay(pepse), na.rm = TRUE))

pepseMod <- scpModelWorkflow(pepse, formula = ~ 1 + patient + day + raw.file + lc_batch)

dayRes <- scpDifferentialAnalysis(
  pepseMod,
  contrasts = list(c("day", "d0", "d"))
)[[1]] %>%
  filter(padj <= 0.05,
         Estimate >= 1) %>%
  rownames()


