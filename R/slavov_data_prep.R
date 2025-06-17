library(tidyverse)
library(QFeatures)

pep <- read.csv("/mnt/disk3/slavov_data/scope_pep-001.csv")
coldata <- read.csv("/mnt/disk3/slavov_data/scprot_metadata.csv")
coldata$quantCols <- coldata$sc_id
coldata$runCols <- coldata$raw.file
rownames(coldata) <- coldata$sc_id
subpep <- pep[, c("modseq", "sc_id", "pep_quant")]
subpep <- subpep[subpep$sc_id %in% rownames(coldata),]
wide <- pivot_wider(subpep, names_from = "sc_id", values_from = "pep_quant")

coldata
readQFeatures(pep, coldata, quantCols = "pep_quant", fnames = "unique_id")

se <- readSummarizedExperiment(wide, quantCols = 2:13145, fnames = "modseq")

