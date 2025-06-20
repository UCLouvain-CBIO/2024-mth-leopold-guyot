library(tidyverse)
library(QFeatures)
library(parallel)
library(lme4)
library(pcaMethods)
library(msqrob2)
source(file = "R/utils.R")

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


# Implement known changes

shift_ratio <- 0.1
shift_coef <- 0.8
shift_sd <- 0.4

is_d2 <- colData(pepse)$day == "d2"

mat <- assay(pepse)

shift_flag <- runif(nrow(mat)) < shift_ratio

peptide_shifts <- numeric(nrow(mat))
peptide_shifts[shift_flag] <- rnorm(
  sum(shift_flag),
  mean = shift_coef,
  sd = shift_sd
)

mat[, is_d2] <- sweep(
  mat[, is_d2, drop = FALSE],
  MARGIN = 1,
  STATS = peptide_shifts,
  FUN = "+"
)

rowData(pepse)$shift_value <- peptide_shifts

assay(pepse) <- mat

pepsePB <- aggregate_se(pepse,
             group_by_cols = c("cell_type_lowerres", "patient", "day"),
             fun = mean)

scpModel <- scpModelWorkflow(pepse, formula = ~ 1 + cell_type_lowerres + patient + day, verbose = TRUE)
scpRes <- scpDifferentialAnalysis(
  scpModel,
  contrasts = list(c("day", "d0", "d2"))
)[[1]]

colnames(scpRes) <- c("feature", "logFC", "se", "df", "t", "pval", "adjPval")
rownames(scpRes) <- scpRes$feature

L <- makeContrast("dayd2=0", parameterNames = c("dayd2"))
msqModel <- suppressWarnings(msqrob(pepse, formula = ~ cell_type + condition + (1 | patient_id)))
msqRes <- rowData(hypothesisTest(object = msqModel, contrast = L))$conditionB

pbModel <- suppressWarnings(msqrob(pepsePB, ~ 1 + condition + cell_type))
pbRes <- rowData(hypothesisTest(object = pbModel, contrast = L))$conditionB

fdrtpr <- compute_performance(list("scp" = scpRes, "msqrob2" = msqRes, "pseudobulk" = pbRes), rowdata = rowData(daSimData))

plot <- ggplot(fdrtpr, aes(x = FDR, y = TPR, color = method)) +
  geom_vline(
    xintercept = c(0.01, 0.05, 0.1),
    linetype = "dashed", color = "grey50", linewidth = 0.3
  ) +
  geom_point(size = 0.5, alpha = 0.8) +
  geom_line(size = 0.7) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0.01, 0.05, 0.2, 0.4, 0.8, 1),
    labels = scales::label_number()
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  )

ggsave("Figs/fdrtpr.pdf")
