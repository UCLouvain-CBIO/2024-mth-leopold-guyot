library(msqrob2)
library(scp)
library(iCOBRA)
source("R/multipatient_simulation.R")
library(patchwork)

compute_performance <- function(modRes, rowdata) {
  df <- data.frame()
  for (name in names(modRes)) {
    curr <- data.frame(
      sid = rownames(modRes[[name]]),
      pval = modRes[[name]]$pval,
      adjPval = modRes[[name]]$adjPval,
      is_da = rowdata[rownames(modRes[[name]]), "shifted"],
      method = rep(name, nrow(modRes[[name]]))
    )
    df <- rbind(df, curr)
  }
  rownames(df) <- NULL
  truth_df <- df %>%
    select(sid, is_da) %>%
    distinct() %>%
    column_to_rownames("sid")

  adjPval_df <- df %>%
    select(sid, method, adjPval) %>%
    pivot_wider(names_from = method, values_from = adjPval) %>%
    column_to_rownames("sid")

  cobdata <- COBRAData(
    padj = as.data.frame(adjPval_df),
    truth = as.data.frame(truth_df)
  )

  perf <- calculate_performance(cobdata,
                                binary_truth = "is_da",
                                aspects = "fdrtpr",
                                maxsplit = Inf,
                                thrs = seq(from = 0.0001,  to = 0.2, by = 0.0001)
  ) %>%
    fdrtpr() %>%
    mutate(thr = as.numeric(sub("thr", "", thr)))
  return(perf)
}

aggregate_se <- function(se, group_by_cols, fun = mean) {
  if (!all(group_by_cols %in% colnames(colData(se)))) {
    stop("Some specified group_by_cols do not exist in colData.")
  }

  colData(se)$Group <- apply(colData(se)[, group_by_cols, drop = FALSE], 1, paste, collapse = "_")
  group_factor <- factor(colData(se)$Group)
  assay_matrix <- assay(se)

  aggregated_assay <- do.call(cbind, lapply(split(seq_along(group_factor), group_factor),
                                            function(idx) apply(assay_matrix[, idx, drop = FALSE], 1, fun, na.rm = TRUE)))

  new_colData <- colData(se) %>%
    as.data.frame() %>%
    group_by(Group) %>%
    summarise(across(all_of(group_by_cols), dplyr::first), .groups = "drop") %>%
    as.data.frame()


  rownames(new_colData) <- new_colData$Group

  new_colData <- new_colData[colnames(aggregated_assay),]
  new_se <- SummarizedExperiment(assays = list(counts = aggregated_assay), colData = new_colData)

  return(new_se)
}

## Explore model coefficients

sigma <- lapply(mod, function(x) x$sigma)
sigmadf <- data.frame(name = names(sigma), sigma = as.numeric(sigma))
coef <- lapply(mod, function(x) x$coefficients)
coefdf <- as.data.frame(do.call(rbind, coef))
ggplot(sigmadf, aes(x = sigma)) + geom_density()

coefdf %>%
  filter(grepl("cell_type", rownames(coefdf))) %>%
  ggplot(aes(x = Estimate)) + geom_density()

coefdf %>%
  filter(grepl("patient", rownames(coefdf))) %>%
  ggplot(aes(x = Estimate)) + geom_density()


## Simulate data

simData <- simulate_peptide_data(
  mod,
  n_cells_per_comb = 50,
  cell_types = c("CT1", "lowerresCD8T", "lowerresmonocyte", "lowerresNK"),
  n_synthetic_patients = 16
)

rownames(simData) <- rowData(simData)$peptide

modscp <- scpModelWorkflow(simData, formula = ~ patient_id + cell_type)
pca <- scpComponentAnalysis(modscp, method = "APCA", residuals = FALSE, unmodelled = FALSE)
colnames(colData(simData))[[1]] <- "cell"
bySamplePCs <- scpAnnotateResults(
  pca$bySample, colData(simData), by = "cell"
)

scpComponentPlot(
  bySamplePCs,
  pointParams = list(
    aes(colour = patient_id, shape = cell_type),
    alpha = 0.6
  )
) |>
  wrap_plots(guides = "collect")

var <- scpVarianceAnalysis(modscp)
scpVariancePlot(var)

daSimData <- add_patient_group_shift_SE(simData,
                                        group_a = unique(colData(simData)$patient_id)[1:8],
                                        group_b = unique(colData(simData)$patient_id)[9:16],
                                        shift = 3,
                                        sd = 0.2,
                                        ratio = 0.1,
                                        seed = 123)

daSimDataAggregated <- aggregate_se(daSimData,
                                    group_by_cols = c("cell_type", "patient_id", "condition"),
                                    fun = mean)

scpModel <- scpModelWorkflow(daSimData, formula = ~ 1 + cell_type + condition, verbose = TRUE)
scpRes <- scpDifferentialAnalysis(
  scpModel,
  contrasts = list(c("condition", "A", "B"))
)[[1]]

colnames(scpRes) <- c("feature", "logFC", "se", "df", "t", "pval", "adjPval")
rownames(scpRes) <- scpRes$feature

L <- makeContrast("conditionB=0", parameterNames = c("conditionB"))
msqModel <- suppressWarnings(msqrob(daSimData, formula = ~ cell_type + condition + (1 | patient_id)))
msqRes <- rowData(hypothesisTest(object = msqModel, contrast = L))$conditionB

pbModel <- suppressWarnings(msqrob(daSimDataAggregated, ~ 1 + condition + cell_type))
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
