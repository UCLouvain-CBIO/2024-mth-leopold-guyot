library(msqrob2)
library(scp)
source("R/multipatient_simulation.R")

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

simData <- simulate_peptide_data(
  mod,
  n_cells_per_comb = 150,
  cell_types = c("CT1", "lowerresCD8T", "lowerresmonocyte", "lowerresNK"),
  n_synthetic_patients = 16
)


daSimData <- add_patient_group_shift_SE(simData,
                                        group_a = unique(colData(simData)$patient_id)[1:8],
                                        group_b = unique(colData(simData)$patient_id)[9:16],
                                        shift = 1,
                                        sd = 0.5,
                                        ratio = 0.1,
                                        seed = 123)

daSimDataAggregated <- aggregate_se(daSimData,
                                    group_by_cols = c("cell_type", "patient_id"),
                                    fun = mean)

scpModel <- scpModelWorkflow(daSimData, formula = ~ 1 + cell_type + condition, verbose = TRUE)
scpRes <- scpDifferentialAnalysis(
  scpSIM,
  contrasts = list(c("condition", "A", "B"))
)[[1]]

L <- makeContrast("conditionB=0", parameterNames = c("conditionB"))
msqModel <- suppressWarnings(msqrob(daSimData, formula = ~ cell_type + condition + (1 | patient_id)))
msqRes <- rowData(hypothesisTest(object = msqSIM, contrast = L))$conditionB

pbSIM <- suppressWarnings(msqrob(daSimDataAggregated, ~ 1 + condition + cell_type))
pbRes <- rowData(hypothesisTest(object = pbSIM, contrast = L))$conditionB

createCobraData <- function(modRes, rowdata) {
  data.frame(
    sid = modRes$sid,
    pval = modRes$p_val,
    padj.loc = modRes$p_adj.loc,
    padj.glb = modRes$p_adj.glb,
    is_de = rowdata$is_de,
    method = method_name
  )
}

perf_metrics <- function(metadata, selectedData, selectedNC, selectedProps, padjLoc) {
  subset <- metadata %>%
    filter(data == selectedData, NC == selectedNC, props == selectedProps, model != "MMdream2")

  df_list <- lapply(seq_len(nrow(subset)), function(i) {
    file <- subset$path[i]
    print(file)
    df <- readRDS(file) %>%
      mutate(E = (sim_mean.A + sim_mean.B) / 2) %>%
      filter(E > 0.1) %>%
      mutate(sid = paste(gene, cluster_id, sep = "_"))

    method_name <- paste(subset$model[i], subset$aggregationType[i], subset$valueType[i], sep = "_")

    # Return a data frame with method column names
    data.frame(
      sid = df$sid,
      pval = df$p_val,
      padj.loc = df$p_adj.loc,
      padj.glb = df$p_adj.glb,
      is_de = df$is_de,
      method = method_name
    )
  })

  combined_df <- bind_rows(df_list)

  pval_df <- combined_df %>%
    select(sid, method, pval) %>%
    pivot_wider(names_from = method, values_from = pval) %>%
    column_to_rownames("sid")

  padj.loc_df <- combined_df %>%
    select(sid, method, padj.loc) %>%
    pivot_wider(names_from = method, values_from = padj.loc) %>%
    column_to_rownames("sid")
  padj.glb_df <- combined_df %>%
    select(sid, method, padj.glb) %>%
    pivot_wider(names_from = method, values_from = padj.glb) %>%
    column_to_rownames("sid")

  truth_df <- combined_df %>%
    select(sid, is_de) %>%
    distinct() %>%
    column_to_rownames("sid")

  cobdata <- COBRAData(
    padj = as.data.frame(if (padjLoc) padj.loc_df else padj.glb_df),
    truth = as.data.frame(truth_df)
  )
  perf <- calculate_performance(cobdata,
                                binary_truth = "is_de",
                                aspects = "fdrtpr",
                                maxsplit = Inf,
                                thrs = c(0.01, 0.05, 0.1, 0.2)
  ) %>%
    fdrtpr() %>%
    mutate(thr = as.numeric(sub("thr", "", thr)))
  return(perf)
}
