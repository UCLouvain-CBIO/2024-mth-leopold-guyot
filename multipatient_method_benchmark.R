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
