library(SummarizedExperiment)
library(QFeatures)
library(tidyverse)
# load model

mod <- readRDS("dataOutput/slavovModels/pepMod.RDS")
mod <- mod[!is.na(mod)]

get_patient_effect_variance <- function(model) {
    patient_coefs <- grep("^patient_id", names(model$coefficients[, "Estimate"]), value = TRUE)
    var(c(0, model$coefficients[patient_coefs, "Estimate"])) # add 0 for the reference level
}

simulate_peptide <- function(model_info, metadata_df, synthetic_patient_effects = NULL) {
    coefs <- model_info$coefficients[, "Estimate"]
    #sigma <- model_info$sigma
    sigma <- 0.1
    #intercept <- coefs["(Intercept)"]
    intercept <- 0
    predicted <- rep(intercept, nrow(metadata_df))

    predicted <- predicted + sapply(metadata_df$cell_type, function(ct) {
        coef_name <- paste0("cell_type_", ct)
        if (coef_name %in% names(coefs)) coefs[coef_name] else 0
    })

    predicted <- predicted + sapply(metadata_df$patient_id, function(pid) {
        synthetic_patient_effects[[pid]]
    })

    intensity <- rnorm(nrow(metadata_df), mean = predicted, sd = sigma)
    return(intensity)
}


simulate_peptide_data <- function(mod,
                                  n_cells_per_comb = 10,
                                  cell_types,
                                  n_synthetic_patients = 6) {

    synthetic_patient_ids <- paste0("SynthP", seq_len(n_synthetic_patients))

    synthetic_meta <- expand.grid(
        cell_type = cell_types,
        patient_id = synthetic_patient_ids,
        stringsAsFactors = FALSE
    )
    synthetic_meta <- synthetic_meta[rep(seq_len(nrow(synthetic_meta)), each = n_cells_per_comb), ]
    synthetic_meta$cell_id <- paste0("SynCell_", seq_len(nrow(synthetic_meta)))

    sim_results <- lapply(names(mod), function(peptide) {
        model_info <- mod[[peptide]]

        # Estimate variance of patient effects
        patient_coefs <- grep("^patient_id", names(model_info$coefficients[, "Estimate"]), value = TRUE)
        var_pat <- get_patient_effect_variance(model_info)

        # Simulate effects for new patients
        synthetic_effects <- rnorm(n_synthetic_patients, mean = 0, sd = sqrt(var_pat))
        names(synthetic_effects) <- synthetic_patient_ids

        intensities <- simulate_peptide(model_info, synthetic_meta, synthetic_patient_effects = synthetic_effects)

        data.frame(
            cell_id = synthetic_meta$cell_id,
            peptide = peptide,
            intensity = intensities,
            cell_type = synthetic_meta$cell_type,
            patient_id = synthetic_meta$patient_id,
            stringsAsFactors = FALSE
        )
    })

    simData <- do.call(rbind, sim_results)

    simQuant <- simData[, c("cell_id", "peptide", "intensity")]
    simQuant <- pivot_wider(simQuant, names_from = "cell_id", values_from = "intensity")

    simColData <- simData[, c("cell_id", "cell_type", "patient_id")] %>%
        group_by(cell_id) %>%
        summarize("cell_type" = unique(cell_type), "patient_id" = unique(patient_id)) %>%
        as("DataFrame")
    rownames(simColData) <- simColData$cell_id
    rownames(simQuant) <- simQuant$peptide
    simse <- readSummarizedExperiment(simQuant, quantCols = 2:ncol(simQuant))
    colData(simse) <- simColData
    return(simse)
}

add_patient_group_shift_SE <- function(se,
                                       group_a,
                                       group_b,
                                       shift = 1,
                                       sd = 0.5,
                                       ratio = 0.1,
                                       seed = 123) {
    set.seed(seed)

    peptides <- rownames(se)
    n_peptides <- length(peptides)

    n_shift <- ceiling(n_peptides * ratio)
    shifted_peptides <- sample(peptides, n_shift)

    shift_values <- rnorm(n_shift, mean = shift, sd = sd)
    names(shift_values) <- shifted_peptides

    patient_ids <- colData(se)$patient_id
    group_b_cells <- colnames(se)[patient_ids %in% group_b]

    assay_mat <- assay(se)
    for (pep in shifted_peptides) {
        assay_mat[pep, group_b_cells] <- assay_mat[pep, group_b_cells] + shift_values[pep]
    }

    assay(se) <- assay_mat

    rowData(se)$shifted <- rownames(se) %in% shifted_peptides
    rowData(se)$shift_value <- NA
    rowData(se)$shift_value[match(shifted_peptides, rownames(se))] <- shift_values

    colData(se)$condition <- "A"
    colData(se)[group_b_cells, "condition"] <- "B"
    return(se)
}
