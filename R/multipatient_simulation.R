library(SummarizedExperiment)
library(QFeatures)
library(tidyverse)
# load model

mod <- readRDS("dataOutput/slavovModels/pepMod.RDS")
mod <- mod[!is.na(mod)]

simulate_peptide <- function(model_info, metadata_df) {
    coefs <- model_info$coefficients[, "Estimate"]
    sigma <- model_info$sigma

    intercept <- coefs["(Intercept)"]

    predicted <- rep(intercept, nrow(metadata_df))

    predicted <- predicted + sapply(metadata_df$cell_type, function(ct) {
        coef_name <- paste0("cell_type_", ct)
        if (coef_name %in% names(coefs)) coefs[coef_name] else 0
    })

    predicted <- predicted + sapply(metadata_df$patient_id, function(pid) {
        coef_name <- paste0("patient_id", pid)
        if (coef_name %in% names(coefs)) coefs[coef_name] else 0
    })

    # Add Gaussian noise
    intensity <- rnorm(nrow(metadata_df), mean = predicted, sd = sigma)
    return(intensity)
}

simulate_peptide_data <- function(mod,
                                  n_cells_per_comb = 10,
                                  cell_types,
                                  patient_ids) {
    synthetic_meta <- expand.grid(
        cell_type = cell_types,
        patient_id = patient_ids,
        stringsAsFactors = FALSE
    )

    synthetic_meta <- synthetic_meta[rep(seq_len(nrow(synthetic_meta)), each = n_cells_per_comb), ]
    synthetic_meta$cell_id <- paste0("SynCell_", seq_len(nrow(synthetic_meta)))

    sim_results <- lapply(names(mod), function(peptide) {
        intensities <- simulate_peptide(mod[[peptide]], synthetic_meta)
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


    simse <- readSummarizedExperiment(simQuant, quantCols = 2:ncol(simQuant))
    colData(simse) <- simColData
    return(simse)
}

simulate_peptide_data(mod, n_cells_per_comb = 500,
                      cell_types = c("CT1", "lowerresCD8T", "lowerresmonocyte", "lowerresNK"),
                      patient_ids = c("P1", "3430861_d0", "3431438_d0", "3431452_d0")  )
