#' Generate Replicated Proteomic Data
#'
#' This function generates replicated proteomic datasets with varying numbers of features,
#' runs, and missing data rates. The generated datasets are saved to the specified folder
#' and can be created for TMT, LFDIA, or plexDIA data types.
#'
#' @param nFeaturesRange A numeric vector specifying the range of feature counts to generate.
#' @param nRunsRange A numeric vector specifying the range of run counts to generate.
#' @param naRateRange A numeric vector specifying the range of missing data rates (NA rates)
#'        to apply during data generation.
#' @param nReplicates An integer specifying the number of replicates to generate for each
#'        combination of parameters.
#' @param type A character string specifying the data type to generate. Accepted values
#'        are "TMT", "LFDIA", or "plexDIA".
#' @param folder A character string specifying the directory path where the generated datasets
#'        will be saved.
#'
#' @return No return value. This function generates files and saves them to the specified folder.
#'
#' @examples
#' # Generate replicated data for TMT with specified parameter ranges
#' # replicateData(nFeaturesRange = 50:100, nRunsRange = 5:10,
#' #               naRateRange = c(0, 0.1), nReplicates = 3,
#' #               type = "TMT", folder = "data_folder")
#'
#' @import progress
#'
#' @export
replicateData <- function(
        nFeaturesRange,
        nRunsRange,
        naRateRange,
        nReplicates,
        type,
        folder) {
    if (!dir.exists(folder)) {
        dir.create(folder)
    }
    generateFunction <- switch(type,
        "TMT" = generateTMTtoFile,
        "LFDIA" = generateLFDIAtoFile,
        "plexDIA" = generatePlexDIAtoFile,
        stop("Unknown type in replicateData: ", type)
    )
    steps <- length(nFeaturesRange) *
        length(nRunsRange) *
        length(naRateRange) *
        nReplicates

    cat("Starting", type, " Generation ...")

    pb <- progress_bar$new(
        format = paste0(
            type,
            " Generation: (:percent) [:bar] ",
            ":current/:total | Elapsed: :elapsed"
        ),
        total = steps,
        clear = FALSE,
        width = 80
    )
    pb$tick(0)

    for (nFeatures in nFeaturesRange) {
        for (nRuns in nRunsRange) {
            for (naRate in naRateRange) {
                for (nReplicate in 1:nReplicates) {
                    .createReplicate(
                        nFeatures = nFeatures,
                        nRuns = nRuns,
                        naRate = naRate,
                        nReplicate = nReplicate,
                        type = type,
                        folder = folder,
                        generateFunction = generateFunction
                    )
                    pb$tick()
                }
            }
        }
    }
    cat("TMT Generation Finished.\n")
}

#' Create a Single Replicate of Proteomic Data
#'
#' This internal helper function generates a single replicate of proteomic data
#' based on the specified parameters and saves it to the appropriate subfolder.
#'
#' @param nFeatures An integer specifying the number of features to generate in the dataset.
#' @param nRuns An integer specifying the number of runs or experiments to generate.
#' @param naRate A numeric value specifying the rate of missing data (NA rate) to apply.
#' @param nReplicate An integer indicating the replicate number for the given parameter combination.
#' @param type A character string specifying the type of data to generate ("TMT", "LFDIA", or "plexDIA").
#' @param folder A character string specifying the main folder where the dataset will be saved.
#' @param generateFunction A function that generates the dataset based on the type of data.
#'
#' @return No return value. This function generates files and saves them to the specified folder.
#'
#' @keywords internal
.createReplicate <- function(nFeatures,
    nRuns,
    naRate,
    nReplicate,
    type,
    folder,
    generateFunction) {
    # Define the subfolder name based on the current combination
    subfolder <- paste0(
        folder, "/",
        type,
        "_nFeatures_", nFeatures,
        "_nRuns_", nRuns,
        "_naRate_", naRate,
        "_replicate_", nReplicate
    )

    if (!dir.exists(subfolder)) {
        dir.create(subfolder)
    }

    # Generate the file paths for the files
    quantPath <- paste0(subfolder, "/quant.csv")
    designPath <- paste0(subfolder, "/design.csv")

    generateFunction(
        nRuns = nRuns,
        nFeatures = nFeatures,
        naRate = naRate,
        quantPath = quantPath,
        designPath = designPath
    )
}

#' Generate 16-TMT Quantitative Data and Save to File
#'
#' This function generates synthetic quantitative 16-TMT data along with experimental design data
#' and saves them as CSV files. The generated data includes random values with customizable
#' missing data rates.
#'
#' @param nRuns Integer. The number of runs to generate.
#' @param nFeatures Integer. The number of features per run.
#' @param naRate Numeric. The rate of missing values (0-1) to introduce into the quantitative data.
#' @param quantPath Character. The file path where the quantitative data will be saved.
#' @param designPath Character. The file path where the design data will be saved.
#'
#' @return No return value. This function generates and saves files.
#'
#' @examples
#' # Generate data for 5 runs, 100 features per run, with 10% missing data
#' # generateTMTtoFile(nRuns = 5, nFeatures = 100, naRate = 0.1,
#' #                  quantPath = "quant.csv", designPath = "design.csv")
#'
#' @export
generateTMTtoFile <- function(nRuns, nFeatures, naRate, quantPath, designPath) {
    quantTable <- .generateQuantTMT(nRuns, nFeatures, naRate)
    designTable <- .generateDesignTMT(nRuns)

    write.csv(quantTable, file = quantPath, row.names = FALSE)
    write.csv(designTable, file = designPath, row.names = FALSE)
}

#' Generate Quantitative Data for 16-TMT
#'
#' This function creates synthetic quantitative data for a 16-TMT experiment.
#' Data includes reporter ion intensities and various metadata columns.
#'
#' @param nRuns Integer. The number of experimental runs.
#' @param nFeatures Integer. The number of features per run.
#' @param naRate Numeric. The proportion of missing data to introduce.
#'
#' @return A data.frame containing the synthetic quantitative data for 16-TMT.
#'
#' @keywords internal
.generateQuantTMT <- function(nRuns, nFeatures, naRate) {
    base <- .generateQuantBaseTMT(nRuns, nFeatures, naRate)
    meta <- .generateQuantMetaTMT(nRuns, nFeatures)

    cbind(base, meta)
}

#' Generate Experimental Design Data for 16-TMT
#'
#' This function creates synthetic experimental design data for a 16-TMT experiment.
#'
#' @param nRuns Integer. The number of experimental runs.
#'
#' @return A data.frame containing the synthetic design data for 16-TMT.
#'
#' @keywords internal
.generateDesignTMT <- function(nRuns) {
    .generateDesignBaseTMT(nRuns)
}

#' Generate Base Quantitative Data for 16-TMT
#'
#' This function generates the core quantitative data, including reporter ion intensities
#' and feature-specific metadata, for a 16-TMT experiment.
#'
#' @param nRuns Integer. The number of runs in the experiment.
#' @param nFeatures Integer. The number of features per run.
#' @param naRate Numeric. The rate of missing data (0-1) to apply to the quantitative data.
#'
#' @return A data.frame containing the base quantitative data for 16-TMT.
#'
#' @keywords internal
.generateQuantBaseTMT <- function(nRuns, nFeatures, naRate) {
    set.seed(123)
    totalRows <- nRuns * nFeatures
    res <- data.frame(
        PSM = paste0("PSM", rep(1:nFeatures, nRuns)),
        run = paste0("run", rep(1:nRuns, each = nFeatures)),
        `Reporter.ion.1` = runif(totalRows, max = 150000),
        `Reporter.ion.2` = runif(totalRows, max = 1500),
        `Reporter.ion.3` = runif(totalRows, max = 1500),
        `Reporter.ion.4` = runif(totalRows, max = 1500),
        `Reporter.ion.5` = runif(totalRows, max = 1500),
        `Reporter.ion.6` = runif(totalRows, max = 1500),
        `Reporter.ion.7` = runif(totalRows, max = 1500),
        `Reporter.ion.8` = runif(totalRows, max = 1500),
        `Reporter.ion.9` = runif(totalRows, max = 1500),
        `Reporter.ion.10` = runif(totalRows, max = 1500),
        `Reporter.ion.11` = runif(totalRows, max = 1500),
        `Reporter.ion.12` = runif(totalRows, max = 1500),
        `Reporter.ion.13` = runif(totalRows, max = 1500),
        `Reporter.ion.14` = runif(totalRows, max = 1500),
        `Reporter.ion.15` = runif(totalRows, max = 1500),
        `Reporter.ion.16` = runif(totalRows, max = 1500)
    )
    reporter_cols <- grep("Reporter.ion", colnames(res))
    for (col in reporter_cols[-1]) {
        n <- nrow(res)
        # Randomly sample naRate of the indices
        zero_indices <- sample(1:n, size = floor(naRate * n))
        res[zero_indices, col] <- 0
    }

    # Add peptidesId
    proteinIds <- paste0("Peptide", sample(floor(0.5 * totalRows),
        totalRows,
        replace = TRUE
    ))
    res$peptidesId <- proteinIds

    # Add Leading.razor.protein column
    proteinIds <- paste0("Protein", sample(floor(0.2 * totalRows),
        totalRows,
        replace = TRUE
    ))
    res$Leading.razor.protein <- proteinIds

    # Add PEP column with random probabilities
    res$PEP <- pmin(pmax(rexp(totalRows, 30), 0), 1)
    # Add PIF column with random fractions

    pifValues <- rnorm(totalRows, mean = 0.95, sd = 0.25)
    pifValues <- pmin(pmax(pifValues, 0), 1) # Clip values to [0, 1]

    # Set 10% of PIF values to NA
    naPIF <- sample(1:totalRows, size = floor(0.1 * totalRows))
    pifValues[naPIF] <- NA
    res$PIF <- pifValues

    # Add Reverse column with 5% of '+' and rest empty
    res$Reverse <- ""
    reverse_indices <- sample(1:(totalRows), size = floor(0.05 * totalRows))
    res$Reverse[reverse_indices] <- "+"

    # Add Potential.contaminant column with 5% of '+' and rest empty
    res$Potential.contaminant <- ""
    contaminant_indices <- sample(1:(totalRows), size = floor(0.05 * totalRows))
    res$Potential.contaminant[contaminant_indices] <- "+"


    return(res)
}

#' Generate Additional Metadata for Quantitative Data in 16-TMT
#'
#' This function generates supplementary columns for the quantitative data in a 16-TMT experiment.
#'
#' @param nRuns Integer. The number of runs in the experiment.
#' @param nFeatures Integer. The number of features per run.
#'
#' @return A data.frame containing metadata for the quantitative data in 16-TMT.
#'
#' @keywords internal
.generateQuantMetaTMT <- function(nRuns, nFeatures) {
    set.seed(123)
    totalRows <- nRuns * nFeatures
    doubleCols <- replicate(85, rnorm(totalRows))
    charCols <- replicate(21, sample(LETTERS, totalRows, replace = TRUE))
    intCols <- replicate(3, sample(1:100, totalRows, replace = TRUE))
    logicalCols <- replicate(5, sample(c(TRUE, FALSE),
        totalRows,
        replace = TRUE
    ))

    data <- data.frame(doubleCols, charCols, intCols, logicalCols)

    # Set column names
    colnames(data) <- c(
        paste0("Double_", seq(85)),
        paste0("Char_", seq(21)),
        paste0("Int_", seq(3)),
        paste0("Logical_", seq(5))
    )

    return(data)
}

#' Generate Experimental Design Data for 16-TMT
#'
#' This function generates base columns for design data used in a 16-TMT experiment.
#'
#' @param nRuns Integer. The number of runs in the experiment.
#'
#' @return A data.frame containing the base design data for 16-TMT.
#'
#' @keywords internal
.generateDesignBaseTMT <- function(nRuns) {
    set.seed(123)
    x_column <- paste0("run", rep(1:nRuns, each = 16))
    reporter_ion_column <- rep(paste0("Reporter.ion.", 1:16), times = nRuns)
    SampleType <- rep(c(
        "Carrier", rep("Macrophage", 5),
        rep("Monocyte", 5), rep("Blank", 5)
    ), nRuns)

    data.frame(
        runCol = x_column,
        quantCols = reporter_ion_column,
        SampleType = SampleType
    )
}

#' Generate Additional Metadata for Experimental Design in 16-TMT
#'
#' This function generates supplementary columns for the design data in a 16-TMT experiment.
#'
#' @param nRuns Integer. The number of runs in the experiment.
#'
#' @return A data.frame containing metadata for the design data in 16-TMT.
#'
#' @keywords internal
.generateDesignMetaTMT <- function(nRuns) {
    set.seed(123)
    totalRows <- nRuns * 16 # 16-TMT
    doubleCols <- replicate(8, rnorm(totalRows))
    charCols <- replicate(10, sample(LETTERS, totalRows, replace = TRUE))

    data <- data.frame(doubleCols, charCols)

    # Set column names
    colnames(data) <- c(
        paste0("Double_", seq(8)),
        paste0("Char_", seq(10))
    )

    return(data)
}
