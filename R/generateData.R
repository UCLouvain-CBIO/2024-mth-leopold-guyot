### Replicate Data Function ###
replicateData <- function(nFeaturesRange,
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
        format = paste0(type,
                        " Generation: (:percent) [:bar] ",
                        ":current/:total | Elapsed: :elapsed"),
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


### Generate random 16-TMT data of desired size ###

### Package Loading ###

library(scp)

### Functions ###

generateTMTtoFile <- function(nRuns, nFeatures, naRate, quantPath, designPath) {
    quantTable <- generateQuantTMT(nRuns, nFeatures, naRate)
    designTable <- generateDesignTMT(nRuns)

    write.csv(quantTable, file = quantPath, row.names = FALSE)
    write.csv(designTable, file = designPath, row.names = FALSE)
}

generateQuantTMT <- function(nRuns, nFeatures, naRate) {
    base <- .generateQuantBaseTMT(nRuns, nFeatures, naRate)
    meta <- .generateQuantMetaTMT(nRuns, nFeatures)

    cbind(base, meta)
}

generateDesignTMT <- function(nRuns) {
    .generateDesignBaseTMT(nRuns)
}

.generateQuantBaseTMT <- function(nRuns, nFeatures, naRate) {
    set.seed(123)
    totalRows <- nRuns * nFeatures
    res <- data.frame(
        PSM = paste0("PSM", rep(1:nFeatures, nRuns)),
        run = paste0("run", rep(1:nRuns, each = nFeatures)),
        `Reporter.ion.1` = runif(totalRows,  max = 150000),
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
                                           totalRows, replace = TRUE))
    res$peptidesId <- proteinIds

    # Add Leading.razor.protein column
    proteinIds <- paste0("Protein", sample(floor(0.2 * totalRows),
                                            totalRows, replace = TRUE))
    res$Leading.razor.protein <- proteinIds

    # Add PEP column with random probabilities
    res$PEP <- pmin(pmax(rexp(totalRows, 30), 0), 1)
    # Add PIF column with random fractions

    pifValues <- rnorm(totalRows, mean = 0.95, sd = 0.25)
    pifValues <- pmin(pmax(pifValues, 0), 1)  # Clip values to [0, 1]

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

# This function create filler column for the quantitative data. The goal is to
# mimic the size of a typical rowData object in a TMT experience. We use here as
# reference the leduc dataset. This dataset has a rowData that contains 114
# variables. There are 85 `double`, 21 `character`, 3 `integer` and 5 `logical`.
# We will mimic this distribution in the filler columns.
# We generate 85 `double` columns, 21 `character` columns, 3 `integer` columns.

.generateQuantMetaTMT <- function(nRuns, nFeatures) {
    set.seed(123)
    totalRows <- nRuns * nFeatures
    doubleCols <- replicate(85, rnorm(totalRows))
    charCols <- replicate(21, sample(LETTERS, totalRows, replace = TRUE))
    intCols <- replicate(3, sample(1:100, totalRows, replace = TRUE))
    logicalCols <- replicate(5, sample(c(TRUE, FALSE),
                                       totalRows, replace = TRUE))

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

.generateDesignBaseTMT <- function(nRuns) {
    set.seed(123)
    x_column <- paste0("run", rep(1:nRuns, each = 16))
    reporter_ion_column <- rep(paste0("Reporter.ion.", 1:16), times = nRuns)
    SampleType <- rep(c("Carrier", rep("Macrophage", 5),
                        rep("Monocyte", 5), rep("Blank", 5)), nRuns)

    data.frame(
        runCol = x_column,
        quantCols = reporter_ion_column,
        SampleType = SampleType
    )
}

# This function create filler column for the design data. The goal is to
# mimic the size of a typical colData object in a TMT experience. We use here as
# reference the leduc dataset. This dataset has a colData that contains 18
# variables. There are 8 `double`, 10 `character`.
# We will mimic this distribution in the filler columns.
# We generate 8 `double` columns, 10 `character` columns.

.generateDesignMetaTMT <- function(nRuns) {
    set.seed(123)
    totalRows <- nRuns*16 # 16-TMT
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
