### Replicate Data Function ###
replicateData <- function(nFeaturesRange, nRunsRange, naRateRange, nReplicates, type, folder) {
    if (!dir.exists(folder)) {
        dir.create(folder)
    }
    generateFunction <- switch(type,
        "TMT" = generateTMTtoFile,
        "LFDIA" = generateLFDIAtoFile,
        "plexDIA" = generatePlexDIAtoFile,
        stop("Unknown type in replicateData: ", type)
    )
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
                }
            }
        }
    }
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

    cat(
        "Generated", type, "data for nFeatures =", nFeatures,
        "nRuns =", nRuns, "naRate =", naRate, "in", subfolder, "\n"
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
    res <- data.frame(
        PSM = paste0("PSM", rep(1:nFeatures, nRuns)),
        run = paste0("run", rep(1:nRuns, each = nFeatures)),
        `Reporter.ion.1` = runif(nRuns * nFeatures, max = 10000),
        `Reporter.ion.2` = runif(nRuns * nFeatures, max = 10000),
        `Reporter.ion.3` = runif(nRuns * nFeatures, max = 10000),
        `Reporter.ion.4` = runif(nRuns * nFeatures, max = 10000),
        `Reporter.ion.5` = runif(nRuns * nFeatures, max = 10000),
        `Reporter.ion.6` = runif(nRuns * nFeatures, max = 10000),
        `Reporter.ion.7` = runif(nRuns * nFeatures, max = 10000),
        `Reporter.ion.8` = runif(nRuns * nFeatures, max = 10000),
        `Reporter.ion.9` = runif(nRuns * nFeatures, max = 10000),
        `Reporter.ion.10` = runif(nRuns * nFeatures, max = 10000),
        `Reporter.ion.11` = runif(nRuns * nFeatures, max = 10000),
        `Reporter.ion.12` = runif(nRuns * nFeatures, max = 10000),
        `Reporter.ion.13` = runif(nRuns * nFeatures, max = 10000),
        `Reporter.ion.14` = runif(nRuns * nFeatures, max = 10000),
        `Reporter.ion.15` = runif(nRuns * nFeatures, max = 10000),
        `Reporter.ion.16` = runif(nRuns * nFeatures, max = 10000)
    )
    reporter_cols <- grep("Reporter.ion", colnames(res))
    for (col in reporter_cols) {
        n <- nrow(res)
        # Randomly sample naRate of the indices
        zero_indices <- sample(1:n, size = floor(naRate * n))
        res[zero_indices, col] <- 0
    }

    return(res)
}

.generateQuantMetaTMT <- function(nRuns, nFeatures) {
    set.seed(123)
    metaRef <- read.csv(file = "data/refTMTQuantTable.csv", row.names = 1)
    metaCols <- grep("Reporter.intensity",
        names(metaRef),
        invert = TRUE,
        value = TRUE
    )
    metaCols <- setdiff(metaCols, c("Raw.file", "uid"))
    metaRef <- metaRef[, metaCols]

    metaRef[sample(1:nrow(metaRef), nFeatures * nRuns, replace = TRUE), ]
}

.generateDesignBaseTMT <- function(nRuns) {
    set.seed(123)
    x_column <- paste0("run", rep(1:nRuns, each = 16))
    reporter_ion_column <- rep(paste0("Reporter.ion.", 1:16), times = nRuns)
    SampleType <- rep(c("Carrier", rep("Macrophage", 5), rep("Monocyte", 5), rep("Blank", 5)), nRuns)

    data.frame(
        runCol = x_column,
        quantCols = reporter_ion_column,
        SampleType = SampleType
    )
}
