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

generateDesignTMT <- function(nRuns){
    base <- .generateDesignBaseTMT(nRuns)
    meta <- .generateDesignMetaTMT(nRuns)

    cbind(base, meta)
}

.generateQuantBaseTMT <- function(nRuns, nFeatures, naRate) {
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
    metaRef <- read.csv(file = "data/refTMTQuantTable.csv", row.names = 1)
    metaCols <- grep("Reporter.intensity",
                     names(metaRef),
                     invert = TRUE,
                     value = TRUE)
    metaCols <- setdiff(metaCols, c("Raw.file", "uid"))
    metaRef <- metaRef[, metaCols]

    metaRef[sample(1:nrow(metaRef), nFeatures*nRuns, replace = TRUE), ]
}

.generateDesignBaseTMT <- function(nRuns) {
    x_column <- paste0("run", rep(1:nRuns, each = 16))
    reporter_ion_column <- rep(paste0("Reporter.ion.", 1:16), times = nRuns)

    data.frame(
        runCol = x_column,
        quantCols = reporter_ion_column
    )
}

.generateDesignMetaTMT <- function(nRuns) {
    # should maybe be created manually to avoid inconstistencies with design
    metaRef <- read.csv(file = "data/refTMTDesignTable.csv", row.names = 1)
    metaRef <- metaRef[sample(1:nrow(metaRef), nRuns * 16, replace = TRUE),]

    subset(metaRef, select = -c(runCol, quantCols))
}
