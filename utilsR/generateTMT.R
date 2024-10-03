### Generate random 16-TMT data of desired size ###

### Package Loading ###

library(scp)

# Params:
# - `nRuns` (`integer`) the number of runs
# - `nPSM` (`integer`) the number of PSM by run
# - `naRate` (`numeric`) the rate of missing values

# Output: two tables, one with the quantitative data and
# one with the experimental design

generateTMT <- function(nRuns, nPSM, naRate) {
    # Generate the quantitative data
    generateQuantTMT(nRuns, nPSM, naRate)

    # WIP
}

generateQuantTMT <- function(nRuns, nPSM, naRate) {
    base <- .generateQuantBaseTMT(nRuns, nPSM, naRate)
    meta <- .generateQuantMetaTMT(nRuns, nPSM)

    cbind(base, meta)
}

generateDesignTMT <- function(nRuns){
    # WIP
}

.generateQuantBaseTMT <- function(nRuns, nPSM, naRate) {
    res <- data.frame(
        PSM = paste0("PSM", rep(1:nPSM, nRuns)),
        run = rep(1:nRuns, each = nPSM),
        `Reporter.ion.1` = runif(nRuns * nPSM, max = 10000),
        `Reporter.ion.2` = runif(nRuns * nPSM, max = 10000),
        `Reporter.ion.3` = runif(nRuns * nPSM, max = 10000),
        `Reporter.ion.4` = runif(nRuns * nPSM, max = 10000),
        `Reporter.ion.5` = runif(nRuns * nPSM, max = 10000),
        `Reporter.ion.6` = runif(nRuns * nPSM, max = 10000),
        `Reporter.ion.7` = runif(nRuns * nPSM, max = 10000),
        `Reporter.ion.8` = runif(nRuns * nPSM, max = 10000),
        `Reporter.ion.9` = runif(nRuns * nPSM, max = 10000),
        `Reporter.ion.10` = runif(nRuns * nPSM, max = 10000),
        `Reporter.ion.11` = runif(nRuns * nPSM, max = 10000),
        `Reporter.ion.12` = runif(nRuns * nPSM, max = 10000),
        `Reporter.ion.13` = runif(nRuns * nPSM, max = 10000),
        `Reporter.ion.14` = runif(nRuns * nPSM, max = 10000),
        `Reporter.ion.15` = runif(nRuns * nPSM, max = 10000),
        `Reporter.ion.16` = runif(nRuns * nPSM, max = 10000)
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

.generateQuantMetaTMT <- function(nRuns, nPSM) {
    metaRef <- read.csv(file = "data/refTMTQuantTable.csv")
    metaCols <- grep("Reporter.intensity",
                     names(metaRef),
                     invert = TRUE,
                     value = TRUE)
    metaCols <- setdiff(metaCols, c("Raw.file", "uid"))
    metaRef <- metaRef[, metaCols]
    res <- metaRef[sample(1:nrow(metaRef), nPSM*nRuns, replace = TRUE), ]

    return(res)
}
