library(visreg)
library(scpdata)

source("R/generate_data.R")

leduc <- scpdata::leduc2022_pSCoPE()

PSM2000 <- generateTMTPSM(leduc, 2000)


.fillEmptyExpsCopy <- function (exps, subr)
{
    if (!any(names(subr) %in% names(exps)))
        stop("No matching experiment names in subset list", call. = FALSE)
    if (!all(names(exps) %in% names(subr))) {
        outnames <- setdiff(names(exps), names(subr))
        names(outnames) <- outnames
        subr <- c(subr, lapply(outnames, function(x) NULL))
    }
    subr[names(exps)]
}

.matchReorderSubCopy <- function(assayMap, identifiers)
{
    positions <- unlist(lapply(identifiers, function(ident) {
        which(!is.na(match(assayMap[["primary"]], ident)))
    }))
    assayMap[positions, ]
}

subsetByColDataCopy <- function(x, y)
{
    coldata <- colData(x)
    if (length(y) > nrow(coldata))
        stop("subscript vector 'j' in 'mae[i, j, k]' is out-of-bounds",
             call. = FALSE)
    newcoldata <- coldata[y, , drop = FALSE]
    listMap <- mapToList(sampleMap(x), "assay")
    listMap <- lapply(listMap, function(elementMap, keepers) {
        .matchReorderSubCopy(elementMap, keepers)
    }, keepers = rownames(newcoldata))
    newMap <- listToMap(listMap, fill = FALSE)
    columns <- lapply(listMap, function(mapChunk) {
        mapChunk[, "colname", drop = TRUE]
    })
    columns <- .fillEmptyExpsCopy(experiments(x), columns)
    newSubset <- Map(function(x, j) {
        x[, j, drop = FALSE]
    }, x = experiments(x), j = columns)
    newSubset <- ExperimentList(newSubset)
    BiocBaseUtils::setSlots(x, ExperimentList = newSubset, colData = newcoldata,
                            sampleMap = newMap, check = FALSE)
}
profvis::profvis(subsetByColDataCopy(PSM2000, PSM2000$filterBench > 1))
