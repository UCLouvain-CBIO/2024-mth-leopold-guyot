library(visreg)
library(scpdata)
library(microbenchmark)

source("R/generate_data.R")

leduc <- scpdata::leduc2022_pSCoPE()

PSM1000 <- generateTMTPSM(leduc, 1000)


.fillEmptyExpsCopy <- function(exps, subr) {
    if (!any(names(subr) %in% names(exps))) {
        stop("No matching experiment names in subset list", call. = FALSE)
    }
    if (!all(names(exps) %in% names(subr))) {
        outnames <- setdiff(names(exps), names(subr))
        names(outnames) <- outnames
        subr <- c(subr, lapply(outnames, function(x) NULL))
    }
    subr[names(exps)]
}

.matchReorderSubCopy <- function(assayMap, identifiers) {
    positions <- unlist(lapply(identifiers, function(ident) {
        which(!is.na(match(assayMap[["primary"]], ident)))
    }))
    assayMap[positions, ]
}

subsetByColDataCopy <- function(x, y) {
    coldata <- colData(x)
    if (length(y) > nrow(coldata)) {
        stop("subscript vector 'j' in 'mae[i, j, k]' is out-of-bounds",
            call. = FALSE
        )
    }
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
    BiocBaseUtils::setSlots(x,
        ExperimentList = newSubset, colData = newcoldata,
        sampleMap = newMap, check = FALSE
    )
}
profvis::profvis(subsetByColDataCopy(PSM1000, PSM1000$filterBench > 1))

x <- PSM1000
y <- PSM1000$filterBench > 1

coldata <- colData(x)
if (length(y) > nrow(coldata)) {
    stop("subscript vector 'j' in 'mae[i, j, k]' is out-of-bounds",
        call. = FALSE
    )
}
newcoldata <- coldata[y, , drop = FALSE]
listMap <- mapToList(sampleMap(x), "assay")

# Bottleneck
profvis::profvis(
    lapply(listMap, function(elementMap, keepers) {
        .matchReorderSubCopy(elementMap, keepers)
    }, keepers = rownames(newcoldata))
)

assayMap <- listMap[[1]]
identifiers <- rownames(newcoldata)

profvis::profvis({
    positions <- unlist(lapply(identifiers, function(ident) {
        which(!is.na(match(assayMap[["primary"]], ident)))
    }))
    assayMap[positions, ]
})

profvis::profvis({
    positions <- unlist(lapply(identifiers, function(ident) {
        which(!is.na(match(assayMap[["primary"]], ident)))
    }))
    assayMap[positions, ]
})

profvis::profvis({
    positions <- lapply(identifiers, function(ident) {
        which(!is.na(match(assayMap[["primary"]], ident)))
    })
})

# I think the problem is that we iterate over all the samples (across all the assays)
# But we know that only the sample that belong to this assay can be found on this assay
# If we use only the sample names coming from the assay the runtime is now so little that
# it cannot be profiled
# This is also why the complexity is not linear bc as we increase the number of assay:
# - We increase the number of time we need to call .matchReorderSub
# - We increase the iteration of the inner loop of .matchReorderSub
# The complexity is in O(n²)
identifiers_run_1 <- grep("run_1_", identifiers)

profvis::profvis({
    positions <- lapply(identifiers, function(ident) {
        which(!is.na(match(assayMap[["primary"]], ident)))
    })
    profvis::pause(1)
})



subsetByColDataCopy2 <- function(x, y) {
    coldata <- colData(x)
    if (length(y) > nrow(coldata)) {
        stop("subscript vector 'j' in 'mae[i, j, k]' is out-of-bounds",
            call. = FALSE
        )
    }
    newcoldata <- coldata[y, , drop = FALSE]
    listMap <- mapToList(sampleMap(x), "assay")
    assaysName <- names(listMap)
    listMap <- lapply(assaysName, function(assayName) {
        .matchReorderSubCopy(
            listMap[[assayName]],
            rownames(newcoldata[newcoldata$Set == assayName, ])
        )
    })
    names(listMap) <- assaysName
    newMap <- listToMap(listMap, fill = FALSE)
    columns <- lapply(listMap, function(mapChunk) {
        mapChunk[, "colname", drop = TRUE]
    })
    columns <- .fillEmptyExpsCopy(experiments(x), columns)
    newSubset <- Map(function(x, j) {
        x[, j, drop = FALSE]
    }, x = experiments(x), j = columns)
    newSubset <- ExperimentList(newSubset)
    BiocBaseUtils::setSlots(x,
        ExperimentList = newSubset, colData = newcoldata,
        sampleMap = newMap, check = FALSE
    )
}

profvis::profvis(subsetByColDataCopy2(PSM1000, PSM1000$filterBench > 1))


subsetByColDataCopy3 <- function(x, y) {
    coldata <- colData(x)
    if (length(y) > nrow(coldata)) {
        stop("subscript vector 'j' in 'mae[i, j, k]' is out-of-bounds",
            call. = FALSE
        )
    }
    newcoldata <- coldata[y, , drop = FALSE]
    listMap <- mapToList(sampleMap(x), "assay")
    listMap <- lapply(listMap, function(elementMap) {
        .matchReorderSubCopy(
            elementMap,
            intersect(
                rownames(newcoldata),
                elementMap$primary
            )
        )
    })
    newMap <- listToMap(listMap, fill = FALSE)
    columns <- lapply(listMap, function(mapChunk) {
        mapChunk[, "colname", drop = TRUE]
    })
    columns <- .fillEmptyExpsCopy(experiments(x), columns)
    newSubset <- Map(function(x, j) {
        x[, j, drop = FALSE]
    }, x = experiments(x), j = columns)
    newSubset <- ExperimentList(newSubset)
    BiocBaseUtils::setSlots(x,
        ExperimentList = newSubset, colData = newcoldata,
        sampleMap = newMap, check = FALSE
    )
}

profvis::profvis(subsetByColDataCopy3(PSM1000, PSM1000$filterBench > 1))

PSM4000 <- generateTMTPSM(leduc, 4000)
microbenchmark(subsetByColDataCopy(PSM4000, PSM4000$filterBench > 1),
    subsetByColDataCopy2(PSM4000, PSM4000$filterBench > 1),
    subsetByColDataCopy3(PSM4000, PSM4000$filterBench > 1),
    times = 3
)
