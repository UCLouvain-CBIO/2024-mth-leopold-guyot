library(QFeatures)
library(scp)
library(SingleCellExperiment)
library(tidyverse)
# Generate TMT PSM data using the TMT-18 leduc dataset provided as base
# scpdata::leduc2022_pSCoPE()
generateTMTPSM <- function(base, nCell) {
    base <- base[, , -(135:138)]
    nRun <- nCell %/% 18
    nAssays <- length(base)
    if (nRun <= nAssays) {
        sampledAssays <- sample(seq_len(nAssays), nRun, replace = FALSE)
    } else {
        sampledAssays <- sample(seq_len(nAssays), nRun, replace = TRUE)
    }
    new_qfeatures <- QFeatures()
    psmCounter <- 0
    for (i in seq_len(nRun)) {
        assay_idx <- sampledAssays[i]
        original_se <- getWithColData(base, assay_idx)
        newSampleNames <- paste0("run_", i, "_RI", seq_len(18))
        newAssay <- assay(original_se)
        newFeaturesNames <- paste0(
            "PSM_",
            seq(
                from = psmCounter + 1,
                length.out = nrow(newAssay)
            )
        )
        colnames(newAssay) <- newSampleNames
        rownames(newAssay) <- newFeaturesNames

        noise <- pmin(matrix(rnorm(n = length(newAssay), mean = 0, sd = 5),
                             nrow = nrow(newAssay),
                             ncol = ncol(newAssay)
        ), 0)
        noisyNewAssay <- newAssay + noise

        newColData <- colData(original_se)
        newColData$Set <- paste0("run_", i)
        newColData$Channel <- as.character(newColData$Channel)
        newColData$IsolationTimeStamp <- as.character(newColData$IsolationTimeStamp)
        newColData$filterBench <- rnorm(nrow(newColData), mean = 1, sd = 1)
        rownames(newColData) <- newSampleNames

        newRowData <- rowData(original_se)
        rownames(newRowData) <- newFeaturesNames
        newRowData$filterBench <- rnorm(nrow(newRowData), mean = 1, sd = 1)

        psmCounter <- psmCounter + nrow(newRowData)
        noisy_se <- SingleCellExperiment(
            assays = list(noisyNewAssay),
            rowData = newRowData,
            colData = newColData
        )
        new_qfeatures <- addAssay(new_qfeatures, noisy_se, name = paste0("run_", i))
    }
    return(new_qfeatures)
}
# using a QFeatures that only contains PSM assays, aggregate and joins these assays
generateTMTPeptides <- function(TMTPSM) {
    removeDuplicates <- function(x) {
        apply(x, 2, function(xx) xx[which(!is.na(xx))[1]])
    }
    peptidesAssays <- paste0("peptides_", names(TMTPSM))
    aggregated <- aggregateFeaturesOverAssays(
        object = TMTPSM,
        i = names(TMTPSM),
        fcol = "modseq",
        name = peptidesAssays,
        fun = removeDuplicates
    )
    ## Generate a list of DataFrames with the information to modify
    ppMap <- rbindRowData(aggregated, i = grep("^pep", names(aggregated))) %>%
        data.frame() %>%
        group_by(modseq) %>%
        ## The majority vote happens here
        mutate(
            Leading.razor.protein.symbol =
                names(sort(table(Leading.razor.protein),
                    decreasing = TRUE
                ))[1]
        ) %>%
        select(modseq, Leading.razor.protein.symbol) %>%
        filter(!duplicated(modseq, Leading.razor.protein.symbol))
    consensus <- lapply(peptidesAssays, function(i) {
        ind <- match(rowData(aggregated[[i]])$modseq, ppMap$modseq)
        DataFrame(
            Leading.razor.protein.symbol =
                ppMap$Leading.razor.protein.symbol[ind]
        )
    })
    ## Name the list
    names(consensus) <- peptidesAssays
    ## Modify the rowData
    rowData(aggregated) <- consensus
    aggregated <- infIsNA(aggregated, i = peptidesAssays)
    aggregated <- zeroIsNA(aggregated, i = peptidesAssays)
    joinAssays(aggregated, i = peptidesAssays, name = "peptides")
}

generateTMTProteins <- function(TMTPeptides) {
    aggregateFeatures(
        object = TMTPeptides,
        i = "peptides",
        fcol = "Leading.razor.protein.symbol",
        name = "proteins",
        fun = matrixStats::colMedians,
        na.rm = TRUE
    )
}

# using scpdata::brunner2022()
generateLFPSM <- function(base, nCell) {
    base <- base[, , -435]
    return(NULL) #WIP
}

# Generate QFeatures with 4 different variable:
#   - Number of cells by a assay
#   - Number of features (PSM) by assay
#   - Number of assay
#   - Number of cols in rowData
#
# Using a LF dataset as base (scpdata::brunner2022())
generate4VarData <- function(base,
                             nCell,
                             nFeat,
                             nAssay,
                             nCol,
                             seed = 123) {
    set.seed(seed)
    sampledAssays <- sample(length(base) - 1, nAssay, replace = TRUE)

    new_qfeatures <- QFeatures()
    i <- 0
    for (assay in sampledAssays) {
        i <- i + 1
        sampledFeatures <- sample(nrow(base[[assay]]), nFeat, replace = TRUE)
        newAssay <- assay(base[[assay]])
        newAssay <- as.matrix(newAssay[sampledFeatures, rep(1, nCell)])

        newRowData <- rowData(base[[assay]])
        sampledRowDataCols <- sample(ncol(newRowData), nCol, replace = TRUE)
        newRowData <- newRowData[sampledFeatures, sampledRowDataCols]

        newColData <- colData(getWithColData(base, i = assay))
        newColData <- newColData[rep(1, nCell), ]

        sce <- SingleCellExperiment(
            assays = list(newAssay),
            rowData = newRowData,
            colData = newColData
        )
        rownames(sce) <- paste0("run_", i, "_PSM", 1:nrow(sce))
        colnames(sce) <- paste0("run_", i, "_Cell", 1:ncol(sce))
        new_qfeatures <- addAssay(new_qfeatures, sce, name = paste0("run_", i))
    }
    for (assay in names(new_qfeatures)) colData(new_qfeatures[[assay]]) <- NULL
    return(new_qfeatures)
}
