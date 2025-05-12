library(scp)
library(tidyverse)
library(UpSetR)
library(muscat)
library(msqrob2)

sceBalanced <- readRDS("data/simulatedData/simulatedData_sampling_balanced_leduc.rds")

sceRes <- scpModelWorkflow(sceBalanced, ~ Set + Mock)

daRes <- scpDifferentialAnalysis(
    sceRes,
    contrasts = list(c("Mock", "mock2", "mock1"))
)[[1]] %>%
    as.data.frame() %>%
    filter(padj <= 0.05)

resNames <- daRes$feature

trueNames <- row.names(rowData(sceRes)[rowData(sceRes)$is_DA_mock1_vs_mock2, ])

#####

aggregate_se <- function(se, group_by_cols, fun = mean) {
    if (!all(group_by_cols %in% colnames(colData(se)))) {
        stop("Some specified group_by_cols do not exist in colData.")
    }

    colData(se)$Group <- apply(colData(se)[, group_by_cols, drop = FALSE], 1, paste, collapse = "_")
    group_factor <- factor(colData(se)$Group)
    assay_matrix <- assay(se)

    aggregated_assay <- do.call(cbind, lapply(split(seq_along(group_factor), group_factor),
                                              function(idx) apply(assay_matrix[, idx, drop = FALSE], 1, fun, na.rm = TRUE)))

    new_colData <- colData(se) %>%
        as.data.frame() %>%
        group_by(Group) %>%
        summarise(across(all_of(group_by_cols), dplyr::first), .groups = "drop") %>%
        as.data.frame()

    rownames(new_colData) <- new_colData$Group

    new_se <- SummarizedExperiment(assays = list(counts = aggregated_assay), colData = new_colData)

    return(new_se)
}

aggSE <- aggregate_se(sceBalanced, c("Set", "Mock"))

aggSEmd <- scpModelWorkflow(aggSE, ~ Set + Mock)
aggSEmd2 <- scpModelWorkflow(aggSE, ~ Mock)

daAggRes <- scpDifferentialAnalysis(
    aggSEmd,
    contrasts = list(c("Mock", "mock2", "mock1"))
)[[1]] %>%
    as.data.frame() %>%
    filter(padj <= 0.05)

daAggRes2 <- scpDifferentialAnalysis(
    aggSEmd2,
    contrasts = list(c("Mock", "mock2", "mock1"))
)[[1]] %>%
    as.data.frame() %>%
    filter(padj <= 0.05)

resAggNames <- daAggRes$feature
resAggNames2 <- daAggRes2$feature


listUpSet <- list("model" = resNames,
                  "aggModel" = resAggNames,
                  "real" = trueNames)
upset(fromList(listUpSet), order.by = "freq")
saveRDS(listUpSet, "dataOutput/pseudoBulk/listUpSet.rds")

listUpSet <- readRDS("dataOutput/pseudoBulk/listUpSet.rds")

upset(UpSetR::fromList(listUpSet))
