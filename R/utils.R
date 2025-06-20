library(rmarkdown)

renderRmarkdown <- function(rmdFile) {
    rmarkdown::render(file.path("Rmd", rmdFile),
                      output_format = "html_document",
                      output_file = sub(".rmd", ".html", rmdFile),
                      output_dir = "reports",
                      knit_root_dir = getwd()
    )
}

getColDataSize <- function(qfeatures) {
    object.size(qfeatures@colData)
}

getRowDataSize <- function(qfeatures) {
    sizes <- sapply(names(qfeatures), function(x) {
        object.size(rowData(qfeatures[[x]]))
    })

    sum(sizes)
}

getAssaySize <- function(qfeatures) {
    sizes <- sapply(names(qfeatures), function(x) {
        object.size(assay(qfeatures[[x]]))
    })

    sum(sizes)
}

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

    new_colData <- new_colData[colnames(aggregated_assay),]
    new_se <- SummarizedExperiment(assays = list(counts = aggregated_assay), colData = new_colData)

    return(new_se)
}

compute_performance <- function(modRes, rowdata) {
    df <- data.frame()
    for (name in names(modRes)) {
        curr <- data.frame(
            sid = rownames(modRes[[name]]),
            pval = modRes[[name]]$pval,
            adjPval = modRes[[name]]$adjPval,
            is_da = rowdata[rownames(modRes[[name]]), "shifted"],
            method = rep(name, nrow(modRes[[name]]))
        )
        df <- rbind(df, curr)
    }
    rownames(df) <- NULL
    truth_df <- df %>%
        select(sid, is_da) %>%
        distinct() %>%
        column_to_rownames("sid")

    adjPval_df <- df %>%
        select(sid, method, adjPval) %>%
        pivot_wider(names_from = method, values_from = adjPval) %>%
        column_to_rownames("sid")

    cobdata <- COBRAData(
        padj = as.data.frame(adjPval_df),
        truth = as.data.frame(truth_df)
    )

    perf <- calculate_performance(cobdata,
                                  binary_truth = "is_da",
                                  aspects = "fdrtpr",
                                  maxsplit = Inf,
                                  thrs = seq(from = 0.0001,  to = 0.2, by = 0.0001)
    ) %>%
        fdrtpr() %>%
        mutate(thr = as.numeric(sub("thr", "", thr)))
    return(perf)
}
