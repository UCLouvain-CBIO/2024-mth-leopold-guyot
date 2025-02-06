library("SCP.replication")

## Core packages of this workflow
library(SingleCellExperiment)
library(scp)
library(scpdata)
library(limma)
# library(scater)

## Utility packages for data manipulation and visualization
library(tidyverse)
library(patchwork)

leduc2022script <- function(leduc) {
    assaysNames <- names(leduc)
    # Remove contaminant, noisy and low-confidence spectra
    leduc <- leducFilterFeatures(leduc)
    # Sample to carrier filter
    leduc <- leducSCR(leduc)
    # Filter on summed single-cell signal
    leduc <- leducScSums(leduc)
    # Normalize to reference
    leduc <- leducNormToRef(leduc)
    # Aggregate PSM data to peptide data
    leduc <- leducAggPSM(leduc)
    # Consensus mapping of peptides to proteins
    leduc <- leducConsensus(leduc, assaysNames)
    # Cleaning missing data
    leduc <- leducMissingData(leduc, assaysNames)
    # Join assays
    leduc <- leducJoin(leduc, assaysNames)
    # Filter single-cells based on median CV
    leduc <- leducFilterCV(leduc, assaysNames)
    # Normalization peptides
    leduc <- leducNormPep(leduc)
    # Missing data filtering
    leduc <- leducFilteringNA(leduc)
    # Log-transformation
    leduc <- leducLogTransfo(leduc)
    # Aggregate peptide data to protein data
    leduc <- leducAggPep(leduc)
    # Normalization proteins
    leduc <- leducNormPro(leduc)
    # Imputation
    leduc <- leducImpute(leduc)
    # Batch correction
    leduc <- leducBatch(leduc)
    # Normalization batch corrected proteins
    leduc <- leducNormBatch(leduc)
    # PCA
    # leduc <- leducPCA(leduc)

    return(object.size(leduc))
}

leducFilterFeatures <- function(leduc) {
    rowDataNames(leduc)


    rd <- data.frame(rbindRowData(leduc, i = names(leduc)))
    ggplot(rd) +
        aes(x = dart_PEP) +
        geom_histogram() +
        geom_vline(xintercept = 0.01) +
        ggplot(rd) +
        aes(x = PIF) +
        geom_histogram() +
        geom_vline(xintercept = 0.6)

    filterFeatures(leduc, ~ Potential.contaminant != "+" &
        !grepl("CON", Proteins) &
        Reverse != "+" &
        !grepl("REV", Leading.razor.protein) &
        (is.na(PIF) | PIF > 0.6) &
        dart_qval < 0.01)
}

leducSCR <- function(leduc) {
    table(leduc$SampleType)

    leduc <- computeSCR(leduc, names(leduc),
        colvar = "SampleType",
        samplePattern = "Mel|Macro",
        carrierPattern = "Carrier",
        sampleFUN = "mean",
        rowDataName = "MeanSCR"
    )

    rbindRowData(leduc, i = names(leduc)) %>%
        data.frame() %>%
        ggplot(aes(x = MeanSCR)) +
        geom_histogram() +
        geom_vline(xintercept = 0.1) +
        scale_x_log10()

    filterFeatures(leduc, ~
        !is.na(MeanSCR) & !is.infinite(MeanSCR) &
            MeanSCR < 0.05)
}

leducScSums <- function(leduc) {
    sums <- lapply(names(leduc), function(i) {
        sce <- leduc[[i]]
        sel <- grep("Mel|Macro|Neg", colData(leduc)[colnames(sce), "SampleType"])
        x <- assay(sce)[, sel, drop = FALSE]
        rs <- rowSums(x, na.rm = TRUE)
        DataFrame(ScSums = rs)
    })

    names(sums) <- names(leduc)
    rowData(leduc) <- sums

    rbindRowData(leduc, i = names(leduc)) %>%
        data.frame() %>%
        ggplot(aes(x = ScSums)) +
        geom_histogram() +
        scale_x_log10()

    filterFeatures(leduc, ~ ScSums != 0)
}

leducNormToRef <- function(leduc) {
    divideByReference(leduc,
        i = names(leduc),
        colvar = "SampleType",
        samplePattern = ".",
        refPattern = "Reference"
    )
}

leducAggPSM <- function(leduc) {
    remove.duplicates <- function(x) {
        apply(x, 2, function(xx) xx[which(!is.na(xx))[1]])
    }


    peptideAssays <- paste0("peptides_", names(leduc))


    leduc <- aggregateFeaturesOverAssays(leduc,
        i = names(leduc),
        fcol = "modseq",
        name = peptideAssays,
        fun = remove.duplicates
    )

    leduc
}

leducConsensus <- function(leduc, assaysNames) {
    peptideAssays <- paste0("peptides_", assaysNames)

    ## Generate a list of DataFrames with the information to modify
    rbindRowData(leduc, i = grep("^pep", names(leduc))) %>%
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
        filter(!duplicated(modseq, Leading.razor.protein.symbol)) ->
    ppMap
    consensus <- lapply(peptideAssays, function(i) {
        ind <- match(rowData(leduc[[i]])$modseq, ppMap$modseq)
        DataFrame(
            Leading.razor.protein.symbol =
                ppMap$Leading.razor.protein.symbol[ind]
        )
    })
    ## Name the list
    names(consensus) <- peptideAssays
    ## Modify the rowData
    rowData(leduc) <- consensus

    leduc
}

leducMissingData <- function(leduc, assaysNames) {
    peptideAssays <- paste0("peptides_", assaysNames)

    leduc <- infIsNA(leduc, i = peptideAssays)
    zeroIsNA(leduc, i = peptideAssays)
}

leducJoin <- function(leduc, assaysNames) {
    peptideAssays <- paste0("peptides_", assaysNames)

    leduc <- joinAssays(leduc,
        i = peptideAssays,
        name = "peptides"
    )

    leduc
}

leducFilterCV <- function(leduc, assaysNames) {
    peptideAssays <- paste0("peptides_", assaysNames)

    leduc <- medianCVperCell(leduc,
        i = peptideAssays,
        groupBy = "Leading.razor.protein.symbol",
        nobs = 3,
        na.rm = TRUE,
        colDataName = "MedianCV",
        norm = "SCoPE2"
    )

    colData(leduc) %>%
        data.frame() %>%
        filter(grepl("Mono|Mel|Neg", SampleType)) %>%
        mutate(control = ifelse(grepl("Neg", SampleType), "ctl", "sc")) %>%
        ggplot() +
        aes(
            x = MedianCV,
            fill = control
        ) +
        geom_density(alpha = 0.5, adjust = 1) +
        geom_vline(xintercept = 0.42) +
        xlim(0.2, 0.8) +
        theme_minimal() +
        scale_fill_manual(values = c("black", "purple2")) +
        xlab("Quantification variability") +
        ylab("Fraction of cells")


    subsetByColData(
        leduc,
        !is.na(leduc$MedianCV) &
            leduc$MedianCV < 0.42 &
            grepl("Mono|Mel", leduc$SampleType)
    )
}

leducNormPep <- function(leduc) {
    ## Scale column with median
    leduc <- QFeatures::normalize(leduc,
        i = "peptides",
        method = "div.median",
        name = "peptides_norm1"
    )
    ## Scale rows with median
    sweep(leduc,
        i = "peptides_norm1",
        name = "peptides_norm2",
        MARGIN = 1,
        FUN = "/",
        STATS = rowMedians(assay(leduc[["peptides_norm1"]]),
            na.rm = TRUE
        )
    )
}

leducFilteringNA <- function(leduc) {
    leduc <- filterNA(leduc,
        i = "peptides_norm2",
        pNA = 0.99
    )

    nnaRes <- nNA(leduc, "peptides_norm2")
    sel <- nnaRes$nNAcols$pNA < 99
    leduc[["peptides_norm2"]] <- leduc[["peptides_norm2"]][, sel]
    leduc
}

leducLogTransfo <- function(leduc) {
    logTransform(leduc,
        base = 2,
        i = "peptides_norm2",
        name = "peptides_log"
    )
}

leducAggPep <- function(leduc) {
    aggregateFeatures(leduc,
        i = "peptides_log",
        name = "proteins",
        fcol = "Leading.razor.protein.symbol",
        fun = matrixStats::colMedians,
        na.rm = TRUE
    )
}

leducNormPro <- function(leduc) {
    ## Center columns with median
    leduc <- QFeatures::normalize(leduc,
        i = "proteins",
        method = "center.median",
        name = "proteins_norm1"
    )
    ## Scale rows with median
    sweep(leduc,
        i = "proteins_norm1",
        name = "proteins_norm2",
        MARGIN = 1,
        FUN = "-",
        STATS = rowMedians(assay(leduc[["proteins_norm1"]]),
            na.rm = TRUE
        )
    )
}

leducImpute <- function(leduc) {
    data.frame(pNA = nNA(leduc, "proteins_norm2")$nNAcols$pNA) %>%
        ggplot(aes(x = pNA)) +
        geom_histogram() +
        xlab("Percentage missingnes per cell")


    leduc <- imputeKnnSCoPE2(leduc,
        i = "proteins_norm2",
        name = "proteins_impd",
        k = 3
    )

    impute(leduc,
        i = "proteins_norm2",
        method = "knn",
        k = 3, rowmax = 1, colmax = 1,
        maxp = Inf, rng.seed = 1234
    )
}

leducBatch <- function(leduc) {
    sce <- getWithColData(leduc, "proteins_impd")

    model <- model.matrix(~SampleType, data = colData(sce))
    assay(sce) <- removeBatchEffect(
        x = assay(sce),
        batch = sce$lcbatch,
        batch2 = sce$Channel,
        design = model
    )

    leduc <- addAssay(leduc, y = sce, name = "proteins_batchC")
    leduc <- addAssayLinkOneToOne(leduc,
        from = "proteins_impd",
        to = "proteins_batchC"
    )
    leduc
}

leducNormBatch <- function(leduc) {
    ## Center columns with median
    leduc <- QFeatures::normalize(leduc,
        i = "proteins_batchC",
        method = "center.median",
        name = "proteins_batchC_norm1"
    )
    ## Scale rows with median
    sweep(leduc,
        i = "proteins_batchC_norm1",
        name = "proteins_processed",
        MARGIN = 1,
        FUN = "-",
        STATS = rowMedians(assay(leduc[["proteins_batchC_norm1"]]),
            na.rm = TRUE
        )
    )
}

leducPCA <- function(leduc) {
    runPCA(sce,
        ncomponents = 50,
        ntop = Inf,
        scale = TRUE,
        exprs_values = 1,
        name = "PCA"
    ) %>%
        ## Plotting is performed in a single line of code
        plotPCA(colour_by = "SampleType")
}
