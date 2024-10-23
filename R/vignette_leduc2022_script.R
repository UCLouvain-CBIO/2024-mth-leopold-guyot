library("SCP.replication")

## Core packages of this workflow
library(SingleCellExperiment)
library(scp)
library(scpdata)
library(limma)
## Utility packages for data manipulation and visualization
library(tidyverse)
library(patchwork)

leduc2022script <- function(leduc){
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

    leduc <- filterFeatures(leduc, ~ Potential.contaminant != "+" &
                                !grepl("CON", Proteins) &
                                Reverse != "+" &
                                !grepl("REV", Leading.razor.protein) &
                                (is.na(PIF) | PIF > 0.6) &
                                dart_qval < 0.01)

    table(leduc$SampleType)

    leduc <- computeSCR(leduc, names(leduc),
                        colvar = "SampleType",
                        samplePattern = "Mel|Macro",
                        carrierPattern = "Carrier",
                        sampleFUN = "mean",
                        rowDataName = "MeanSCR")

    rbindRowData(leduc, i = names(leduc)) %>%
        data.frame %>%
        ggplot(aes(x = MeanSCR)) +
        geom_histogram() +
        geom_vline(xintercept = 0.1) +
        scale_x_log10()

    leduc <- filterFeatures(leduc, ~
                                !is.na(MeanSCR) & !is.infinite(MeanSCR) &
                                MeanSCR < 0.05)

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
        data.frame %>%
        ggplot(aes(x = ScSums)) +
        geom_histogram() +
        scale_x_log10()

    leduc <- filterFeatures(leduc, ~ ScSums != 0)


    leduc <- divideByReference(leduc, i = names(leduc),
                               colvar = "SampleType",
                               samplePattern = ".",
                               refPattern = "Reference")

    remove.duplicates <- function(x)
        apply(x, 2, function(xx) xx[which(!is.na(xx))[1]] )


    peptideAssays <- paste0("peptides_", names(leduc))


    leduc <- aggregateFeaturesOverAssays(leduc,
                                         i = names(leduc),
                                         fcol = "modseq",
                                         name = peptideAssays,
                                         fun = remove.duplicates)

    leduc

    ## Generate a list of DataFrames with the information to modify
    rbindRowData(leduc, i = grep("^pep", names(leduc))) %>%
        data.frame %>%
        group_by(modseq) %>%
        ## The majority vote happens here
        mutate(Leading.razor.protein.symbol =
                   names(sort(table(Leading.razor.protein),
                              decreasing = TRUE))[1]) %>%
        select(modseq, Leading.razor.protein.symbol) %>%
        filter(!duplicated(modseq, Leading.razor.protein.symbol)) ->
        ppMap
    consensus <- lapply(peptideAssays, function(i) {
        ind <- match(rowData(leduc[[i]])$modseq, ppMap$modseq)
        DataFrame(Leading.razor.protein.symbol =
                      ppMap$Leading.razor.protein.symbol[ind])
    })
    ## Name the list
    names(consensus) <- peptideAssays
    ## Modify the rowData
    rowData(leduc) <- consensus

    ## Cleaning missing data

    leduc <- infIsNA(leduc, i = peptideAssays)
    leduc <- zeroIsNA(leduc, i = peptideAssays)

    ## Join assays

    leduc <- joinAssays(leduc,
                        i = peptideAssays,
                        name = "peptides")

    leduc

    # Filter single-cells based on median CV


    leduc <- medianCVperCell(leduc,
                             i = peptideAssays,
                             groupBy = "Leading.razor.protein.symbol",
                             nobs = 3,
                             na.rm = TRUE,
                             colDataName = "MedianCV",
                             norm = "SCoPE2")

    colData(leduc) %>%
        data.frame %>%
        filter(grepl("Mono|Mel|Neg", SampleType)) %>%
        mutate(control = ifelse(grepl("Neg", SampleType), "ctl", "sc")) %>%
        ggplot() +
        aes(x = MedianCV,
            fill = control) +
        geom_density(alpha = 0.5, adjust = 1) +
        geom_vline(xintercept = 0.42) +
        xlim(0.2, 0.8) +
        theme_minimal() +
        scale_fill_manual(values = c( "black", "purple2")) +
        xlab("Quantification variability") +
        ylab("Fraction of cells")

    leduc <-
        subsetByColData(leduc,
                        !is.na(leduc$MedianCV) &
                            leduc$MedianCV < 0.42 &
                            grepl("Mono|Mel", leduc$SampleType))


    # Normalization


    ## Scale column with median
    leduc <- normalize(leduc,
                       i = "peptides",
                       method = "div.median",
                       name = "peptides_norm1")
    ## Scale rows with median
    leduc <- sweep(leduc,
                   i = "peptides_norm1",
                   name = "peptides_norm2",
                   MARGIN = 1,
                   FUN = "/",
                   STATS = rowMedians(assay(leduc[["peptides_norm1"]]),
                                      na.rm = TRUE))

    leduc <- filterNA(leduc,
                      i = "peptides_norm2",
                      pNA = 0.99)

    nnaRes <- nNA(leduc, "peptides_norm2")
    sel <- nnaRes$nNAcols$pNA < 99
    leduc[["peptides_norm2"]] <- leduc[["peptides_norm2"]][, sel]

    # Log-transformation

    leduc <- logTransform(leduc,
                          base = 2,
                          i = "peptides_norm2",
                          name = "peptides_log")


    # Aggregate peptide data to protein data

    leduc <- aggregateFeatures(leduc,
                               i = "peptides_log",
                               name = "proteins",
                               fcol = "Leading.razor.protein.symbol",
                               fun = matrixStats::colMedians,
                               na.rm = TRUE)

    # Normalization

    ## Center columns with median
    leduc <- normalize(leduc,
                       i = "proteins",
                       method = "center.median",
                       name = "proteins_norm1")
    ## Scale rows with median
    leduc <- sweep(leduc,
                   i = "proteins_norm1",
                   name = "proteins_norm2",
                   MARGIN = 1,
                   FUN = "-",
                   STATS = rowMedians(assay(leduc[["proteins_norm1"]]),
                                      na.rm = TRUE))

    # Imputation

    data.frame(pNA = nNA(leduc, "proteins_norm2")$nNAcols$pNA) %>%
        ggplot(aes(x = pNA)) +
        geom_histogram() +
        xlab("Percentage missingnes per cell")


    leduc <- imputeKnnSCoPE2(leduc,
                             i = "proteins_norm2",
                             name = "proteins_impd",
                             k = 3)

    leduc <- impute(leduc,
                    i = "proteins_norm2",
                    method = "knn",
                    k = 3, rowmax = 1, colmax= 1,
                    maxp = Inf, rng.seed = 1234)

    # Batch correction

    sce <- getWithColData(leduc, "proteins_impd")

    model <- model.matrix(~ SampleType, data = colData(sce))
    assay(sce) <- removeBatchEffect(x = assay(sce),
                                    batch = sce$lcbatch,
                                    batch2 = sce$Channel,
                                    design = model)

    leduc <- addAssay(leduc, y = sce, name = "proteins_batchC")
    leduc <- addAssayLinkOneToOne(leduc, from = "proteins_impd",
                                  to = "proteins_batchC")

    # Normalization

    ## Center columns with median
    leduc <- normalize(leduc,
                       i = "proteins_batchC",
                       method = "center.median",
                       name = "proteins_batchC_norm1")
    ## Scale rows with median
    leduc <- sweep(leduc,
                   i = "proteins_batchC_norm1",
                   name = "proteins_processed",
                   MARGIN = 1,
                   FUN = "-",
                   STATS = rowMedians(assay(leduc[["proteins_batchC_norm1"]]),
                                      na.rm = TRUE))

    # PCA

    sce <- getWithColData(leduc, "proteins_processed")
    pcaRes <- pcaSCoPE2(sce)
    ## Compute percent explained variance
    pcaPercentVar <- round(pcaRes$values[1:2] / sum(pcaRes$values) * 100)
    ## Plot PCA
    data.frame(PC = pcaRes$vectors[, 1:2],
               colData(sce)) %>%
        ggplot() +
        aes(x = PC.1,
            y = PC.2,
            colour = SampleType) +
        geom_point(alpha = 0.5) +
        xlab(paste0("PC1 (", pcaPercentVar[1], "%)")) +
        ylab(paste0("PC2 (", pcaPercentVar[2], "%)"))+
        ggtitle("PCA on scp processed protein data")

    library(scater)
    ## Perform PCA, see ?runPCA for more info about arguments
    runPCA(sce, ncomponents = 50,
           ntop = Inf,
           scale = TRUE,
           exprs_values = 1,
           name = "PCA") %>%
        ## Plotting is performed in a single line of code
        plotPCA(colour_by = "SampleType")
}



sd <- benchmarkme::get_sys_details()
cat("Machine: ", sd$sys_info$sysname, " (", sd$sys_info$release, ")\n",
    "R version: R.", sd$r_version$major, ".", sd$r_version$minor,
    " (svn: ", sd$r_version$`svn rev`, ")\n",
    "RAM: ", round(sd$ram / 1E9, 1), " GB\n",
    "CPU: ", sd$cpu$no_of_cores, " core(s) - ", sd$cpu$model_name, "\n",
    sep = "")

### Timing

timing <- Sys.time() - timeStart
cat(timing[[1]], attr(timing, "units"))

### Memory

format(object.size(leduc), units = "GB")
