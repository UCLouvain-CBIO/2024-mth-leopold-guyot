############### Packages ################

library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggpmisc)

############### Vignette Benchmark ############

inputDir <- "dataOutput/vignetteBenchmark"
objectSize <- read.table(
    file = file.path(inputDir, "size_report.tsv"),
    header = TRUE
)
objectSize$nCell <- as.factor(objectSize$nCell)
objectSize$rep <- as.factor(objectSize$rep)

dataFrames <- lapply(
    list.files(
        path = file.path(
            inputDir,
            "memoryOutput"
        ),
        full.names = TRUE
    ),
    function(x) {
        df <- read.csv(x)
        total_time <- sum(df$Elapsed_Time_sec, na.rm = TRUE)
        df %>%
            mutate(ratioTime = Elapsed_Time_sec / total_time)
    }
)
names(dataFrames) <- list.files(path = file.path(inputDir, "memoryOutput"), full.names = FALSE)
names(dataFrames) <- gsub(pattern = ".csv", replacement = "", x = names(dataFrames))

benchmarkDF <- bind_rows(dataFrames, .id = "run")
benchmarkDF <- separate_wider_delim(
    data = benchmarkDF, cols = c("run"), delim = "_",
    names = c("nCell", "replicate")
)

benchmarkDF$vignette_step <- gsub(
    pattern = "(.*leduc)(.*)([.,/(]leduc.*)",
    replacement = "\\2",
    x = benchmarkDF$Function_Call
)
irrelevantIndices <- c(
    grep("write.table", benchmarkDF$vignette_step),
    grep("assaysNames", benchmarkDF$vignette_step)
)
benchmarkDF <- benchmarkDF[-irrelevantIndices, ]


plotSize <- objectSize %>%
    group_by(nCell) %>%
    summarise(sizeGB = mean(size) / (1024^3)) %>% # Divide by 1024^3 to get GB
    ggplot(aes(x = nCell, y = sizeGB, fill = nCell)) +
    geom_col() + xlab(label = "Number of cells") + ylab(label = "Final object size in GB") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
    ) +
    guides(color = guide_legend(title = NULL))

plotTime <- benchmarkDF %>%
    # Convert `nCell` to a factor with levels ordered numerically
    mutate(nCell = factor(nCell, levels = sort(as.numeric(unique(nCell))))) %>%
    group_by(nCell, replicate) %>%
    summarise(time = sum(Elapsed_Time_sec) / 60) %>%
    ggplot(aes(x = nCell, y = time, color = nCell)) +
    geom_boxplot() + xlab(label = "Number of cells") + ylab(label = "Global pipeline runtime (minutes)") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )

plotPeak <- benchmarkDF %>%
    group_by(nCell, replicate) %>%
    summarise(peak = max(Peak_RAM_Used_MiB) / 1024, .groups = "drop") %>% # Convert to GB
    group_by(nCell) %>%
    summarise(peak = mean(peak), .groups = "drop") %>%
    mutate(nCell = factor(nCell, levels = sort(as.numeric(unique(nCell))))) %>%
    ggplot(aes(x = nCell, y = peak, fill = nCell)) +
    geom_col() +
    labs(y = "Peak RAM Used (GB)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

plotPeakBoxPlot <- benchmarkDF %>%
    group_by(nCell, replicate) %>%
    summarise(peak = max(Peak_RAM_Used_MiB) / 1024, .groups = "drop") %>%  # Convert to GB
    mutate(nCell = factor(nCell, levels = sort(as.numeric(unique(nCell))))) %>%
    ggplot(aes(x = nCell, y = peak, color = nCell)) +
    geom_boxplot() +
    labs(y = "Peak RAM Used (GB)", x = "Number of cells") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )

ggsave("Figs/report/plotTimeVignette.pdf", plotTime, width = 10, height = 5)
ggsave("Figs/report/plotSizeVignette.pdf", plotSize, width = 10, height = 6)
ggsave("Figs/report/plotPeakVignette.pdf", plotPeakBoxPlot, width = 10, height = 5)


plotRatio <- benchmarkDF %>%
    group_by(vignette_step) %>%
    summarise(meanTimePercentage = mean(ratioTime) * 100) %>%
    arrange(desc(meanTimePercentage)) %>%
    mutate(vignette_step = factor(vignette_step, levels = vignette_step)) %>%
    ggplot(aes(x = vignette_step, y = meanTimePercentage)) +
    geom_bar(stat = "identity") +
    labs(
        title = "Proportion of vignette_step by meanTimePercentage",
        y = "Percentage of time taken (average)", x = "Vignette Step"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


plotRatio1000 <- benchmarkDF %>%
    filter(nCell == "1000") %>%
    group_by(vignette_step) %>%
    summarise(meanTimePercentage = mean(ratioTime) * 100) %>%
    arrange(desc(meanTimePercentage)) %>%
    mutate(vignette_step = factor(vignette_step, levels = vignette_step)) %>%
    ggplot(aes(x = vignette_step, y = meanTimePercentage)) +
    geom_bar(stat = "identity") +
    labs(
        title = "Proportion of vignette_step by meanTimePercentage",
        y = "Percentage of time taken (average)", x = "Vignette Step"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


plotRatio16000 <- benchmarkDF %>%
    filter(nCell == "16000") %>%
    group_by(vignette_step) %>%
    summarise(meanTimePercentage = mean(ratioTime) * 100) %>%
    arrange(desc(meanTimePercentage)) %>%
    mutate(vignette_step = factor(vignette_step, levels = vignette_step)) %>%
    ggplot(aes(x = vignette_step, y = meanTimePercentage)) +
    geom_bar(stat = "identity") +
    labs(
        title = "Proportion of vignette_step by meanTimePercentage",
        y = "Percentage of time taken (average)", x = "Vignette Step"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

step_order <- benchmarkDF %>%
    filter(nCell == 1000) %>%
    group_by(vignette_step) %>%
    summarise(meanTimePercentage1000 = mean(ratioTime) * 100) %>%
    arrange(desc(meanTimePercentage1000)) %>%  # Order by mean time
    pull(vignette_step)  # Extract ordered vignette_step values

# Process the full dataset
plotRatioComplete <- benchmarkDF %>%
    group_by(vignette_step, nCell) %>%
    summarise(meanTimePercentage = mean(ratioTime) * 100) %>%
    ungroup() %>%
    mutate(nCell = as.numeric(as.character(nCell))) %>%
    mutate(nCell = factor(nCell, levels = sort(unique(nCell)))) %>%  # Order numerically
    mutate(vignette_step = factor(vignette_step, levels = step_order)) %>%
    ggplot(aes(x = vignette_step, y = meanTimePercentage, fill = nCell)) +
    geom_bar(position = "dodge", stat = "identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Your desired mapping
rename_map <- c(
    "AggPSM" = "aggregatePSM",
    "SCR" = "filterSCR",
    "ScSums" = "filterSumSignal",
    "NormToRef" = "divideByReference",
    "Consensus" = "consensusMap",
    "MissingData" = "zeroToNA",
    "Join" = "joinPeptide",
    "FilterCV" = "filterCV",
    "NormPep" = "normPeptide",
    "FilteringNA" = "filterNA",
    "LogTransfo" = "logTransformation",
    "AggPep" = "aggregatePeptide",
    "NormPro" = "normProt1",
    "Impute" = "impute",
    "Batch" = "removeBatch",
    "NormBatch" = "normProt2",
    "FilterFeatures" = "filterPSM"
)

# Your pipeline order
step_order <- c(
    "filterPSM",
    "filterSCR",
    "filterSumSignal",
    "divideByReference",
    "aggregatePSM",
    "consensusMap",
    "zeroToNA",
    "joinPeptide",
    "filterCV",
    "normPeptide",
    "filterNA",
    "logTransformation",
    "aggregatePeptide",
    "normProt1",
    "impute",
    "removeBatch",
    "normProt2"
)

valid_original_levels <- intersect(names(rename_map), unique(benchmarkDF$vignette_step))
rename_map <- rename_map[valid_original_levels]

plotRatioComplete <- benchmarkDF %>%
    mutate(vignette_step = recode(vignette_step, !!!rename_map)) %>%
    filter(vignette_step %in% step_order) %>%   # Remove unneeded steps
    group_by(vignette_step, nCell) %>%
    summarise(meanTimePerStep = mean(Elapsed_Time_sec), .groups = "drop") %>%
    mutate(nCell = as.numeric(as.character(nCell))) %>%
    mutate(nCell = factor(nCell, levels = sort(unique(nCell)))) %>%
    mutate(vignette_step = factor(vignette_step, levels = step_order)) %>%
    ggplot(aes(x = vignette_step, y = meanTimePerStep, fill = nCell)) +
    geom_bar(position = "dodge", stat = "identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Pipeline steps") + ylab("Mean time per step (seconds)")


plotRatioComplete
ggsave("Figs/report/vignetteBenchmark_ratioComplete.pdf", width = 10, height = 6)
ggsave("Figs/vignetteBenchmark_ratioComplete.png", width = 10, height = 5)

plotTime / plotPeak / plotSize
ggsave("Figs/vignetteBenchmark_summary.pdf", width = 7, height = 7)
ggsave("Figs/vignetteBenchmark_summary.png", width = 7, height = 7)

############ 4 var Benchmark ##############

convert_bytes_to_mib <- function(bytes) {
    mib <- bytes / (2^20)
    return(mib)
}

sizeTable <- read.table("dataOutput/4_variables_benchmark/qfeatures_size_report.tsv", header = TRUE)

sizeCols <- grep("size", colnames(sizeTable))
for (sizeCol in sizeCols) {
    sizeTable[[paste0(colnames(sizeTable)[[sizeCol]], "_MiB")]] <- convert_bytes_to_mib(sizeTable[[sizeCol]])
}
sizeTable$totalCells <- sizeTable$nCell * sizeTable$nAssay


sizeLFTMT <- sizeTable %>%
    filter(nFeat == 4000,
           nCol == 100,
           totalCells == 256) %>%
    mutate(totalCells = as.factor(totalCells),
           nCell = as.factor(nCell)) %>%
    ggplot(aes(x = nCell, y = sizeTotal_MiB)) + geom_bar(stat = "identity")

sizeLFTMT

propTable <- sizeTable %>%
    filter(nCol == 50,
           nFeat == 4000) %>%
    filter(totalCells == 256) %>%
    mutate(nCell = as.factor(nCell)) %>%
    pivot_longer(cols = 10:12, names_to = "component", values_to = "size")

sizeWithProp <- ggplot(propTable, aes(fill = component, y = size, x = nCell)) +
    geom_bar(position="stack", stat="identity") +
    ylab("Total Size (MiB)") +
    xlab("Number of Cells by Assay")
ggsave(filename = "sizeByAssaySize.png")
library(rlang)

create_plot <- function(df, x_var, color_var, fixed_filters) {
    df_filtered <- df %>%
        filter(!!!fixed_filters)

    ggplot(df_filtered, aes(x = .data[[x_var]], y = sizeTotal_MiB, color = factor(.data[[color_var]]))) +
        geom_point() +
        geom_line(alpha = 0.8) +
        labs(color = color_var, x = x_var, y = "Total size (MiB)")
}

all_vars <- c("nCell", "nFeat", "nCol", "nAssay")

fixed_values <- list(nCol = 50, nAssay = 128, nFeat = 2000, nCell = 16)

plots <- list()

for (x_var in all_vars) {
    remaining_vars <- setdiff(all_vars, x_var)
    for (color_var in remaining_vars) {

        fixed_vars <- setdiff(remaining_vars, color_var)
        fixed_filters <- list(quo(!!sym(fixed_vars[1]) == fixed_values[[fixed_vars[1]]]),
                              quo(!!sym(fixed_vars[2]) == fixed_values[[fixed_vars[2]]]))

        plots[[paste(x_var, color_var, sep = "_")]] <- create_plot(sizeTable, x_var, color_var, fixed_filters)
    }
}
theme_update(plot.title = element_text(hjust = 0.5))
combined_plot <-
    ((plots[["nCell_nFeat"]] + ggtitle("Size by nCell") + theme(axis.title.x=element_blank(),
                                                             axis.text.x=element_blank(),
                                                             axis.ticks.x=element_blank(),
                                                             plot.title = element_text(size=22))) /
         (plots[["nCell_nCol"]] + theme(axis.title.x=element_blank(),
                                       axis.text.x=element_blank(),
                                       axis.ticks.x=element_blank())) /
         plots[["nCell_nAssay"]]) |
    ((plots[["nFeat_nCol"]] + ggtitle("Size by nFeat") + theme(axis.title.x=element_blank(),
                                                                axis.text.x=element_blank(),
                                                                axis.ticks.x=element_blank(),
                                                                plot.title = element_text(size=22))) /
         (plots[["nFeat_nAssay"]] + theme(axis.title.x=element_blank(),
                                         axis.text.x=element_blank(),
                                         axis.ticks.x=element_blank())) /
         plots[["nFeat_nCell"]]) |
    ((plots[["nCol_nFeat"]] + ggtitle("Size by nCol") + theme(axis.title.x=element_blank(),
                                                               axis.text.x=element_blank(),
                                                               axis.ticks.x=element_blank(),
                                                               plot.title = element_text(size=22))) /
         (plots[["nCol_nAssay"]] + theme(axis.title.x=element_blank(),
                                        axis.text.x=element_blank(),
                                        axis.ticks.x=element_blank())) /
         plots[["nCol_nCell"]]) |
    ((plots[["nAssay_nCol"]] + ggtitle("Size by nAssay") + theme(axis.title.x=element_blank(),
                                                                  axis.text.x=element_blank(),
                                                                  axis.ticks.x=element_blank(),
                                                                  plot.title = element_text(size=22))) /
         (plots[["nAssay_nFeat"]] + theme(axis.title.x=element_blank(),
                                         axis.text.x=element_blank(),
                                         axis.ticks.x=element_blank())) /
         plots[["nAssay_nCell"]])

combined_plot
ggsave("Figs/4var_global.pdf", width = 18, height = 10)
ggsave("Figs/4var_global.png", width = 18, height = 10)

#### STEPS BENCHMARK

inputDir <- "dataOutput/individualStepsBenchmark"
files <- list.files(path = inputDir, full.names = TRUE)
fileSizeReport <- "dataOutput/individualStepsBenchmark/qfeatures_size_report.tsv"
stepResults <- files[!files == fileSizeReport]

stepDataList <- lapply(stepResults, function(x) {
    l <- as.list(read.csv(x)[1, 3:5])
    l$name <- sub(".csv", "", basename(x))
    l
})

sizeTable <- read.table(fileSizeReport, header = TRUE)

peakRamTable <- do.call(rbind, lapply(stepDataList, as.data.frame))
peakRamTable <- separate(peakRamTable, name, into = c("nCell", "step", "replicate"), sep = "_")

nCell_order <- c("500", "1000", "2000", "4000")
step_order <- c(
    "benchmarkFilterFeatures",
    "benchmarkFilterSamples",
    "benchmarkZeroisNA",
    "benchmarkAggPSM",
    "benchmarkJoinPSM",
    "benchmarkNormSampPep",
    "benchmarkNormFeatPep",
    "benchmarkLogPep",
    "benchmarkAggPep",
    "benchmarkImputePro"
)

timeMedian <- peakRamTable %>%
    group_by(nCell, step) %>%
    summarise(medianTimeElasped = median(Elapsed_Time_sec)) %>%
    ungroup() %>%
    mutate(
        nCell = factor(nCell, levels = nCell_order),
        step = factor(step, levels = rev(step_order))
    )


timeMed <- ggplot(timeMedian, aes(nCell, step, fill = medianTimeElasped)) +
             geom_tile() +
             scale_fill_gradient2(
                 low = "white", high = "red",
                 midpoint = median(peakRamTable$Elapsed_Time_sec, na.rm = TRUE),
                 transform = "log2"
             ) + xlab("Number of cells") + ylab("Function names")
ggsave("Figs/report/singleStepTime.pdf", timeMed, width = 9, height = 7)
convert_mib_to_gb <- function(mib) {
    gb <- mib * 0.001024
    return(gb)
}
peakRamMedian <- peakRamTable %>%
    group_by(nCell, step) %>%
    summarise(medianPeakRam = median(Peak_RAM_Used_MiB)) %>%
    ungroup() %>%
    mutate(
        nCell = factor(nCell, levels = nCell_order),
        step = factor(step, levels = rev(step_order))
    )

ggplot(peakRamMedian, aes(nCell, step, fill = medianPeakRam)) +
             geom_tile() +
             scale_fill_gradient2(
                 low = "white", high = "red",
                 midpoint = median(peakRamTable$Peak_RAM_Used_MiB, na.rm = TRUE),
                 transform = "log2"
             )
ggsave("Figs/singleStepRam.png")
convert_bytes_to_mib <- function(bytes) {
    mib <- bytes / (2^20)
    return(mib)
}

sizeComponentDiff <- sizeTable %>%
    filter(nCell == 4000) %>%
    group_by(step) %>%
    summarise(
        medianDiffAssay = convert_bytes_to_mib(median(sizeAssayAfter - sizeAssayBefore)),
        medianDiffRowData = convert_bytes_to_mib(median(sizeRowDataAfter - sizeRowDataBefore)),
        medianDiffColData = convert_bytes_to_mib(median(sizeColDataAfter - sizeColDataBefore))
    ) %>%
    ungroup() %>%
    pivot_longer(
        cols = c(medianDiffAssay, medianDiffRowData, medianDiffColData),
        names_to = "component", values_to = "medianSizeMiB"
    ) %>%
    mutate(step = factor(step, levels = rev(step_order)))

peakRamTableOrdered <- peakRamTable %>%
    mutate(
        step = factor(step, levels = rev(step_order)),
        nCell = as.integer(nCell)
    ) %>%
    group_by(nCell, step) %>%
    summarise("medianTimeSec" = median(Elapsed_Time_sec)) %>%
    filter(step == "benchmarkFilterSamples")

ggplotly(ggplot(peakRamTableOrdered, aes(x = nCell, y = medianTimeSec)) +
             geom_line())
ggsave("Figs/singleStepSample.png")

filterSample <- peakRamTable %>%
    filter(step == "benchmarkFilterSamples") %>%
    mutate(
        nCellFactor = factor(nCell, levels = nCell_order),
        nCellNum = as.numeric(nCell)
    )

aggPSM <- peakRamTable %>%
    filter(step == "benchmarkAggPSM") %>%
    mutate(
        nCellFactor = factor(nCell, levels = nCell_order),
        nCellNum = as.numeric(nCell)
    )

joinPSM <- peakRamTable %>%
    filter(step == "benchmarkJoinPSM") %>%
    mutate(
        nCellFactor = factor(nCell, levels = nCell_order),
        nCellNum = as.numeric(nCell)
    )

subsetCombined <- peakRamTable %>%
    filter(step %in% c("benchmarkJoinPSM", "benchmarkAggPSM", "benchmarkFilterSamples")) %>%
    group_by(step, nCell) %>%
    mutate(
        nCellFactor = factor(nCell, levels = nCell_order),
        nCellNum = as.numeric(nCell)
    )

# Split data
quad_data <- subsetCombined %>% filter(step == "benchmarkFilterSamples")
lin_data  <- subsetCombined %>% filter(step != "benchmarkFilterSamples")

p_time_mixed <- ggplot(subsetCombined, aes(x = nCellNum, y = Elapsed_Time_sec)) +
    geom_boxplot(
        aes(group = interaction(step, nCellNum), color = step),
        alpha = 0.4,
        position = position_dodge(width = 0.75)
    ) +

    # Quadratic smooth + annotation
    geom_smooth(
        data = quad_data,
        aes(color = step),
        method = "lm",
        formula = y ~ poly(x, 2),
        se = TRUE
    ) +
    stat_poly_eq(
        data = quad_data,
        aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
        formula = y ~ poly(x, 2),
        parse = TRUE,
        label.x.npc = "left",
        label.y.npc = "top",
        size = 3
    ) +

    # Linear smooth + annotation
    geom_smooth(
        data = lin_data,
        aes(color = step),
        method = "lm",
        formula = y ~ x,
        se = TRUE
    ) +
    stat_poly_eq(
        data = lin_data,
        aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
        formula = y ~ x,
        parse = TRUE,
        label.x.npc = "left",
        label.y.npc = "top",
        size = 3
    ) +

    facet_wrap(~ step, scales = "fixed") +
    ylim(0, NA) +
    xlab("Number of Cells") +
    ylab("Elapsed Time (seconds)") +
    theme_minimal() +
    scale_x_continuous(breaks = as.numeric(nCell_order)) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )

ggsave("Figs/report/substepTime.pdf", p_time_mixed, width = 8, height = 4)


p_peak <- ggplot(subsetCombined, aes(x = nCellNum, y = Peak_RAM_Used_MiB)) +
    geom_boxplot(
        aes(group = interaction(step, nCellNum), color = step),
        alpha = 0.4,
        position = position_dodge(width = 0.75)
    ) +
    # Linear smooth + annotation
    geom_smooth(
        aes(color = step),
        method = "lm",
        formula = y ~ x,
        se = TRUE
    ) +
    stat_poly_eq(
        aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
        formula = y ~ x,
        parse = TRUE,
        label.x.npc = "left",
        label.y.npc = "top",
        size = 3
    ) +

    facet_wrap(~ step, scales = "fixed") +
    ylim(0, NA) +
    xlab("Number of Cells") +
    ylab("Peak RAM (MiB)") +
    scale_x_continuous(breaks = as.numeric(nCell_order)) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
    )

ggsave("Figs/report/substepPeak.pdf", p_peak, width = 8, height = 4)

## agg and join diff size

sizeComponentDiff %>%
    filter(step %in% c("benchmarkAggPSM", "benchmarkJoinPSM")) %>%
    ggplot(aes(x = component, y = step, fill = medianSizeMiB)) +
    geom_tile() +
    scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = 0)
