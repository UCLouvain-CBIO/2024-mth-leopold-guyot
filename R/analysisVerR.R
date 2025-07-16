library("tidyverse")

noSubSCE <- read.csv(file = "dataOutput/optimisationsBench/noSubsetRowDataSCE.csv")
subSCE <- read.csv(file = "dataOutput/optimisationsBench/subsetRowDataSCE.csv")
noSubSE <- read.csv(file = "dataOutput/optimisationsBench/noSubsetRowDataSE.csv")
subSE <- read.csv(file = "dataOutput/optimisationsBench/subsetRowDataSE.csv")

noSubSE <- noSubSE %>%
    separate(X, into = c("nCell", "rep"), sep = "\\.") %>%
    mutate(subsetRowData = FALSE) %>%
    mutate(SE = TRUE) %>%
    filter(version == "base")

subSE <- subSE %>%
    separate(X, into = c("nCell", "rep"), sep = "\\.") %>%
    mutate(subsetRowData = TRUE) %>%
    mutate(SE = TRUE) %>%
    filter(version == "allOpti")

noSubSCE <- noSubSCE %>%
    separate(X, into = c("nCell", "rep"), sep = "\\.") %>%
    mutate(subsetRowData = FALSE) %>%
    mutate(SE = FALSE) %>%
    filter(version != "allOpti")

subSCE <- subSCE %>%
    separate(X, into = c("nCell", "rep"), sep = "\\.") %>%
    mutate(subsetRowData = TRUE) %>%
    mutate(SE = FALSE) %>%
    filter(version == "base")

combined <- rbind(noSubSCE, subSCE, noSubSE, subSE) %>%
    mutate(versionComplete = paste(version, SE, subsetRowData, sep = "_"))

combined_labeled <- combined %>%
    separate(versionComplete, into = c("function_type", "SE", "subsetRowData"), sep = "_") %>%
    mutate(
        SE = SE == "TRUE",
        subsetRowData = subsetRowData == "TRUE",
        versionLabel = case_when(
            function_type == "aggregation" ~ "Aggregation (no subset, no SE)",
            function_type == "allOpti" & subsetRowData & SE ~ "All Optimisations (subset, SE)",
            function_type == "base" & !subsetRowData & !SE ~ "Base (no subset, no SE)",
            function_type == "base" & subsetRowData ~ "Base (subset)",
            function_type == "base" & SE ~ "Base (SE)",
            function_type == "subsetByColData" ~ "Subset by ColData (no subset, no SE)",
            TRUE ~ function_type
        )
    ) %>%
    mutate(
        versionLabel = factor(
            versionLabel,
            levels = c(
                "Base (no subset, no SE)",
                "Subset by ColData (no subset, no SE)",
                "Aggregation (no subset, no SE)",
                "Base (SE)",
                "Base (subset)",
                "All Optimisations (subset, SE)"
            )
        )
    )

## function of time

total_df <- combined_labeled %>%
    group_by(versionLabel, nCell, replicate) %>%
    summarise(totalRuntime = sum(Elapsed_Time_sec), maxPeakRAM = max(Peak_RAM_Used_MiB), .groups = "drop")

summary_total <- total_df %>%
    group_by(versionLabel, nCell) %>%
    summarise(medianTotalRuntime = median(totalRuntime),
              medianMaxPeakRAM = max(maxPeakRAM),
              .groups = "drop")

perCellTime <- ggplot(summary_total, aes(x = as.numeric(nCell), y = medianTotalRuntime, color = versionLabel)) +
    geom_line(size = 1) +
    geom_point() +
    labs(
        x = "Number of Cells (nCell)",
        y = "Median Total Runtime (sec)",
        color = "Version"
    ) +
    theme_minimal() +
    xlim(0, NA)

ggsave("Figs/report/perCellTime.pdf", perCellTime)
ggplot(summary_total, aes(x = as.numeric(nCell), y = medianMaxPeakRAM, color = versionLabel)) +
    geom_line(size = 1) +
    geom_point() +
    labs(
        x = "Number of Cells (nCell)",
        y = "Median Total Runtime (sec)",
        color = "Version"
    ) +
    theme_minimal()

#### MAX cells

summary_df_time <- combined_labeled %>%
    filter(nCell == 4000) %>%
    group_by(versionLabel, Function_Call) %>%
    summarise(medianRuntime = median(Elapsed_Time_sec))

summary_df_ram <- combined_labeled %>%
    filter(nCell == 4000) %>%
    group_by(versionLabel, Function_Call) %>%
    summarise(medianRAM = median(Total_RAM_Used_MiB))

summary_df_peak <- combined_labeled %>%
    filter(nCell == 4000) %>%
    group_by(versionLabel) %>%
    summarise(maxPeakRAM = max(Peak_RAM_Used_MiB))

time4000 <- ggplot(summary_df_time, aes(x = versionLabel, y = medianRuntime, colour = Function_Call)) +
    geom_col(aes(fill = Function_Call)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

ggsave("Figs/report/time4000.pdf", time4000, width = 10, height = 4)

RAM4000 <- ggplot(summary_df_ram, aes(x = versionLabel, y = medianRAM, colour = Function_Call)) +
    geom_col(aes(fill = Function_Call)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")

ggsave("Figs/report/RAM4000.pdf", RAM4000, width = 10, height = 5)

peak4000 <- ggplot(summary_df_peak, aes(x = versionLabel, y = maxPeakRAM)) +
    geom_col(aes(fill = versionLabel)) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Figs/report/peak4000.pdf", peak4000, width = 10, height = 5)

combined_labeled %>%
    filter(nCell == 4000) %>%
    group_by(versionLabel, Function_Call) %>%
    summarise(medianRuntime = median(Elapsed_Time_sec)) %>%
    ungroup() %>%
    group_by(versionLabel) %>%
    summarise(medianRuntime = sum(medianRuntime))

combined_labeled %>%
    filter(nCell == 4000) %>%
    group_by(versionLabel, Function_Call) %>%
    summarise(medianRAM = median(Total_RAM_Used_MiB)) %>%
    ungroup() %>%
    group_by(versionLabel) %>%
    summarise(ram = sum(medianRAM))
