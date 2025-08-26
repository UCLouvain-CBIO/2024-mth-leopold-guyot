library("tidyverse")

noSubSCE <- read.csv(file = "dataOutput/optimisationsBench/noSubsetRowDataSCE.csv")
subSCE <- read.csv(file = "dataOutput/optimisationsBench/subsetRowDataSCE.csv")
noSubSE <- read.csv(file = "dataOutput/optimisationsBench/noSubsetRowDataSE.csv")
subSE <- read.csv(file = "dataOutput/optimisationsBench/subsetRowDataSE.csv")

noSubSE <- noSubSE %>%
    separate(X, into = c("nCell", "rep"), sep = "\\.") %>%
    mutate(subsetRowData = FALSE) %>%
    mutate(SE = TRUE)

subSE <- subSE %>%
    separate(X, into = c("nCell", "rep"), sep = "\\.") %>%
    mutate(subsetRowData = TRUE) %>%
    mutate(SE = TRUE)

noSubSCE <- noSubSCE %>%
    separate(X, into = c("nCell", "rep"), sep = "\\.") %>%
    mutate(subsetRowData = FALSE) %>%
    mutate(SE = FALSE)

subSCE <- subSCE %>%
    separate(X, into = c("nCell", "rep"), sep = "\\.") %>%
    mutate(subsetRowData = TRUE) %>%
    mutate(SE = FALSE)

combined <- rbind(noSubSCE, subSCE, noSubSE, subSE) %>%
    mutate(versionComplete = paste(version, SE, subsetRowData, sep = "_"))

summary_df_time <- combined %>%
    filter(nCell == 4000) %>%
    group_by(version, subsetRowData, SE, Function_Call) %>%
    summarise(medianRuntime = median(Elapsed_Time_sec)) %>%
    ungroup() %>%
    group_by(version, subsetRowData, SE) %>%
    summarize(medianRuntime = sum(medianRuntime)) %>%
    mutate(version = factor(version, levels = c("base", "subsetByColData", "aggregation", "allOpti")))

ggplot(summary_df_time, aes(x = version, y = medianRuntime)) +
    geom_col() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(rows = vars(subsetRowData), cols = vars(SE), labeller = label_both)

