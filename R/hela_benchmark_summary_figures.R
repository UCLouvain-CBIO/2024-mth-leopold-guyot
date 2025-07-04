library(tidyverse)

hela <- read.table("dataOutput/helaBenchmark/size_report.tsv", header = TRUE)
helaOpti <- read.table("dataOutput/helaBenchmark_opti/size_report.tsv", header = TRUE)

df_combined <- bind_rows(
    hela %>% mutate(type = "Vanilla"),
    helaOpti %>% mutate(type = "Optimized")
) %>%
    mutate(state = factor(state, levels = c("before", "after"))) %>%
    mutate(sizeTotal = sizeTotal / 1e9,
           sizeAssay = sizeAssay / 1e9,
           sizeRowData = sizeRowData / 1e9,
           sizeColData = sizeColData / 1e9)
df_long <- df_combined %>%
    pivot_longer(cols = starts_with("size"),
                 names_to = "SizeType",
                 values_to = "SizeValue")

# Keep only assay, rowdata, coldata
df_long_filtered <- df_long %>%
    filter(SizeType %in% c("sizeAssay", "sizeRowData", "sizeColData"))

# Compute medians
df_median <- df_long_filtered %>%
    group_by(state, type, SizeType) %>%
    summarise(medianSize = median(SizeValue), .groups = "drop")

df_median <- df_median %>%
    mutate(
        x_label = factor(paste(state, type, sep = " (") %>% paste0(")"),)
    ) %>%
    mutate(
        x_label = factor(x_label, levels = c("before (Vanilla)", "before (Optimized)",
                   "after (Vanilla)", "after (Optimized)"))
    )


size <- ggplot(df_median, aes(
    x = x_label,
    y = medianSize,
    fill = SizeType,
    color = type
)) +
    geom_bar(stat = "identity", alpha = 0.5, size = 1) +
    scale_fill_manual(values = c(
        "sizeAssay" = "darkorange",
        "sizeRowData" = "forestgreen",
        "sizeColData" = "steelblue"
    )) +
    labs(x = "Object state and optimisation",
         y = "Median Size (GB)",
         fill = "Component") +
    theme_minimal()



ggsave("Figs/report/HeLa_medianSizeComponents.pdf", size, width = 10, height = 5)



# Define the directories
dirs <- c("dataOutput/helaBenchmark/memoryOutput/", "dataOutput/helaBenchmark_opti/memoryOutput/")

# Function to read a file and add metadata columns
read_file_with_metadata <- function(file_path, dir_name) {
    read_csv(file_path) %>%
        mutate(
            Directory = dir_name,
            Filename = basename(file_path)
        )
}
df <- data.frame()
for (dir in dirs) {
    files <- list.files(path = dir, full.names = TRUE)
    for (file in files) {
        df <- rbind(df, read_file_with_metadata(file, dir))
    }
}

df <- df %>%
    mutate(
        Step = case_when(
            grepl("joining", Function_Call) ~ "Joining",
            grepl("aggregationToPep", Function_Call) ~ "Aggregation",
            TRUE ~ Function_Call
        )
    ) %>%
    filter(Step %in% c("Joining", "Aggregation")) %>%
    mutate(
        Version = case_when(
            grepl("/helaBenchmark/", Directory) ~ "Vanilla",
            grepl("/helaBenchmark_opti/", Directory) ~ "Optimised",
            TRUE ~ Directory
        ))


# Plot all steps
stepTime <- ggplot(df, aes(x = Step, y = Elapsed_Time_sec, color = Version)) +
    geom_boxplot() +
    labs(
         x = "Operation",
         y = "Elapsed Time (s)",
         fill = "Type") +
    theme_minimal()

ggsave("Figs/report/HeLa_runtimeSteps.pdf", stepTime, width = 10, height = 5)
