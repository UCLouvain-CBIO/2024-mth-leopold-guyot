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


ggplot(df_combined, aes(x = state, y = sizeTotal, fill = type)) +
    geom_bar(position = "dodge", stat = "identity") +
    ylab("Total Size (GB)")
ggsave(filename = "Figs/HeLa_sizeTotal.png")
df_long <- df_combined %>%
    pivot_longer(cols = starts_with("size"),
                 names_to = "SizeType",
                 values_to = "SizeValue")


# Create the plot
ggplot(df_long, aes(x = state, y = SizeValue, fill = type)) +
    geom_bar(position = "dodge", stat = "identity") +
    facet_wrap(~ SizeType, scales = "free_y") +
    theme_minimal() +
    labs(title = "Comparison of Different Size Metrics Before and After",
         x = "State",
         y = "Size Value",
         fill = "Type")


### Steps


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

# Convert numeric columns
df <- df %>%
    mutate(
        Elapsed_Time_sec = as.numeric(Elapsed_Time_sec),
        Total_RAM_Used_MiB = as.numeric(Total_RAM_Used_MiB),
        Peak_RAM_Used_MiB = as.numeric(Peak_RAM_Used_MiB),
        Version = as.factor(Directory)
    )
levels(df$Version) <- c("Optimized", "Vanilla")
df %>%
    group_by(Filename, Version) %>%
    summarise(totalTime = sum(Elapsed_Time_sec)) %>%
    ggplot(aes(x = Version, y = totalTime)) +
        geom_boxplot() +
        ylab("Total Time (s)")
ggsave(filename = "Figs/HeLa_timeTotal.png")
df %>%
    filter(Function_Call == "hela<-joining(hela,assaysNames)") %>%
    group_by(Filename, Version) %>%
    ggplot(aes(x = Version, y = Elapsed_Time_sec)) +
    geom_boxplot() +
    ylab("Joining Time (s)")
ggsave(filename = "Figs/HeLa_timeJoin.png")

df %>%
    filter(Function_Call == "hela<-aggregationToPep(hela,assaysNames)") %>%
    group_by(Filename, Version) %>%
    ggplot(aes(x = Version, y = Elapsed_Time_sec)) +
    geom_boxplot() +
    ylab("Aggregation Time (s)")
ggsave(filename = "Figs/HeLa_timeAgg.png")
