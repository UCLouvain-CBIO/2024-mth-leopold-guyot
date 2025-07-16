library(tidyverse)
library(patchwork)

readData <- function(path) {
    inputDir <- path
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
    benchmarkDF[-irrelevantIndices, ]
}

SCE <- readData("dataOutput/vignetteBenchmarkSCE/")
SE <- readData("dataOutput/vignetteBenchmarkSE/")

SCE <- mutate(SCE, type = "sce")
SE <- mutate(SE, type = "se")

benchRes <- rbind(SCE, SE)
funcOrder <- benchRes %>%
    filter(nCell == 1000) %>%
    group_by(Function_Call) %>%
    summarize(median_time_sec = median(Elapsed_Time_sec)) %>%
    arrange(desc(median_time_sec)) %>%
    pull(Function_Call)

time <- benchRes %>%
    filter(nCell == 1000) %>%
    mutate(Function_Call = factor(Function_Call, levels = funcOrder)) %>%
    ggplot(aes(x = Function_Call, y = Elapsed_Time_sec, colour = type)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

perct <- benchRes %>%
    filter(nCell == 1000) %>%
    mutate(Function_Call = factor(Function_Call, levels = funcOrder)) %>%
    group_by(Function_Call, type) %>%
    summarize(median_time_sec = median(Elapsed_Time_sec), .groups = "drop") %>%
    pivot_wider(names_from = type, values_from = median_time_sec) %>%
    mutate(percentage_diff = ((se - sce) / sce) * 100) %>%
    ggplot(aes(x = Function_Call, y = percentage_diff)) +
    geom_col() +
    ylab("Median decrease in % when using SE") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

time <- time + theme(axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank())
combined_plot <- (time / perct) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
ggsave("Figs/SCEvsSE_vignette.pdf", combined_plot,
       height = 10,
       width = 15)
