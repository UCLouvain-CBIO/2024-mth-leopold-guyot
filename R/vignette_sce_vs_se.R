library(tidyverse)

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

benchRes %>%
    filter(nCell == 1000) %>%
    ggplot(aes(x = Function_Call, y = Elapsed_Time_sec, colour = type)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
