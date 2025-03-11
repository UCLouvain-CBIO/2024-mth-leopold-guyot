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

agg <- readData("dataOutput/optimisationsBench/agg") %>%
    mutate(optimisation = "agg")
all <- readData("dataOutput/optimisationsBench/all") %>%
    mutate(optimisation = "all")
subset <- readData("dataOutput/optimisationsBench/subset") %>%
    mutate(optimisation = "subset")
base <- readData("dataOutput/optimisationsBench/base") %>%
    mutate(optimisation = "base")
allSubset <- readData("dataOutput/optimisationsBench/allSubset") %>%
    mutate(optimisation = "allSubset")
subsetRowData <- readData("dataOutput/optimisationsBench/subsetRowData") %>%
    mutate(optimisation = "subsetRowData")
df <- rbind(base, agg, subset, subsetRowData, allSubset, all)

df %>%
    mutate(optimisation = factor(optimisation, levels = c("base", "agg", "subset", "subsetRowData", "all", "allSubset"))) %>%
    filter(Function_Call == "leduc<-leducAggPSM(leduc)") %>%
    group_by(optimisation, replicate) %>%
    summarise(TotalTime = sum(Elapsed_Time_sec)) %>%
    ggplot(aes(x = optimisation, y = TotalTime, color = optimisation)) +
        geom_boxplot()
