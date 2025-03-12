library(tidyverse)

readData <- function(path) {
    inputDir <- path
    res <- read.csv(
        file = file.path(inputDir, "microbenchres.csv")
    )
    res
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
    ggplot(aes(x = optimisation, y = time, color = optimisation)) +
        geom_boxplot() + ylim(0, max(df$time))
