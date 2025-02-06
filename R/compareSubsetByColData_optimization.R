library(callr)
library(dplyr)

runBench <- function(git, rep = 3L) {
    remotes::install_github(git, force = TRUE)
    library(bench)
    library(MultiAssayExperiment)
    library(QFeatures)
    library(scpdata)
    qfeatures <- leduc2022_pSCoPE()
    colData(qfeatures)$filterBench <- rnorm(nrow(colData(qfeatures)), mean = 1, sd = 1)
    bench::mark(
        result = qfeatures[, qfeatures$filterBench > 1, ],
        iterations = rep
    )
}

old <- callr::r(
    func = runBench,
    args = list(git = "waldronlab/MultiAssayExperiment@325c85ca26582027bcdd643df673a06ae649a36b", rep = 3L)
)

new <- callr::r(
    func = runBench,
    args = list(git = "leopoldguyot/MultiAssayExperiment@subsetByColData_optimization", rep = 3L)
)

print(old)
print(new)

## More advanced:

runPress <- function(git, nAssays, rep = 3L) {
    remotes::install_github(git, force = TRUE)
    library(bench)
    library(MultiAssayExperiment)
    library(QFeatures)
    library(scpdata)
    qfeatures <- leduc2022_pSCoPE()
    colData(qfeatures)$filterBench <- rnorm(nrow(colData(qfeatures)), mean = 1, sd = 1)
    bench::press(
        nAssays = nAssays,
        rep = 1:3,
        {
            subqfeatures <- qfeatures[, , 1:nAssays]
            results <- bench::mark(
                bracket = subqfeatures[, subqfeatures$filterBench > 1, ],
                subset = subsetByColData(subqfeatures, subqfeatures$filterBench > 1)
            )

            tibble::tibble(
                expression = as.character(results$expression),
                median_time = results$median,
                mem_alloc = results$mem_alloc
            )
        }
    )
}

new <- callr::r(
    func = runPress,
    args = list(
        git = "leopoldguyot/MultiAssayExperiment@subsetByColData_optimization",
        nAssays = c(5, 20, 50, 80, 120),
        rep = 3L
    )
)

old <- callr::r(
    func = runPress,
    args = list(
        git = "waldronlab/MultiAssayExperiment@325c85ca26582027bcdd643df673a06ae649a36b",
        nAssays = c(5, 20, 50, 80, 120),
        rep = 3L
    )
)


old <- old %>%
    mutate(version = "Old")

new <- new %>%
    mutate(version = "New")

combined <- bind_rows(old, new) %>%
    mutate(
        median_time = as.numeric(median_time), # in seconds
        mem_alloc = as.numeric(sub(
            "MB", "",
            as.character(mem_alloc)
        ))
    )


write.csv(combined, file = "dataOutput/benchmark_subsetByColData_opti.csv")


combined <- read.csv("dataOutput/benchmark_subsetByColData_opti.csv")

library(tidyr)
library(ggplot2)
library(plotly)
plotTime <- combined %>%
    mutate(
        nCol = nAssays * 18,
        nAssays = as.factor(nAssays)
    ) %>%
    ggplot(aes(x = nCol, y = median_time, color = version)) +
    xlab("Total number of columns") +
    ylab("RunTime (s)") +
    geom_point() +
    geom_smooth() +
    facet_wrap(~expression)

ggplotly(plotTime)

plotMem <- combined %>%
    mutate(
        nCol = nAssays * 18,
        nAssays = as.factor(nAssays),
        mem_alloc_MB = mem_alloc / (1024**2)
    ) %>%
    ggplot(aes(x = nCol, y = mem_alloc_MB, color = version)) +
    xlab("Total number of columns") +
    ylab("RAM used (MB)") +
    geom_point() +
    geom_smooth() +
    facet_wrap(~expression)

ggplotly(plotMem)
