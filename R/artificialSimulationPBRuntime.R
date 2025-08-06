library(tidyverse)

res <- readRDS("dataOutput/artificialSimPB/tprfdrRes.rds")
runtime <- list()
for (name in names(res)) {
    curr <- res[[name]][["timings"]]
    splited <- str_split_1(name, "_")
    runtime[[name]] <- data.frame(parameters = name,
                                  cellPerComb = splited[[1]],
                                  introducedShift = splited[[2]])
    for (step in names(curr)) {
        runtime[[name]][[step]] <- curr[[step]]
    }

}

combined <- do.call(rbind, runtime) %>%
    mutate(pseudobulkMean = pseudobulkMean + aggMean,
           pseudobulkMedian = pseudobulkMedian + aggMedian,
           pseudobulkRobustSummary = pseudobulkRobustSummary + aggRobust,
           pseudobulkSum = pseudobulkSum + aggSum) %>%
    pivot_longer(cols = 4:16, names_to = "methods", values_to = "runtime") %>%
    mutate(cellPerCombNum = as.numeric(sub("nCell", "", cellPerComb))) %>%
    mutate(nCell = cellPerCombNum * 16 * 5)


combined %>%
    filter(!(methods %in% c("simulateCellPatientData", "addTreatmentEffect",
                          "compute_performance",
                          "aggMean", "aggMedian", "aggRobust", "aggSum"))) %>%
    ggplot(aes(x = nCell, y = runtime, color = methods)) +
        geom_boxplot(
            aes(group = interaction(nCell, methods), color = methods),
            alpha = 0.7,
            position = position_dodge(width = 1),
            size = 0.8,
            width = 2000
        ) +
        stat_summary(
            fun = median,
            geom = "line",
            aes(group = methods),
            size = 1
        ) +
        ylim(0, NA)+
        ylab("Time (s)") +
        xlab("Number of cells")

ggsave("Figs/report/runtimeSim.pdf")
