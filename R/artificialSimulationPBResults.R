library(tidyverse)

res <- readRDS("dataOutput/artificialSimPB/tprfdrRes.rds")
for (name in names(res)) {
    splited <- str_split_1(name, "_")
    res[[name]]$performance$parameters <- name
    res[[name]]$performance$cellPerComb <- splited[[1]]
    res[[name]]$performance$introducedShift <- splited[[2]]
}

combined <- do.call(rbind, lapply(res, function(x) x$performance))

combined <- combined %>%
    filter(thr %in% c(0.01, 0.05, 0.1)) %>%
    mutate(cellPerComb = factor(cellPerComb, levels = c("nCell10", "nCell25", "nCell50", "nCell100")))

## Per cell per combination
combinedThr01 <- combined %>%
    filter(thr %in% c(0.05),
           introducedShift == "shift0.15")

combined %>%
    filter(introducedShift == "shift0.15") %>%
    ggplot(aes(x = FDR, y = TPR, color = cellPerComb)) +
    geom_vline(
        xintercept = c(0.01, 0.05, 0.1),
        linetype = "dashed", color = "grey50", linewidth = 0.3
    ) +
    geom_point(size = 1, alpha = 0.8) +
    geom_line(size = 0.5) +
    scale_x_continuous(
        trans = "sqrt",
        limits = c(0, 1),
        breaks = c(0.01, 0.2, 0.6, 1),
        labels = scales::label_number()
    ) +
    scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.2)
    ) +
    theme_minimal() +
    labs(
        x = "False Discovery Rate (sqrt scale)",
        y = "True Positive Rate",
        color = "Number of Cells per Combination"
    ) + facet_wrap(~method) +
    theme(
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.position = "bottom"
    )

ggsave("Figs/artiSimPerCell.pdf", width = 7, height = 5)


## Per shift plot

combinedThr50 <- combined %>%
    filter(thr %in% c(0.05),
           cellPerComb == "nCell50")

combined %>%
    filter(cellPerComb == "nCell50") %>%
    ggplot(aes(x = FDR, y = TPR, color = method)) +
        geom_vline(
            xintercept = c(0.01, 0.05, 0.1),
            linetype = "dashed", color = "grey50", linewidth = 0.3
        ) +
        geom_point(size = 1, alpha = 0.8) +
        geom_line(size = 0.5) +
        scale_x_continuous(
            trans = "sqrt",
            limits = c(0, 1),
            breaks = c(0.01, 0.2, 0.6, 1),
            labels = scales::label_number()
        ) +
        scale_y_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, 0.2)
        ) +
        theme_minimal() +
        labs(
            x = "False Discovery Rate (sqrt scale)",
            y = "True Positive Rate",
            color = "Method"
        ) +
        facet_wrap(~introducedShift) +
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            legend.position = "bottom"
        )

ggsave("Figs/artiSimPerShift.pdf", width = 7, height = 5)
