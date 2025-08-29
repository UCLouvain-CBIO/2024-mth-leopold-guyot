library(tidyverse)

res <- readRDS("dataOutput/artificialSimPB/tprfdrRes.rds")
for (name in names(res)) {
    splited <- str_split_1(name, "_")
    res[[name]]$performance$parameters <- name
    res[[name]]$performance$cellPerComb <- splited[[1]]
    res[[name]]$performance$introducedShift <- splited[[2]]
    res[[name]]$performance$replicate <- splited[[3]]
}

combined <- do.call(rbind, lapply(res, function(x) x$performance))

combined <- combined %>%
    filter(thr %in% c(0.01, 0.05, 0.1)) %>%
    mutate(cellPerComb = factor(cellPerComb, levels = c("nCell10", "nCell25", "nCell50", "nCell100")))

subCombined <- combined %>%
    filter(thr == 0.05,
           cellPerComb == "nCell50",
           introducedShift == "shift0.15")
fdrMod <- lm(formula = FDR ~ method, data = subCombined)
write.csv(x = summary(fdrMod)$coefficients, "dataOutput/artificialSimPB/modFdr.csv")

tprMod <- lm(formula = TPR ~ method, data = subCombined)
write.csv(x = summary(tprMod)$coefficients, "dataOutput/artificialSimPB/modTpr.csv")

combined_summary <- combined %>%
    group_by(method, thr, cellPerComb, introducedShift) %>%
    summarise(
        TPR_mean = mean(TPR),
        TPR_sd = sd(TPR),
        FDR_mean = mean(FDR),
        FDR_sd = sd(FDR),
        .groups = "drop"
    )

combined_summary %>%
    filter(introducedShift == "shift2") %>%
    ggplot(aes(x = FDR_mean, y = TPR_mean, color = cellPerComb)) +
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

ggsave("Figs/artiSimPerCell_errorbar.pdf", width = 7, height = 5)


combined_summary %>%
    filter(cellPerComb == "nCell50") %>%
    ggplot(aes(x = FDR_mean, y = TPR_mean, color = method)) +
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

ggsave("Figs/artiSimPerShift_errorbar.pdf", width = 7, height = 5)

combined_latex_shift <- combined_summary %>%
    filter(thr == 0.05,
           cellPerComb == "nCell50") %>%
    mutate(
        TPR = sprintf("%.2f ± %.2f", TPR_mean, TPR_sd),
        FDR = sprintf("%.2f ± %.2f", FDR_mean, FDR_sd)
    ) %>%
    select(method, thr, cellPerComb, introducedShift, TPR, FDR) %>%
    arrange(introducedShift, thr)

write.csv(x = combined_latex_shift, "dataOutput/artificialSimPB/combined_latex_shift.csv")

combined_latex_cell <- combined_summary %>%
    filter(thr == 0.05,
           introducedShift == "shift0.15") %>%
    mutate(
        TPR = sprintf("%.2f ± %.2f", TPR_mean, TPR_sd),
        FDR = sprintf("%.2f ± %.2f", FDR_mean, FDR_sd)
    ) %>%
    select(method, thr, cellPerComb, introducedShift, TPR, FDR) %>%
    arrange(cellPerComb, thr)

write.csv(x = combined_latex, "dataOutput/artificialSimPB/combined_latex_cell.csv")
