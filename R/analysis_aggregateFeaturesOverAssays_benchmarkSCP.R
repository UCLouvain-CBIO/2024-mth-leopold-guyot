library(tidyverse)
library(plotly)

folder <- "dataOutput/aggregateBenchSCP"
files <- list.files(folder, full.names = TRUE)

df <- map_df(files, function(file) {
    curr <- read.csv(file)

    base <- basename(file)

    numbers <- str_extract_all(base, "\\d+")[[1]]
    first_num <- as.integer(numbers[1])
    second_num <- as.integer(numbers[2])

    curr <- curr %>%
        mutate(
            size = first_num,
            replicate = second_num
        )

    return(curr)
})

head(df)

lm_false <- lm(Elapsed_Time_sec ~ size, df[df$Function_Call == 'aggPSM(qfeat,FALSE)',])
lm_true <- lm(Elapsed_Time_sec ~ size, df[df$Function_Call == 'aggPSM(qfeat,TRUE)',])

slope_false <- round(coef(lm_false)[2], 4)
slope_true <- round(coef(lm_true)[2], 4)

plot <- ggplot(data = df, aes(x = size, y = Elapsed_Time_sec, color = Function_Call)) +
    geom_point() +
    geom_smooth(method = "lm", size = 0.5) +
    xlab("size (number of cells, TMT-16)")+
    ylab("Runtime (s)")+
    annotate("text", x = max(df$size) * 0.3, y = max(df$Elapsed_Time_sec) * 0.9,
             label = paste0("Slope (FALSE): ", slope_false), color = "red", hjust = 0) +
    annotate("text", x = max(df$size) * 0.3, y = max(df$Elapsed_Time_sec) * 0.85,
             label = paste0("Slope (TRUE): ", slope_true), color = "blue", hjust = 0)

ggsave(filename = "Figs/benchmarkAggregateFeaturesOverAssays.png", plot)
