library(tidyverse)
library(plotly)

folder <- "dataOutput/aggregateBench"
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
            replicate = second_num,
            Function_Call = str_extract(Function_Call, "^[^(]+")  # Extract part before '('
        )

    return(curr)
})

head(df)

ggplotly(ggplot(data = df, aes(x = size, y = Elapsed_Time_sec, color = Function_Call))+
    geom_point() +geom_smooth(method = "lm", size = 0.5))
