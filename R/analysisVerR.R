library("tidyverse")

noSub <- read.csv(file = "dataOutput/optimisationsBench/noSubsetRowData.csv")
sub <- read.csv(file = "dataOutput/optimisationsBench/subsetRowData.csv")

noSub <- noSub %>%
    separate(X, into = c("nAssay", "rep"), sep = "\\.") %>%
    mutate(numNAssay = as.numeric(nAssay)) %>%
    mutate(nAssay = factor(nAssay, levels = sort(unique(as.numeric(nAssay)))))

sub <- sub %>%
    separate(X, into = c("nAssay", "rep"), sep = "\\.") %>%
    mutate(numNAssay = as.numeric(nAssay)) %>%
    mutate(nAssay = factor(nAssay, levels = sort(unique(as.numeric(nAssay)))))

ggplot(noSub, aes(x = numNAssay, y = Elapsed_Time_sec, colour = version)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE, aes(group = version)) +
    labs(x = "nAssay (numeric)", y = "Elapsed Time (sec)", title = "Performance Growth for noSub") +
    facet_wrap(~Function_Call)

ggplot(sub, aes(x = numNAssay, y = Elapsed_Time_sec, colour = version)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, aes(group = version)) +
    labs(x = "nAssay (numeric)", y = "Elapsed Time (sec)", title = "Performance Growth for sub")+
    facet_wrap(~Function_Call)

