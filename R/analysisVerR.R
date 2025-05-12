library("tidyverse")

noSub <- read.csv(file = "dataOutput/optimisationsBench/noSubsetRowData.csv")
sub <- read.csv(file = "dataOutput/optimisationsBench/subsetRowData.csv")

noSub <- noSub %>%
    separate(X, into = c("nAssay", "rep"), sep = "\\.") %>%
    mutate(numNAssay = as.numeric(nAssay)) %>%
    mutate(nAssay = factor(nAssay, levels = sort(unique(as.numeric(nAssay))))) %>%
    mutate(sub = "noSub")

sub <- sub %>%
    separate(X, into = c("nAssay", "rep"), sep = "\\.") %>%
    mutate(numNAssay = as.numeric(nAssay)) %>%
    mutate(nAssay = factor(nAssay, levels = sort(unique(as.numeric(nAssay))))) %>%
    mutate(sub = "sub")

combined <- rbind(sub, noSub) %>%
    filter(nAssay == 60) %>%
    group_by(sub, Function_Call, version, numNAssay) %>%
    summarise("medianTime" = median(Elapsed_Time_sec)) %>%
    mutate(opti = paste(sub, version, sep = "_")) %>%
    filter(opti %in% c("noSub_aggregation",
                            "noSub_base",
                            "noSub_subsetByColData",
                            "sub_base",
                            "sub_allOpti")) %>%
    mutate(opti = factor(opti, levels = c(
        "noSub_base",
        "noSub_subsetByColData",
        "noSub_aggregation",
        "sub_base",
        "sub_allOpti"
    )))
combinedRAM <- rbind(sub, noSub) %>%
    filter(nAssay == 60) %>%
    group_by(sub, Function_Call, version, numNAssay) %>%
    summarise("peakRam" = max(Peak_RAM_Used_MiB)) %>%
    mutate(opti = paste(sub, version, sep = "_")) %>%
    filter(opti %in% c("noSub_aggregation",
                       "noSub_base",
                       "noSub_subsetByColData",
                       "sub_base",
                       "sub_allOpti")) %>%
    mutate(opti = factor(opti, levels = c(
        "noSub_base",
        "noSub_subsetByColData",
        "noSub_aggregation",
        "sub_base",
        "sub_allOpti"
    )))

combinedUsedRAM <- rbind(sub, noSub) %>%
    filter(nAssay == 60) %>%
    group_by(sub, Function_Call, version, numNAssay) %>%
    summarise("usedRam" = median(Total_RAM_Used_MiB)) %>%
    mutate(opti = paste(sub, version, sep = "_")) %>%
    filter(opti %in% c("noSub_aggregation",
                       "noSub_base",
                       "noSub_subsetByColData",
                       "sub_base",
                       "sub_allOpti")) %>%
    mutate(opti = factor(opti, levels = c(
        "noSub_base",
        "noSub_subsetByColData",
        "noSub_aggregation",
        "sub_base",
        "sub_allOpti"
    )))

ggplot(combined, aes(x = opti, y = medianTime, colour = Function_Call))+
    geom_col(aes(fill = Function_Call))
ggplot(combinedRAM, aes(x = opti, y = peakRam, colour = Function_Call))+
    geom_col(aes(fill = Function_Call))
ggplot(combinedUsedRAM, aes(x = opti, y = usedRam, colour = Function_Call))+
    geom_col(aes(fill = Function_Call))
