library(VerR)
library(ggplot2)

envCreate("expList", packages = c("leopoldguyot/QFeatures@master",
                                  "bioc::scpdata",
                                  "bioc::scp"))

envCreate("list", packages = c("leopoldguyot/QFeatures@listBench",
                               "bioc::scpdata",
                               "bioc::s@Manual{,
    title = {ExperimentHub: Client to access ExperimentHub resources},
    author = {Martin Morgan and Lori Shepherd},
    year = {2025},
    note = {R package version 2.16.0},
    url = {https://bioconductor.org/packages/ExperimentHub},
    doi = {10.18129/B9.bioc.ExperimentHub},
  }cp"))
envCreate("base", packages = c("rformassspectrometry/QFeatures@master",
                               "bioc::scpdata",
                               "bioc::scp"))

res <- benchInEnv(expr = {
    aggregateFeatures(data, fcol = "peptide",i = 1:40, name = paste0("pep", 1:40),
                                                                     fun = colMeans)
    }, setup = {
        library("QFeatures")
        library("scpdata")
        data <- scpdata::gregoire2023_mixCTRL()
})
df <- data.frame(time = unlist(res),
                 version = c(rep(names(res)[[1]], 3),
                             rep(names(res)[[2]], 3),
                             rep(names(res)[[3]], 3)))
ggplot2::ggplot(df, aes(x = version, y = time)) + geom_boxplot() + ylim(0,200)


data <- scpdata::gregoire2023_mixCTRL()

library(microbenchmark)
library(scpdata)
library(QFeatures)
data <- scpdata::gregoire2023_mixCTRL()[,, 1:40]

microbenchmark(as.list(experiments(data)),
               as.list(as(data, "List")),
               lapply(names(data), function(i) data[[i]]),
               times = 1)
