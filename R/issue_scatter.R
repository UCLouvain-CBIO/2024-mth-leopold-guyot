# Issue:
# When scatter is loaded normalize with scp qfeatures is not working properly

# load data

library(scp)
data("scp1", package = "scp")

# Working
test1 <- QFeatures::normalize(scp1,
                              i = "peptides",
                              method = "div.median",
                              name = "peptides_norm1")



library(scater)

# Working
test2 <- QFeatures::normalize(scp1,
                              i = "peptides",
                              method = "div.median",
                              name = "peptides_norm1")

detach("package:scater", unload = TRUE)
library(scater)

# Not working
test3 <- QFeatures::normalize(scp1,
                              i = "peptides",
                              method = "div.median",
                              name = "peptides_norm1")
