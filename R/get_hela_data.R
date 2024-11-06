library(curl)
library(rpx)
library(SummarizedExperiment)
px <- PXDataset("PXD042233")
url <- pxurl(px)
files <- grep("Orbitrap.*zip", pxfiles(px), value = TRUE)
files <- files[1:10]
sucess <- multi_download(file.path(url, files), destfiles = file.path("data/hela/rawData", files))
