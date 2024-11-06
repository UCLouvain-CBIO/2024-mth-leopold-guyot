library(curl)
library(rpx)
library(SummarizedExperiment)
px <- PXDataset("PXD042233")
url <- pxurl(px)
files <- grep("Orbitrap.*zip", pxfiles(px), value = TRUE)
files <- files[1:10]
for (file in files) {
    curl_download(url = file.path(url, file), destfile = file)
    baseName <- sub(".zip", "", file)
    unzip(file, files = file.path(baseName, "evidence.txt"), exdir = "data/hela/rawData")
    file.remove(file)
}
