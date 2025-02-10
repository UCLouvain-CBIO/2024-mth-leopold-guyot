library(curl)
library(rpx)
library(QFeatures)
px <- PXDataset("PXD042233")
url <- pxurl(px)
files <- grep("Orbitrap.*zip", pxfiles(px), value = TRUE)
files <- files[1:500]

sampleAnnot <- pxget(px, "Experimental-Design.sdrf.tsv")
sampleAnnot <- DataFrame(read.delim(sampleAnnot))
sampleAnnot$runCol <- sampleAnnot$source.name
sampleAnnot$quantCols <- rep("Intensity", nrow(sampleAnnot))
rownames(sampleAnnot) <- sampleAnnot$source.name

if (dir.exists("data/hela/rawData")) {
    unlink("data/hela/rawData", recursive = TRUE)
}

dir.create("data/hela/rawData")

sucess <- multi_download(file.path(url, files), destfiles = file.path("data/hela/rawData", files))
for (zip in sucess$destfile) {
    baseName <- basename(sub(".zip", "", zip))
    unzip(zip, files = file.path(baseName, "evidence.txt"), exdir = "data/hela/rawData")
    file.remove(zip)
}
SEList <- lapply(list.files("data/hela/rawData"), function(dir) {
    
    fullPath <- file.path("data/hela/rawData", dir, "evidence.txt")
    quantTable <- read.delim(fullPath)
    quantTable$PSMid <- paste0(dir, "_PSM", seq(nrow(quantTable)))
    quantTable$runCol <- rep(dir, nrow(quantTable))
    rownames(quantTable) <- quantTable$PSMid
    SE <- readSummarizedExperiment(quantTable,
        quantCols = "Intensity",
        fnames = "PSMid"
    )
    colnames(SE) <- dir
    SE
})

names(SEList) <- list.files("data/hela/rawData")
qfeatures <- QFeatures(SEList, colData = sampleAnnot)

saveRDS(qfeatures, file = "dataOutput/hela/qfeatures_PSM.rds")
