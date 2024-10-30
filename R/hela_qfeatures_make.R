library(rpx)
library(QFeatures)

px <- PXDataset("PXD042233")
test <- read.delim("2013_04_03_16_54_Q-Exactive-Orbitrap_1/evidence.txt")
maxzip <- grep("Orbitrap.*zip", pxfiles(px), value = TRUE)
sampleAnnot <- pxget(px, "Experimental-Design.sdrf.tsv")
sampleAnnot <- read.delim(sampleAnnot)
sampleAnnot$runCol <- sampleAnnot$source.name
sampleAnnot$quantCols <- rep("Intensity", nrow(sampleAnnot))
maxzip <- maxzip[1:10]


matrixList <- lapply(maxzip, function(x){
    zipFile <- pxget(px, x)
    baseName <- sub(".zip", "", x)
    unzip(zipFile, files = file.path(baseName, "evidence.txt"), exdir = "tmp")
    quantTable <- read.delim(file.path("tmp", baseName, "evidence.txt"))
    quantTable$runCol <- rep(baseName, nrow(quantTable))
    quantTable
})

unlink("tmp", recursive = TRUE)
quantTable <- do.call(rbind, matrixList)

qfeatures <- readQFeatures(quantTable, colData = sampleAnnot, runCol = "runCol")

saveRDS(qfeatures, file = "dataOutput/hela/qfeatures_PSM.rds")
