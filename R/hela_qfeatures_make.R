library(rpx)
library(QFeatures)
library(SummarizedExperiment)

px <- PXDataset("PXD042233")
maxzip <- grep("Orbitrap.*zip", pxfiles(px), value = TRUE)
sampleAnnot <- pxget(px, "Experimental-Design.sdrf.tsv")
sampleAnnot <- DataFrame(read.delim(sampleAnnot))
sampleAnnot$runCol <- sampleAnnot$source.name
sampleAnnot$quantCols <- rep("Intensity", nrow(sampleAnnot))
rownames(sampleAnnot) <- sampleAnnot$source.name
maxzip <- maxzip[1:100]

SEList <- lapply(maxzip, function(x){
    zipFile <- pxget(px, x)
    baseName <- sub(".zip", "", x)
    unzip(zipFile, files = file.path(baseName, "evidence.txt"), exdir = "tmp")
    quantTable <- read.delim(file.path("tmp", baseName, "evidence.txt"))
    quantTable$PSMid <- paste0(baseName, "_PSM", seq(nrow(quantTable)))
    quantTable$runCol <- rep(baseName, nrow(quantTable))
    rownames(quantTable) <- quantTable$PSMid
    SE <- readSummarizedExperiment(quantTable,
                                   quantCols = "Intensity",
                                   fnames = "PSMid")
    # colnames(SE) <- baseName
    colData(SE) <- sampleAnnot[sampleAnnot$source.name == baseName, ]
    colData(SE) <- NULL
    SE
})
unlink("tmp", recursive = TRUE)

names(SEList) <- sub(".zip", "", maxzip)
qfeatures <- QFeatures(SEList, colData = sampleAnnot)

saveRDS(qfeatures, file = "dataOutput/hela/qfeatures_PSM.rds")
