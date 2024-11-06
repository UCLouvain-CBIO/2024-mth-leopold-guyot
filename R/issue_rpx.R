library(rpx)
px <- PXDataset("PXD042233")
maxzip <- grep("Orbitrap.*zip", pxfiles(px), value = TRUE)
maxzip <- maxzip[1:100]
base <- pxurl(px)
test <- lapply(maxzip, function(x) RCurl::url.exists(file.path(base, x)))
print(sum(as.logical(test))) # should be 100

BiocFileCache::removebfc()
