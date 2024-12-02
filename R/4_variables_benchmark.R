library(scpdata)
source(file.path("R", "generate_data.R"))
source(file.path("R", "utils.R"))

# Variables:

destPath <- file.path("dataOutput", "4_variables_benchmark")

nCells <- c(1, 16, 32)

nFeats <- c(1000, 2000, 4000)

nAssays <- c(8, 16, 32, 64, 128, 256)

nCols <- c(10, 50, 100)

base <- scpdata::brunner2022()
set.seed(123)

unlink(destPath, recursive = TRUE)
dir.create(destPath)
write.table(data.frame(nCell = integer(),
                       nFeat = integer(),
                       nAssay = integer(),
                       nCol = integer(),
                       sizeTotal = numeric(),
                       sizeAssay = numeric(),
                       sizeRowData = numeric(),
                       sizeColData = numeric()),
            file = file.path(destPath, "qfeatures_size_report.tsv"),
            append = FALSE,
            col.names = TRUE,
            row.names = FALSE)

for (nCell in nCells) {
    for (nFeat in nFeats) {
        for (nAssay in nAssays) {
            for (nCol in nCols) {
                print(paste0("Starting: nCell ",
                             nCell, " nFeat ",
                             nFeat, " nAssay ",
                             nAssay, " nCol ",
                             nCol))
                qfeatures <- generate4VarData(base,
                                 nCell = nCell,
                                 nFeat = nFeat,
                                 nAssay = nAssay,
                                 nCol = nCol
                                 )

                write.table(data.frame(nCell = nCell,
                                       nFeat = nFeat,
                                       nAssay = nAssay,
                                       nCol = nCol,
                                       sizeTotal = object.size(qfeatures),
                                       sizeAssay = getAssaySize(qfeatures),
                                       sizeRowData = getRowDataSize(qfeatures),
                                       sizeColData = getColDataSize(qfeatures)),
                            file = file.path(destPath, "qfeatures_size_report.tsv"),
                            append = TRUE,
                            col.names = FALSE,
                            row.names = FALSE)
            }
        }
    }
}
