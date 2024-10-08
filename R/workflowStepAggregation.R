stepAggregate <- function(qfeatures, assays){
    suppressMessages(
        aggregateFeaturesOverAssays(qfeatures,
                                i = assays,
                                fcol = "peptidesId",
                                name = paste0("peptides_",
                                              names(qfeatures[assays])),
                                fun = matrixStats::colMedians, na.rm = TRUE)
    )
}
