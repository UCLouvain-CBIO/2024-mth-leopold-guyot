library(peakRAM)
library(QFeatures)
library(scp)
library(benchmarkme)

source(file.path("R", "utils.R"))
helaBenchmark <- function(nreplicates,
                               outputDir = "dataOutput/helaBenchmark") {
    if (file.exists(file.path(outputDir, "size_report.tsv"))) {
        unlink(file.path(outputDir, "size_report.tsv"))
        cat(file.path(outputDir, "size_report.tsv"), " has been deleted..\n")
    }
    if (file.exists(file.path(outputDir, "memoryOutput"))) {
        unlink(file.path(outputDir, "memoryOutput"), recursive = TRUE)
        cat(file.path(outputDir, "memoryOutput"), " has been deleted\n")
    }
    dir.create(file.path(outputDir, "memoryOutput"))
    write.table(data.frame(rep = integer(),
                           state = character(),
                           sizeTotal = numeric(),
                           sizeAssay = numeric(),
                           sizeRowData = numeric(),
                           sizeColData = numeric()),
                file = file.path(outputDir, "size_report.tsv"),
                append = FALSE,
                col.names = TRUE,
                row.names = FALSE)

    for (rep in 1:nreplicates){
        hela <- readRDS("dataOutput/hela/qfeatures_PSM.rds")
        cat("Starting replicate: ", rep)
        write.table(
            data.frame(
                rep = rep,
                state = "before",
                sizeTotal = object.size(hela),
                sizeAssay = getAssaySize(hela),
                sizeRowData = getRowDataSize(hela),
                sizeColData = getColDataSize(hela)
            ),
            file = file.path(outputDir, "size_report.tsv"),
            append = TRUE,
            col.names = FALSE,
            row.names = FALSE
        )
        res <- peakRAM(
            assaysNames <- names(hela),
            hela <- filterFeat(hela),
            hela <- filterSamples(hela),
            hela <- zeroIsNA(hela, assaysNames),
            hela <- aggregationToPep(hela, assaysNames),
            hela <- joining(hela, assaysNames),
            hela <- normalisation(hela)
            #hela <- aggregationToPro(hela)
        )

        write.csv(res, file = file.path(outputDir, "memoryOutput",
                                        paste0(paste(
                                            rep, sep = "_"),".csv")))
        write.table(
            data.frame(
                rep = rep,
                state = "after",
                sizeTotal = object.size(hela),
                sizeAssay = getAssaySize(hela),
                sizeRowData = getRowDataSize(hela),
                sizeColData = getColDataSize(hela)
            ),
            file = file.path(outputDir, "size_report.tsv"),
            append = TRUE,
            col.names = FALSE,
            row.names = FALSE
        )
    }


    sd <- benchmarkme::get_sys_details()
    if (file.exists(file.path(outputDir, "hardware_software.txt"))) {
        unlink(file.path(outputDir, "hardware_software.txt"))
    }
    sink(file.path(outputDir, "hardware_software.txt"))
    cat("Machine: ", sd$sys_info$sysname, " (", sd$sys_info$release, ")\n",
        "R version: R.", sd$r_version$major, ".", sd$r_version$minor,
        " (svn: ", sd$r_version$`svn rev`, ")\n",
        "RAM: ", round(sd$ram / 1E9, 1), " GB\n",
        "CPU: ", sd$cpu$no_of_cores, " core(s) - ", sd$cpu$model_name, "\n",
        sep = "")
    cat(" ----------------------------------------------------------------- \n")
    print(sessionInfo())
    sink()
}

filterFeat <- function(qfeatures){
    filterFeatures(qfeatures,
                   ~ Potential.contaminant != "+" &
                       Reverse != "+" &
                       PEP <= 0.01)
}


filterSamples <- function(qfeatures) {
    subsetByColData(qfeatures, !grepl("2013_06_", qfeatures$runCol))
}



aggregationToPep <- function(qfeatures, assayName) {
    aggregateFeaturesOverAssays(qfeatures, assayName,
        fcol = "Modified.sequence", name = paste0(assayName, "_peptides"),
        fun = colMedians)
}

joining <- function(qfeatures, assayName) {
    joinAssays(qfeatures, paste0(assayName, "_peptides"), name = "peptides")
}

normalisation <- function(qfeatures) {
    QFeatures::normalize(qfeatures,
                            i = "peptides",
                            method = "div.median",
                            name = "peptides_norm"
    )
}

aggregationToPro <- function(qfeatures) {
    print(names(rowData(qfeatures)[["peptides_norm"]]))
    aggregateFeatures(qfeatures, i = "peptides_norm",
                      fcol = "Leading.razor.protein",
                      name = "proteins", fun = colMedians)
}


# MAIN

helaBenchmark(3)
