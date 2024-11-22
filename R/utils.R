library(rmarkdown)

renderRmarkdown <- function(rmdFile){
    rmarkdown::render(file.path("Rmd", rmdFile),
                      output_file = sub(".rmd", ".html", rmdFile),
                      output_dir = "reports",
                      knit_root_dir = getwd())
}

getColDataSize <- function(qfeatures) {
    object.size(qfeatures@colData)
}

getRowDataSize <- function(qfeatures) {
    sizes <- sapply(names(qfeatures), function(x) {
        object.size(rowData(qfeatures[[x]]))
    })

    sum(sizes)
}

getAssaySize <- function(qfeatures) {
    sizes <- sapply(names(qfeatures), function(x) {
        object.size(assay(qfeatures[[x]]))
    })

    sum(sizes)
}
