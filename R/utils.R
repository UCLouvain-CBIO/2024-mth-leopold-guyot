library(rmarkdown)

renderRmarkdown <- function(rmdFile){
    rmarkdown::render(file.path("Rmd", rmdFile),
                      output_file = sub(".rmd", ".html", rmdFile),
                      output_dir = "reports",
                      knit_root_dir = getwd())
}
