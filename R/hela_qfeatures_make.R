library(curl)
url <- "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2023/12/PXD042233"
files <- c(
    "precursors_aggregated.zip",
    "peptides_aggregated.zip",
    "geneGroups_aggregated.zip"
)

for (file in files) {
    curl_download(url = file.path(url, file), destfile = file)
    unzip(file, exdir = sub("\\..[^\\.]*$", "", file))
    file.remove(file)
}

# Load precursor quantification data
PSM <-  read.csv("precursors_aggregated/precursors/intensities_wide_selected_N49339_M07444.csv")

peptides <- read.csv("peptides_aggregated/peptides/intensities_wide_selected_N42723_M07444.csv")

geneGroups <- read.csv("geneGroups_aggregated/geneGroups/intensities_wide_selected_N04547_M07444.csv")
