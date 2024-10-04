### Step for data import of quantitative proteomic data ###

### Packages ###

library("scp")

### Main function ###

# Output: returns a QFeatures object
#    (scp -> with SingleCellExperiment assays)
#
# Params:
# - `quantitativeTable` (`character`) the path to the table with quantitative data
# - `designTable` (`character`) the path to the table with experimental design.
#   This table should contain a "runCol" and a "quantCols" column
# - `runCol` (`character`) the name of the column in the design table that
#   contains the runs/batches
# - `dataDIA` (`character`) the path to the DIA data (only for plexed DIA-NN)
# - `type` (`character`) the type of quantitative data to import
#   possible values => c("TMT", "LFDIA", "plexDIA")

stepDataImport <- function(quantitativeTable,
                           designTable,
                           runCol,
                           dataDIA = NULL,
                           type) {
    switch(type,
           "TMT" = .readQFeaturesTMT(quantitativeTable, designTable, runCol),
           "LFDIA" = .readQFeaturesLFDIA(quantitativeTable, designTable),
           "plexDIA" = .readQFeaturesPlexDIA(quantitativeTable, designTable, dataDIA),
           stop("Invalid type")
    )

}

### Utility functions ###

# Output: a `QFeatures` object with singleCellExperiment assays
# Params:
# - `quantitativeTable` (`character`) the path to the table with quantitative data
# - `designTable` (`character`) the path to the table with experimental design.
#   This table should contain a "runCol" and a "quantCols" column
# - `runCol` (`character`) the name of the column in the design table that
#   contains the runs/batches

.readQFeaturesTMT <- function(quantitativeTable, designTable, runCol) {
    scp::readSCP(assayData = quantitativeTable,
                 colData = designTable,
                 runCol = runCol,
                 verobse = FALSE)
}

# Output: a `QFeatures` object with singleCellExperiment assays
# Input:
# - `quantitativeTable` (`character`) the path to the table with quantitative data
# - `designTable` (`character`) the path to the table with experimental design.
#   This table should contain a "runCol" and a "quantCols" column

.readQFeaturesLFDIA <- function(quantitativeTable, designTable) {
    scp::readSCPfromDIANN(assayData = quantitativeTable,
                          colData = designTable,
                          multiplexing = "none",
                          verobse = FALSE)
}

# Output: a `QFeatures` object with singleCellExperiment assays
# Input:
# - `quantitativeTable` (`character`) the path to the table with quantitative data
# - `designTable` (`character`) the path to the table with experimental design.
#   This table should contain a "runCol" and a "quantCols" column
# - `dataDIA` (`character`) the path to the DIA data (only for plexed DIA-NN)
.readQFeaturesPlexDIA <- function(quantitativeTable, designTable, dataDIA) {
    scp::readSCPfromDIANN(assayData = quantitativeTable,
                          colData = designTable,
                          multiplexing = "mTRAQ",
                          extractedData = dataDIA,
                          verobse = FALSE)
}
