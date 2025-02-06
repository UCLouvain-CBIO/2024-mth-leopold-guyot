#' Data Import for Quantitative Proteomics
#'
#' Imports quantitative proteomics data into a `QFeatures` object, compatible with
#' single-cell and multiplexed data types such as TMT, LFDIA, and plexDIA.
#'
#' @param quantitativePath A character string specifying the file path to the table
#'        containing quantitative data.
#' @param designPath A character string specifying the file path to the table with
#'        experimental design details. The design table must contain columns "runCol"
#'        and "quantCols" as required.
#' @param runCol A character string indicating the column in the design table that
#'        specifies the runs or batches.
#' @param dataDIA A character string specifying the file path to DIA data, only required
#'        if using "plexDIA" as the `type` parameter. Default is NULL.
#' @param type A character string specifying the type of quantitative data to import.
#'        Accepted values are "TMT", "LFDIA", or "plexDIA".
#'
#' @return A `QFeatures` object containing `SingleCellExperiment` assays based on the specified
#'         quantitative data type.
#'
#' @examples
#' # Importing TMT data
#' # qfeatures <- stepDataImport("quantData.csv", "designData.csv", "Run", type = "TMT")
#'
#' @export
stepDataImport <- function(
        quantitativePath,
        designPath,
        runCol,
        dataDIA = NULL,
        type) {
    quantitativeTable <- read.csv(quantitativePath)
    designTable <- read.csv(designPath)
    switch(type,
        "TMT" = .readQFeaturesTMT(quantitativeTable, designTable, runCol),
        "LFDIA" = .readQFeaturesLFDIA(quantitativeTable, designTable),
        "plexDIA" = .readQFeaturesPlexDIA(quantitativeTable, designTable, dataDIA),
        stop("Invalid type")
    )
}

#' Read Quantitative TMT Data into QFeatures
#'
#' Utility function to read TMT data and convert it into a `QFeatures` object.
#'
#' @param quantitativeTable A data frame containing quantitative data.
#' @param designTable A data frame containing experimental design details,
#'        which includes columns "runCol" and "quantCols".
#' @param runCol A character string specifying the name of the column in the design table
#'        that contains the runs or batches.
#'
#' @return A `QFeatures` object containing `SingleCellExperiment` assays.
#'
#' @importFrom scp readSCP
#'
#' @keywords internal
.readQFeaturesTMT <- function(quantitativeTable, designTable, runCol) {
    scp::readSCP(
        assayData = quantitativeTable,
        colData = designTable,
        runCol = runCol,
        verbose = FALSE
    )
}

#' Read Quantitative LFDIA Data into QFeatures
#'
#' Utility function to read LFDIA data and convert it into a `QFeatures` object.
#'
#' @param quantitativeTable A data frame containing quantitative data.
#' @param designTable A data frame containing experimental design details,
#'        which includes columns "runCol" and "quantCols".
#'
#' @return A `QFeatures` object containing `SingleCellExperiment` assays.
#'
#' @importFrom scp readSCPfromDIANN
#'
#' @keywords internal
.readQFeaturesLFDIA <- function(quantitativeTable, designTable) {
    scp::readSCPfromDIANN(
        assayData = quantitativeTable,
        colData = designTable,
        multiplexing = "none",
        verbose = FALSE
    )
}

#' Read Quantitative PlexDIA Data into QFeatures
#'
#' Utility function to read PlexDIA data and convert it into a `QFeatures` object.
#'
#' @param quantitativeTable A data frame containing quantitative data.
#' @param designTable A data frame containing experimental design details,
#'        which includes columns "runCol" and "quantCols".
#' @param dataDIA A character string specifying the file path to the DIA data.
#'
#' @return A `QFeatures` object containing `SingleCellExperiment` assays.
#'
#' @importFrom scp readSCPfromDIANN
#'
#' @keywords internal
.readQFeaturesPlexDIA <- function(quantitativeTable, designTable, dataDIA) {
    scp::readSCPfromDIANN(
        assayData = quantitativeTable,
        colData = designTable,
        multiplexing = "mTRAQ",
        extractedData = dataDIA,
        verbose = FALSE
    )
}
