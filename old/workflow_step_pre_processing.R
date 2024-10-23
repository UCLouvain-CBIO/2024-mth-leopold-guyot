#' Preprocess Quantitative Proteomic Data
#'
#' Preprocesses `QFeatures` data according to the specified type, with support
#' for TMT, LFDIA, and plexDIA data types. Each data type follows a different
#' preprocessing pipeline that includes filtering and normalizing the data.
#'
#' @param qfeatures A `QFeatures` object containing the assays to be
#'  preprocessed.
#' @param assays A character vector specifying the names of the assays to be
#'  processed.
#' @param type A character string specifying the type of proteomic data.
#'  Accepted values are "TMT", "LFDIA", or "plexDIA".
#'
#' @return The `QFeatures` object after preprocessing.
#'
#' @examples
#' # Assuming qfeatures is a valid QFeatures object and assays is a vector of
#' # assay indices:
#' # qfeatures <- stepPreProcessing(qfeatures, assays = c("assay1", "assay2"), type = "TMT")
#'
#' @export
stepPreProcessing <- function(qfeatures, assays, type) {
    switch(type,
        "TMT" = .stepPreProcessingTMT(qfeatures, assays),
        "LFDIA" = .stepPreProcessingLFDIA(qfeatures, assays),
        "plexDIA" = .stepPreProcessingPlexDIA(qfeatures, assays),
        stop("Invalid preprocessing type")
    )
}

#' TMT Data Preprocessing Pipeline
#'
#' Preprocesses TMT-type proteomic data by applying several filtering and normalization steps.
#'
#' @param qfeatures A `QFeatures` object containing the assays to be preprocessed.
#' @param assays A character vector specifying the names of the assays to be processed.
#'
#' @return The `QFeatures` object after preprocessing for TMT data.
#'
#' @importFrom QFeatures zeroIsNA
#'
#' @keywords internal
.stepPreProcessingTMT <- function(qfeatures, assays) {
    qfeatures <- zeroIsNA(qfeatures, i = assays)

    qfeatures <- .filterPSMTMT(qfeatures, assays)
    qfeatures <- .filterAssayTMT(qfeatures)
    qfeatures <- .filterSCRTMT(qfeatures, assays)
    qfeatures <- .filterFDRTMT(qfeatures, assays)
}

.setpPreProcessingLFDIA <- function(qfeatures, assays) {
    stop("not yet implemented")
}

.stepPreProcessingPlexDIA <- function(qfeatures, assays) {
    stop("not yet implemented")
}

#' Filter TMT PSMs
#'
#' Filters PSMs (Peptide-Spectrum Matches) for TMT data based on specific criteria such as reverse sequences,
#' contaminants, and PIF (Precursor Ion Fraction) values.
#'
#' @param qfeatures A `QFeatures` object containing the assays to be filtered.
#' @param assays A character vector specifying the names of the assays to be filtered.
#'
#' @return The `QFeatures` object with filtered features for TMT data.
#'
#' @importFrom QFeatures filterFeatures
#'
#' @keywords internal
.filterPSMTMT <- function(qfeatures, assays) {
    suppressMessages(filterFeatures(
        qfeatures,
        ~ Reverse != "+" &
            Potential.contaminant != "+" &
            !is.na(PIF) & PIF > 0.8,
        i = assays
    ))
}

#' Filter TMT Assays Based on Row Dimensions
#'
#' Filters TMT assays based on the number of features. Only assays with more than
#' 120 features are retained.
#'
#' @param qfeatures A `QFeatures` object containing the assays to be filtered.
#'
#' @return The `QFeatures` object with only the assays that meet the row dimension criteria.
#'
#' @keywords internal
.filterAssayTMT <- function(qfeatures) {
    keepAssay <- dims(qfeatures)[1, ] > 120
    qfeatures[, , keepAssay]
}

#' Filter TMT Samples Based on Sample Carrier Ratio (SCR)
#'
#' Computes and filters based on the Sample Carrier Ratio (SCR) for TMT data.
#' Samples are filtered out if their SCR exceeds 0.1.
#'
#' @param qfeatures A `QFeatures` object containing the assays to be filtered.
#' @param assays A character vector specifying the names of the assays to be filtered.
#'
#' @return The `QFeatures` object after filtering based on the SCR.
#'
#' @importFrom QFeatures filterFeatures
#' @importFrom BiocGenerics dims
#' @importFrom scp computeSCR
#'
#' @keywords internal
.filterSCRTMT <- function(qfeatures, assays) {
    qfeatures <- computeSCR(qfeatures,
        i = assays,
        colvar = "SampleType",
        carrierPattern = "Carrier",
        samplePattern = "Blank|Macrophage|Monocyte",
        sampleFUN = "mean",
        rowDataName = "MeanSCR"
    )
    qfeatures <- suppressMessages(filterFeatures(qfeatures,
        ~ !is.na(MeanSCR) &
            MeanSCR < 0.1,
        i = assays
    ))
}

#' Filter TMT Data Based on False Discovery Rate (FDR)
#'
#' Filters TMT data based on q-values for proteins, using a 1% FDR threshold.
#'
#' @param qfeatures A `QFeatures` object containing the assays to be filtered.
#' @param assays A character vector specifying the names of the assays to be filtered.
#'
#' @return The `QFeatures` object with filtered features based on FDR.
#'
#' @importFrom QFeatures filterFeatures
#' @importFrom scp pep2qvalue
#' @keywords internal
.filterFDRTMT <- function(qfeatures, assays) {
    qfeatures <- pep2qvalue(qfeatures,
        i = assays,
        PEP = "PEP",
        groupBy = "Leading.razor.protein",
        rowDataName = "qvalueProteins"
    )
    suppressMessages(filterFeatures(
        qfeatures,
        ~ qvalueProteins < 0.01
    ))
}

#' Normalize TMT Data by Dividing by Reference Samples
#'
#' Normalizes TMT data by dividing each sample by the corresponding reference
#' sample for normalization.
#'
#' @param qfeatures A `QFeatures` object containing the assays to be normalized.
#' @param assays A character vector specifying the names of the assays to be normalized.
#'
#' @return The `QFeatures` object with normalized assays based on reference samples.
#'
#' @importFrom scp divideByReference
#'
#' @keywords internal
.divideReferenceTMT <- function(qfeatures, assays) {
    qfeatures <- divideByReference(qfeatures,
        i = assays,
        colvar = "SampleType",
        samplePattern = ".",
        refPattern = "Reference"
    )
}
