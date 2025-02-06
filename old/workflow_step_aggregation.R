#' Aggregate Features Over Assays
#'
#' Aggregates features over specified assays in a `QFeatures` object.
#' It calculates the median of the values for each group across the
#'  specified assays.
#'
#' @param qfeatures A `QFeatures` object containing the assay data to be aggregated.
#' @param assays A integer vector specifying the indices of the assays to be
#' aggregated within the `QFeatures` object.
#' @param type A character string specifying the type of aggregation to be
#' performed c("PSM", "peptide").
#'
#' @return A modified `QFeatures` object with new assays added. The names of these
#'         new assays are prefixed according to the type of aggregation used.
stepAggregate <- function(qfeatures, assays, type) {
    suppressMessages(
        switch(type,
            "PSM" = .aggregatePSM(qfeatures, assays),
            "peptide" = .aggregatePeptide(qfeatures, assays),
            stop("Invalid aggregation type")
        )
    )
}


#' Aggregate PSM to peptides
#'
#' @param qfeatures A `QFeatures` object containing the assay data to be
#' aggregated.
#' @param assays A integer vector specifying the indices of the assays to be
#' aggregated within the `QFeatures` object.
#'
#' @return A modified `QFeatures` object with new assays added.
#' The names of these new assays are prefixed with "peptides_".
#'
#' @importFrom matrixStats colMedians
#' @importFrom scp aggregateFeaturesOverAssays
#'
#' @keywords internal
.aggregatePSM <- function(qfeatures, assays) {
    aggregateFeaturesOverAssays(qfeatures,
        i = assays,
        fcol = "peptidesId",
        name = paste0(
            "peptides_",
            names(qfeatures[assays])
        ),
        fun = matrixStats::colMedians, na.rm = TRUE
    )
}

.aggregatePeptide <- function(qfeatures, assays) {
    stop("not yet implemented")
}
