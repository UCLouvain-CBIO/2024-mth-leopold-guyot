### Replicate Data Function ###
replicateData <- function(nFeaturesRange, nRunsRange, naRateRange, nReplicates, type, folder) {
    if (!dir.exists(folder)) {
        dir.create(folder)
    }
    generateFunction <- switch(type,
                               "TMT" = generateTMTtoFile,
                               "LFDIA" = generateLFDIAtoFile,
                               "plexDIA" = generatePlexDIAtoFile,
                               stop("Unknown type in replicateData: ", type)
    )
    for (nFeatures in nFeaturesRange) {
        for (nRuns in nRunsRange) {
            for (naRate in naRateRange) {
                for (nReplicate in 1:nReplicates) {
                    .createReplicate(nFeatures = nFeatures,
                                    nRuns = nRuns,
                                    naRate = naRate,
                                    nReplicate = nReplicate,
                                    type = type,
                                    folder = folder,
                                    generateFunction = generateFunction)
                }

            }
        }
    }
}

.createReplicate <- function(nFeatures,
                             nRuns,
                             naRate,
                             nReplicate,
                             type,
                             folder,
                             generateFunction) {
    # Define the subfolder name based on the current combination
    subfolder <- paste0(folder, "/",
                        type,
                        "_nFeatures_", nFeatures,
                        "_nRuns_", nRuns,
                        "_naRate_", naRate,
                        "_replicate_", nReplicate)

    if (!dir.exists(subfolder)) {
        dir.create(subfolder)
    }

    # Generate the file paths for the files
    quantPath <- paste0(subfolder, "/quant", type, ".csv")
    designPath <- paste0(subfolder, "/design", type, ".csv")

    generateFunction(nRuns = nRuns,
                      nFeatures = nFeatures,
                      naRate = naRate,
                      quantPath = quantPath,
                      designPath = designPath)

    cat("Generated", type, "data for nFeatures =", nFeatures,
        "nRuns =", nRuns, "naRate =", naRate, "in", subfolder, "\n")
}
