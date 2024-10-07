################################################################################
############################# WORKFLOW PROFILING ###############################
################################################################################

################################# Parameters ###################################

nfeaturesRange <- c(500, 1000)

nRunsRange <- c(4, 8)

naRateRange <- c(0.5, 0.7)

nReplicates <- 1

############################### Sourcing #######################################

source(file.path("R", "workflowStepOneDataImport.R"))
source(file.path("R", "workflowStepTwoPreProcessing.R"))

source(file.path("R", "generateData.R"))
source(file.path("R", "profilingUtils.R"))

############################# Packages loading #################################

library("QFeatures")
library("scp")

################################## Main ########################################

main <- function() {
    # Ensure temp files are cleaned on exit
    on.exit({
        unlink("profilingTemp", recursive = TRUE)
        cat("Temporary files cleaned up.\n")
    })

    ############################# Data Generation ##############################

    ### TMT ###

    replicateData(
        nFeaturesRange = nfeaturesRange,
        nRunsRange = nRunsRange,
        naRateRange = naRateRange,
        nReplicates = nReplicates,
        type = "TMT",
        folder = "profilingTemp"
    )

    ### label free DIA-NN ###

    ### plexed DIA-NN ###

    ############################# Profiling ####################################
    profOutPath <- "dataOutput/profilingResults"

    if (dir.exists(profOutPath)) {
        unlink(profOutPath, recursive = TRUE)
    }
    profilingWrapper(
        dataFolder = "profilingTemp",
        profOutPath = profOutPath
    )
}

main()
