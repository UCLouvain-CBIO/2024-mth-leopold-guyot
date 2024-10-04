################################################################################
############################# WORKFLOW PROFILING ###############################
################################################################################

################################# Parameters ###################################

nfeaturesRange <- c(500, 1000)

nRunsRange <- c(4, 16)

naRateRange <- c( 0.5, 0.7)

nReplicates <- 1

############################### Sourcing #######################################

source("workflowSriptsR/01_DataImport.R")
source("utilsR/generateTMT.R")
source("utilsR/replicateData.R")
source("utilsR/profilingWrapper.R")

############################# Packages loading #################################

library("profvis")

############################# Data Generation ##################################

### TMT ###

replicateData(nFeaturesRange = nfeaturesRange,
             nRunsRange = nRunsRange,
             naRateRange = naRateRange,
             nReplicates = nReplicates,
             type = "TMT",
             folder = "profTemp")

### label free DIA-NN ###

### plexed DIA-NN ###

############################# Profiling ########################################

profilingWrapper(dataFolder = "profTemp", profOutPath = "profOutPath")

########################## Clean Out Temp Files ################################

unlink("profTemp", recursive = TRUE)
