# Benchmark of QFeatures perfromances

## Installation

```r
if (!require("remotes", quietly = TRUE)){
    install.packages("remotes")
}
# Install the package
remotes::install_github("UCLouvain-CBIO/2024-mth-leopold-guyot",
    build_manual = TRUE
)
# Load the package
library(benchmarkQFeatures)
```

## Launch Profiling

### In R
```r
# Load the package
library(benchmarkQFeatures)

# Run the profiling

profilingExecution(nFeaturesRange = c(500, 1000, 1500),
                        nRunsRange = c(4, 8, 16),
                        naRateRange = c(0.4, 0.6, 0.8),
                        nReplicates = 3,
                        outputFolder = "dataOutput/profilingResults")
```

### In command line

First install the script within R

```r
# Load the package
library(benchmarkQFeatures)

# Install scripts
installScripts("~/dev/scripts")
```

Launch script execution trough terminal

```console
~$ # without parameters
~$ ./workflow_profiling

~$ # with parameters
~$ ./workflow_profiling "list(naRateRange = c(0.4, 0.6), nFeaturesRange = c(500, 1000), nRunsRange = c(4, 8), nReplicates = 3, outputFolder = "profilingOutput")"
```
