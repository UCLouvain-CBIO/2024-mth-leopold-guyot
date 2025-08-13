# Relevant PRs:

## subsetByColData

https://github.com/waldronlab/MultiAssayExperiment/pull/334

## aggregateFeatures

https://github.com/UCLouvain-CBIO/scp/pull/80
and
https://github.com/rformassspectrometry/QFeatures/pull/243
and
https://github.com/rformassspectrometry/QFeatures/pull/233

## SE instead of SE

https://github.com/UCLouvain-CBIO/scp/issues/83
and
https://github.com/rformassspectrometry/QFeatures/pull/229
and
https://github.com/UCLouvain-CBIO/scp/pull/85
and
https://github.com/UCLouvain-CBIO/scp/pull/90


# Benchmark of QFeatures perfromances

## Real pipeline benchmark

### In command line

Run benchmark + generate result report (in reports/leduc2022Results.html)
```console
~$ chmod +x scripts/vignetteBenchmarkScript.R

~$ # without parameters
~$ scripts/vignetteBenchmarkScript.R 

~$ # with parameters
~$ scripts/vignetteBenchmarkScript.R -cellRange "c(1000, 2000, 4000, 8000, 16000)" -nReplicates 3
```

### In R

```R
source(file.path("R", "vignette_benchmark.R"))
source(file.path("R", "utils.R"))

cellRange <- c(1000, 2000, 4000, 8000, 16000)
nReplicates <- 3

leduc2022Benchmark(cellRange, nReplicates)

renderRmarkdown("leduc2022Results.rmd")
```

## QFeatures object memory usage benchmark

``` R
source(file.path("R", "4_variables_benchmark.R"))

source(file.path("R", "utils.R"))

renderRmarkdown("4_var_benchmark_results.rmd")
```

## Individual Components Benchmark

``` R
source(file.path("R", "individual_step_benchmark.R"))

source(file.path("R", "utils.R"))

renderRmarkdown("individual_step_results.rmd")
```

## Challenging dataset Benchmark

``` R
source(file.path("R", "hela_to_qfeatures.R"))
source(file.path("R", "hela_to_qfeatures_opti.R"))
source(file.path("R", "hela_benchmark.R"))

source(file.path("R", "hela_benchmark_analysis.R"))
```
## SCE vs SE

``` R
source(file.path("R", "vignette_benchmark_SCE_SE.R"))
```
## Optimisation summary benchmark

``` R
source(file.path("R", "setupVerR.R"))
source(file.path("R", "analysisVerR.R"))
```

# Pseudobulking simulation

## Quality control

``` R
source(file.path("R", "artificialSimQC.R"))
```

## Simulation + model

``` R
source(file.path("R", "artificialSimulationPB.R"))
source(file.path("R", "artificialSimulationPBResults.R"))
source(file.path("R", "artificialSimulationPBRuntime.R"))
```
# Report Figures

``` R
source(file.path("R", "summary_figures.R"))
```
