# Benchmark of QFeatures perfromances

## Vignette Benchmarking

### In command line

Run benchmark + generate result report (in reports/leduc2022Results.html)
```console
~$ chmod +x scripts/vignetteBenchmarkScript.R

~$ # without parameters
~$ scripts/vignetteBenchmarkScript.R 

~$ # with parameters
~$ scripts/vignetteBenchmarkScript.R -cellRange "c(1000, 2000, 4000)" -nReplicates 3
```

### In R

```R
source(file.path("R", "vignette_benchmark.R"))
source(file.path("R", "utils.R"))

cellRange <- c(1000, 2000, 4000)
nReplicates <- 3

leduc2022Benchmark(cellRange, nReplicates)

renderRmarkdown("leduc2022Results.rmd")
```

## 4 Variables Benchmark

``` R
source(file.path("R", "4_variables_benchmark.R"))

source(file.path("R", "utils.R"))

renderRmarkdown("4_var_benchmark_results.rmd")
```

## Individual Steps Benchmark

``` R
source(file.path("R", "individual_step_benchmark.R"))

source(file.path("R", "utils.R"))

renderRmarkdown("individual_step_results.rmd")
```

## HeLa Benchmark

``` R
source(file.path("R", "hela_to_qfeatures.R"))
source(file.path("R", "hela_to_qfeatures_opti.R"))
source(file.path("R", "hela_benchmark.R"))

source(file.path("R", "utils.R"))

# renderRmarkdown("")
```

## Summary Figures

``` console
~$ Rscript R/summary_figures.R
```
