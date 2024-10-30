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
