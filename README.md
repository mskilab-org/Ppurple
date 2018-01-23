Ppurple
=======
Probabilistic purity ploidy estimation

Installation
------------

1. Install dependent packages and latest Bioconductor (if you haven't already)

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
```

2. Install devtools from CRAN (if you don't have it already)

```{r}
install.packages('devtools')
```

3. Install mskilab R dependencies (gUtils and bamUtils)


```{r}
devtools::install_github('mskilab/gUtils)
devtools::install_github('mskilab/bamUtils)
```

4. Install Ppurple


```{r}
devtools::install_github('mskilab/Ppurple')
```


