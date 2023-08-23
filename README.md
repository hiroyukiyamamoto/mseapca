# mseapca
Metabolite set enrichment analysis for loadings

**mseapca**  provides functions for metabolite set enrichment analysis with loading in principal component anaysis and partial least squares, and single sample enrichment analysis.

- Metabolite set enrichment analysis [1] can be performed using the "msea_ora" or "msea_sub" functions.

- Our mseapca package utilizes the loadings package ( https://cran.r-project.org/web/packages/loadings ) to conduct multivariate analyses like PCA and PLS, along with their respective loadings.

- The "ssea_ora" function can be used to compute single sample enrichment analysis based on over-representation analysis [2].

- Our mseapca package can incorporate a metabolite set list from the PathBank database by referencing the AHPathbankDbs Bioconductor package or through user-generated lists.

**References**  
[1] Yamamoto H. et al., BMC Bioinformatics, (2014) 15(1):51. doi: https://doi.org/10.1186/1471-2105-15-51  
[2] Yamamoto H. , Jxiv, (2023). doi: https://doi.org/10.51094/jxiv.484  

## Installation (in preparation)

The latest stable version can be installed from CRAN:

``` r
install.packages("mseapca")
```

The latest development version can be installed from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("hiroyukiyamamoto/mseapca")
```
