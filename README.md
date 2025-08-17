# mseapca

[![CRAN status](https://www.r-pkg.org/badges/version/mseapca)](https://cran.r-project.org/package=mseapca)
[![Support ? Suzuri Shop](https://img.shields.io/badge/Support-Suzuri%20Shop-blue)](https://suzuri.jp/EigenStudio/18422952/heavyweight-t-shirt/l/white)

Metabolite set enrichment analysis for loadings

**mseapca**  provides functions for metabolite set enrichment analysis with loading in principal component anaysis and partial least squares, and single sample enrichment analysis.

- Metabolite set enrichment analysis [1] can be performed using the "msea_ora" or "msea_sub" functions.

- Our mseapca package utilizes the loadings package ( https://cran.r-project.org/web/packages/loadings ) to conduct multivariate analyses like PCA and PLS, along with their respective loadings.

- The "ssea_ora" function can be used to compute single sample enrichment analysis based on over-representation analysis [2].

- The "msea_ora" function performs over-representation analysis while accounting for undetected metabolites [3].

- The "msea_ora_range" function estimates the possible range of p-values under uncertainty caused by undetected metabolites [3].

**References**  
[1] Yamamoto H. et al., BMC Bioinformatics, (2014) 15(1):51. doi: https://doi.org/10.1186/1471-2105-15-51  
[2] Yamamoto H. , Jxiv, (2023). doi: https://doi.org/10.51094/jxiv.484  
[3] Yamamoto H. , Jxiv, (2024). doi: https://doi.org/10.51094/jxiv.954

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

## Note on ora_est() function

The current ora_est() function in the main branch contains a known bug when calculating enrichment under certain conditions (e.g., when accounting for undetected metabolites).

A fixed version is available in the development directory on GitHub:
dev/ora_est.R
We recommend using this version until the fix is included in the next CRAN release.

