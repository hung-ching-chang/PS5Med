<img src="man/figure/hex-PS5.png" width="200">

## PS5Med: Partial Sum Statistic and Sample Splitting Strategy for MEDiation analysis 
![GitHub release](https://img.shields.io/badge/release-v1.0.0-blue.svg)

PS5Med is an R package designed for high-dimensional mediation analysis, featuring two methods: (`PS5` & `PS5.multisplit`). While `PS5` performs fast operation through a single sample split, `PS5.multisplit` enhances the robustness of mediation analysis by employing multiple sample split. The framework of our method is illustrated in the figure below.

<img src="man/figure/framework-PS5.png" width="4500">

## Installation

The package can be downloaded and installed from GitHub.

```
install.packages('devtools')
devtools::install_github('hung-ching-chang/PS%Med')
```
If the package cannot be installed successfully, it is needed to ensure that the dependency packages are installed. This package is based on R (>= 4.0.5), and the R code for installing the dependent packages is:
```
## from CRAN
CRAN.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
}
CRAN.packages(c("ncvreg", "Matrix", "mvtnorm", "parallel"))
```


## Code for the manuscript
The code for the simulation section in the PS5 manuscript can be found in the folder `inst/simulation`.


## Contact package maintainer:
Hung-Ching Chang

Email: rick.chang@pitt.edu

## Citation:
Chang, H.C., Fang, Y., Gorczyca, M., Batmanghelich, K., Tseng, G.C. (2024+). High-Dimensional Causal Mediation Analysis by Partial Sum Statistic and Sample Splitting Strategy in Imaging Genetics Application.


