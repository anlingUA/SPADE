# SPADE
## Spatial Deconvolution for Domain Specific Cell-type Estimation
![](./SPADEdiagram.png)
SPADE, a reference-based approach, harnesses the power of single-cell RNA sequencing data, spatial location details, and histological information to accurately estimate the proportions of various cell types at each spatial location. A crucial characteristic of SPADE is its ability to account for cell type sparsity across locations, thereby enabling the identification of cell types prior to proportion estimation.

## Installation
```{r}
# install devtools if necessary
install.packages('devtools')

# install the SPADE package
devtools::install_github('anlingUA/SPADE')

# load
library(SPADE)
```

For an example how to use [SPADE](https://yylu5.github.io/SPADE/)
