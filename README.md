# SPADE
## Spatial Deconvolution for Domain Specific Cell-type Estimation
![](./SPADEdiagram.png)
SPADE, a reference-based approach, harnesses the power of single-cell RNA sequencing data, spatial location details, and histological information to accurately estimate the proportions of various cell types at each spatial location. A crucial characteristic of SPADE is its ability to account for cell type sparsity across locations, thereby enabling the identification of cell types prior to proportion estimation.

## Installation

### Dependencies
- R version >= 4.0.0.
- R packages: Biobase, caret, glmnet, tidyverse, NOISeq, EBImage

```{r}
# install devtools if necessary
install.packages('devtools')

# install the SPADE package
devtools::install_github('anlingUA/SPADE')

# load
library(SPADE)
```
- The first step of SPADE utilizes spaGCN to cluster locations into domains. To get the identified spatial domain of your data, refer to the tutorial of [spaGCN](https://github.com/jianhuupenn/SpaGCN)

- For an example how to use [SPADE](https://anlingua.github.io/SPADE/Intro_to_SPADE.html)
- Data used for the tutorial can be downloaded from [here](https://figshare.com/projects/SPADE/168116)
