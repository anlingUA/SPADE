---
layout: default
title: SPADE tutorial
---

``` r
library(tidyverse)
#> ── Attaching packages ──────────────────────────────── tidyverse 1.3.2 ──
#> ✔ ggplot2 3.4.1     ✔ purrr   1.0.1
#> ✔ tibble  3.1.8     ✔ dplyr   1.1.0
#> ✔ tidyr   1.3.0     ✔ stringr 1.5.0
#> ✔ readr   2.1.4     ✔ forcats 1.0.0
#> ── Conflicts ─────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
library(SPADE)
```

Load the example data for simulated spatial gene expression, single cell
RNA-seq data, marker gene lists for cell types, spatial variable genes,
and spatial locations.

``` r
load('simulatedData.RData')
load('scRNAseq.RData')
load('MarkerGenes.RData')
location=read.csv('spa_loc.csv',header = T,row.names = 1)
```

Build cell type reference from scRNAseq

``` r
scref=scRefer(input_data=refData, 
              ct_var="celltype", 
              sample_var="orig.ident")
#> Loading required package: BiocGenerics
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:dplyr':
#> 
#>     combine, intersect, setdiff, union
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter,
#>     Find, get, grep, grepl, intersect, is.unsorted, lapply, Map,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
#>     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
#>     setdiff, sort, table, tapply, union, unique, unsplit,
#>     which.max, which.min
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages
#>     'citation("pkgname")'.
```

Identify cell types for each domain

``` r
marker_list=unlist(sign_list)
pseudo.spatial=simData[[1]]
nlay=max(location$domain)
CTperLayer=CTperDom(loc=location,stcount=pseudo.spatial,scref=scref,sign_list=marker_list,lasso=T)
#> Loading required package: lattice
#> 
#> Attaching package: 'caret'
#> The following object is masked from 'package:purrr':
#> 
#>     lift
#> Loading required package: Matrix
#> 
#> Attaching package: 'Matrix'
#> The following objects are masked from 'package:tidyr':
#> 
#>     expand, pack, unpack
#> Loaded glmnet 4.1-6
#> [1] TRUE
#> [1] TRUE
#> [1] TRUE
#> [1] TRUE
```

based on identified cell type, decompose each domain

``` r
CTest=SPADE(stcount=pseudo.spatial,scref=scref,
              sign_list=sign_list, # sc DEGs, it is in list format 
              loc=location,
              ctData=CTperLayer, # the estimated cell type per domain
              offset=10) # user input tuning parameter to adjust the threshold, if is NULL, then an average value will be used
```

merge domain information

``` r
truep=simData[[2]]
estCT=matrix(0,ncol = ncol(truep),nrow=nrow(truep))
rownames(estCT)=c(rownames(CTest[[1]]),rownames(CTest[[2]]),rownames(CTest[[3]]),rownames(CTest[[4]]))
colnames(estCT)=colnames(truep)
estCT[rownames(CTest[[1]]),colnames(CTest[[1]])] = CTest[[1]]
estCT[rownames(CTest[[2]]),colnames(CTest[[2]])] = CTest[[2]]
estCT[rownames(CTest[[3]]),colnames(CTest[[3]])] = CTest[[3]]
estCT[rownames(CTest[[4]]),colnames(CTest[[4]])] = CTest[[4]]
estCT=estCT[rownames(truep),colnames(truep)]
```

compare results with true proportion

``` r
RMSD <- round(sqrt(mean(as.matrix((truep - estCT)^2), na.rm = T)), digits = 5)
mAD <- round(mean(as.matrix(abs(truep - estCT)), na.rm = T), digits = 5)
pearson <- round(cor(c(as.matrix(truep)), c(as.matrix(estCT))), digits = 5)
print(c('RMSD'=RMSD,'mAD'=mAD,'pearson'=pearson))
#>    RMSD     mAD pearson 
#> 0.03931 0.01680 0.97866
```

visialize the results

``` r
colors=c("PG"="green4","Oligo"='lawngreen',"ET"="gold","mitral"="purple","tufted"="brown2","Astro"="magenta",'Micro'="lightskyblue","GC"="orange",'Immature'="lightcyan3",'OPC'='royalblue1')

dot_plot(truep,estCT,10,colors=colors)
```

![](Intro_to_SPADE_files/figure-markdown/unnamed-chunk-8-1.png)

``` r
library(scatterpie)
fuldat=cbind(location,estCT)
scatter_pie(fuldat=fuldat,dat=estCT,fillcolors=colors)
```

![](Intro_to_SPADE_files/figure-markdown/unnamed-chunk-8-2.png)

``` r

featureplot(estCT,fuldat,c('GC','mitral'))
```

![](Intro_to_SPADE_files/figure-markdown/unnamed-chunk-8-3.png)
