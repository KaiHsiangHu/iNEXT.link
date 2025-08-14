---
output:
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include=FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "",
  fig.path = "README/README-",
  fig.align = 'center',
  fig.retina=2, fig.width = 7, fig.height = 5,
  warning = FALSE,
  message = FALSE
)
options("width"=200)
```

iNEXT.link (R package)
=====

<h5 align="right">Latest version: `r Sys.Date()`</h5>

<font color="394CAE"><h3 color= 394CAE style = "font-weight: bold"> Introduction to iNEXT.link (R package): Excerpt from iNEXT.link UserGuide </h3> </font>
<br>
<h5><b>Anne Chao, Kai-Hsiang Hu, K.W. Chen, C.G. Lo, S.Y. Wang</b>
<br><br>
<i>Institute of Statistics, National Tsing Hua University, Hsin-Chu, Taiwan 30043</i>  
</h5>
<br>
`iNEXT.link` (INterpolation and EXTrapolation in Network diversity) is an R package, available in [Github](https://github.com/AnneChao). Here we provide a quick introduction demonstrating how to run iNEXT.link. An online version of [iNEXT.link Online](https://chao.shinyapps.io/iNEXT_link/) is also available for users without an R background. Detailed information about all functions in iNEXT.link is provided in the iNEXT.link Manual in [iNEXT.link_vignettes](http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/A%20Quick%20Introduction%20to%20iNEXT.link%20via%20Examples.html), which is also available from [Anne Chao's website](http://chao.stat.nthu.edu.tw/wordpress/software_download/).

`iNEXT.link` is an R package that extends the concepts of `iNEXT.3D` (Chao et al., 2021), `iNEXT.4step` (Chao et al., 2020) and `iNEXT.beta3D` (Chao et al., 2023) to ecological networks (Chiu et al., 2023). It is primarily designed to calculate and analyze various measures of diversity in ecological networks. Specifically, the package calculates three Hill numbers of order q (species richness, Shannon diversity, and Simpson diversity) in taxonomic diversity level, as well as phylogenetic and functional diversity levels.

For single ecological networks, `iNEXT.link` provides tools for analyzing diversity. The package provides two types of rarefaction and extrapolation (R/E) sampling curves to estimate diversity and confidence intervals for single ecological network. These include sample-size-based (or size-based) R/E curves and coverage-based R/E curves.

Moreover, `iNEXT.link` offers dissimilarity-turnover curves based on coverage-based R/E for gamma, alpha, and beta diversity measures, which can be used to compare diversity patterns across different ecological networks.


### SOFTWARE NEEDED TO RUN INEXT.3D IN R

-   Required: [R](http://cran.rstudio.com/)
-   Suggested: [RStudio IDE](http://www.rstudio.com/ide/download/)

### HOW TO RUN iNEXT.link:

The `iNEXT.link` package can be downloaded from Anne Chao's [iNEXT.link_github](https://github.com/AnneChao/iNEXT.link) using the following commands. For a first-time installation, additional visualization extension packages (`ggplot2`) from CRAN and (`iNEXT.3D`), (`iNEXT.4steps`), and (`iNEXT.beta3D`) from Anne Chao's github must be installed and loaded. 

```{r eval=FALSE}
## install iNEXT.link package from CRAN
# install.packages("iNEXT.link")  # coming soon

## install the latest version from github
install.packages('devtools')
library(devtools)

# install_github('AnneChao/iNEXT.3D')
# install_github('AnneChao/iNEXT.4steps')
# install_github('AnneChao/iNEXT.beta3D')

install_github('AnneChao/iNEXT.link')

## import packages
library(iNEXT.link)
```



In this document, we provide a quick introduction demonstrating how to run the package `iNEXT.link`(iNterpolation and EXTrapolation in Network diversity). `iNEXT.link` has several main functions: 

## Functions for single network:

- **iNEXT.link** : Computes rarefaction/extrapolation taxonomic, phylogenetic, and functional diversity estimates and sample coverage estimates.

- **DataInfo.link** : exhibits basic data information
    
- **estimateD.link** : computes species diversity with a particular user-specified level of sample size or sample coverage.
    
- **ObsAsy.link**: compute asymptotic or empirical(observed) diversity of order q.

- **Completeness.link** : Calculates estimated sample completeness with order q. 

- **Spec.link.est** :  Computes standardized specialization estimation under specified sample coverage with order q.

- **Spec.link.ObsAsy** :  Computes observed or asymptotic standardized specialization  with order q.


## Function for Multi-networks:

- **iNEXTbeta.link** : Computing standardized gamma, alpha, beta diversity, and four dissimilarity-turnover indices for three dimensions: taxonomic, phylogenetic and functional diversity at specified sample coverage.


## Functions for Visualizing Results:

- **ggCompleteness.link** : Visualizing the output from the function `Completeness.link`

- **ggSpec.link** : Visualizing the output from the function `Spec.link.est` and `Spec.link.ObsAsy`

- **ggObsAsy.link** : Visualizing the output from the function `ObsAsy.link`

- **ggiNEXT.link** : Visualizing the output from the function `iNEXT.link`
    
- **ggiNEXTbeta.link** : Visualizing the output from the function `iNEXTbeta.link`



First, we load data from `iNEXT.link`:

```{r include=FALSE}
beetles = iNEXT.link::beetles_plotA
beetles = list(Closed = beetles$Closed, Open = beetles$Open)
beetles = lapply(beetles, as.data.frame)
beetles_col_tree = iNEXT.link::beetles_col_tree
beetles_col_distM = iNEXT.link::beetles_col_distM
```

#### SINGLE COMMUNITY FUNCTION: iNEXT.link()

We first describe the main function `iNEXT.link()` with default arguments: 

```{r eval=FALSE}
iNEXT.link(data,diversity = "TD", q = c(0, 1, 2), size = NULL, nT = NULL,
           endpoint = NULL, knots = 40, conf = 0.95, nboot = 30, 
           row.tree = NULL, col.tree = NULL, PDtype = "meanPD", 
           row.distM = NULL, col.distM = NULL, FDtype = "AUC", FDtau = NULL)
```

The arguments of this function are briefly described below, and will be explained in more details by illustrative examples in later text.This main function computes diversity estimates of order q, the sample coverage estimates and related statistics for K (if `knots = K`) evenly-spaced knots (sample sizes) between size 1 and the `endpoint`, where the endpoint is described below. Each knot represents a particular sample size for which diversity estimates will be calculated. By default, endpoint = double the reference sample size (total sample size for interaction data (so calles abundance data in iNEXT.3D)). For example, if `endpoint = 10`, `knot = 4`, diversity estimates will be computed for a sequence of samples with sizes (1, 4, 7, 10).  

```{r, echo=FALSE,warning=FALSE}
Des <- c("data"," a list of data.frames, each data.frames represents col.species-by-row.species abundance matrix.",
"diversity"," selection of diversity type: 'TD' = Taxonomic diversity, 'PD' = Phylogenetic diversity, and 'FD' = Functional diversity.",
"q","	a numerical vector specifying the diversity orders. Default is c(0, 1, 2).",
"size","an integer vector of sample sizes for which diversity estimates will be computed. If NULL, then diversity estimates will be calculated for those sample sizes determined by the specified/default endpoint and knots.",
"endpoint","an integer specifying the sample size that is the endpoint for R/E calculation; If NULL, then endpoint=double the reference sample size;",
"knots","	an integer specifying the number of equally-spaced knots between size 1 and the endpoint. Default is 40.",
"conf","	a positive number < 1 specifying the level of confidence interval. Default is 0.95.",
"nboot"," a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 30.",
"row.tree","(required only when diversity = 'PD' a phylogenetic tree of row assemblage in the pooled network row assemblage. ",
"col.tree","	(required only when diversity = 'PD') a phylogenetic tree of column assemblage in the pooled network column assemblage. ",
"PDtype"," (required only when diversity = 'PD'), select PD type: PDtype = 'PD' (effective total branch length) or PDtype = 'meanPD' (effective number of equally divergent lineages). Default is 'meanPD', where meanPD = PD/tree depth.",
"row.distM"," (required only when </code>diversity = 'FD') a species pairwise distance matrix for all species of row assemblage in the pooled network row assemblage. ",
"col.distM"," (required only when diversity = 'FD') a species pairwise distance matrix for all species of column assemblage in the pooled network column assemblage.",
"FDtype"," (required only when diversity = 'FD'), select FD type: FDtype = 'tau_values' for FD under specified threshold values, or FDtype = 'AUC' (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is 'AUC'.",
"FDtau"," (required only when diversity = 'FD' and FDtype = 'tau_values'), a numerical vector between 0 and 1 specifying tau values (threshold levels). If NULL (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy)."
)
output <- 
  matrix(Des, 
         ncol=2, byrow = TRUE)

library(htmlTable)
htmlTable(output,
          header =  c("Argument","Description"),align = "l"
          )
```


## DATA FORMAT/INFORMATION

Supported Data Types:

Individual-based interaction data : Input data matrix for each assemblage/site include samples species interactions in an empirical sample of n total interactions ("reference sample"). When dealing with N networks, the input data consists of N lists of species interaction matrix.


## RAREFACTION/EXTRAPOLATION VIA EXAMPLES

The data set (tree-beetles interaction data) is included in iNEXT.link package. The experiment took place in the Steigerwald forest in Germany, where deadwood objects from six tree species were exposed in open, net, and closed habitats. In each habitat, there are six plots (A, B, C, D, E, F). Saproxilic beetles were sampled using stem emergence traps and classified according to their functional traits. Data from four years were pooled for each plot and habitat, and pairwise distances were computed from the Gower distance. Here, the demonstration only uses data from plot A in each habitat. For these data, the following commands display the sample species interactions and run the `iNEXT.link()` function for three types of diversty  (`"TD"`, `"PD"`, `"FD"` with specified threshold (default is dmean (quadratic entropy)), `"AUC"` which integrates FD from threshold 0 to 1). 

Under taxonomic diversity dimension, `iNEXT.link()` function returns including: `$DataInfo` for summarizing data information; `$iNextEst` for showing diversity estimates along with related statistics for a series of rarefied and extrapolated samples; and `$AsyEst` for showing asymptotic diversity estimates along with related statistics. Result under phylogenetic diversity or functional diversity includes these three parts, too. 
  
`$DataInfo` in TD example, as shown below, returns basic data information.  It can also be presented using function `DataInfo.link()` to get the same result.

Because the three kinds of diversity output are similar, the demo shows TD only.


```{r setup, include=FALSE, echo=FALSE}
knitr::opts_chunk$set(comment="", message=FALSE, warning=FALSE)
library(iNEXT.link)
library(ggplot2)
```

```{r echo=TRUE}
linkoutTD = iNEXT.link(data = beetles, diversity = 'TD', q = c(0,1,2), nboot = 30)
linkoutTD$DataInfo
```

Second part of output from function `iNEXT.link` is diversity estimates and related statistics computed for these 40 knots by default, which locates the reference sample size at the mid-point of the selected knots. The diversity can be based on sample-size-based and sample coverage-based. The first data frame of list `$iNextEst` (as shown below for 'size_based') includes the sample size (`m`), the `Method` (`Rarefaction`, `Observed`, or `Extrapolation`, depending on whether the size `m` is less than, equal to, or greater than the reference sample size), the diversity order (`Order.q`), the diversity estimate of order q (`qD` in TD, `qPD` in PD, `qFD` in FD (under specified thresholds), `qAUC` in FD (area under curve)), the lower and upper confidence limits of diversity (`qD.LCL` and `qD.UCL` in TD, `qPD.LCL` and `qPD.UCL` in PD, `qFD.LCL` and `qFD.UCL` in FD (under specified thresholds), `qAUC.LCL` and `qAUC.UCL` in FD (area under curve)) conditioning on sample size, and the sample coverage estimate (`SC`) along with the lower and upper confidence limits of sample coverage (`SC.LCL`, `SC.UCL`). These sample coverage estimates with confidence intervals are used for plotting the sample completeness curve. It is time consuming for `diversity = FD` and `FDtype = "AUC"`. If the argument `nboot` is greater than zero, then the bootstrap method is applied to obtain the confidence intervals for each diversity and sample coverage estimates. 

Here only show first six rows:

```{r}
head(linkoutTD$iNextEst$size_based)
```

The second data frame of list `$iNextEst` (as shown below for 'coverage_based') includes the sample coverage estimate ('SC'), the sample size (`m`), the `Method` (`Rarefaction`, `Observed`, or `Extrapolation`, depending on whether the size `m` is less than, equal to, or greater than the reference sample size), the diversity order (`Order.q`), the diversity estimate of order q (`qD` in TD, `qPD` in PD, `qFD` in FD (under specified thresholds), `qAUC` in FD (area under curve)), the lower and upper confidence limits of diversity (`qD.LCL` and `qD.UCL` in TD, `qPD.LCL` and `qPD.UCL` in PD, `qFD.LCL` and `qFD.UCL` in FD (under specified thresholds), `qAUC.LCL` and `qAUC.UCL` in FD (area under curve)) conditioning on sample coverage estimate.

Here only show first six rows:

```{r}
head(linkoutTD$iNextEst$coverage_based)
```

The output `$AsyEst` lists the diversity labels (`Diversity` in TD, `Phylogenetic Diversity` in PD, `Functional Diversity` in FD), the observed diversity (`Observed` in TD, `Phylogenetic Observed` in PD, `Functional Observed` in FD), asymptotic diversity estimates  (`Estimator` in TD, `Phylogenetic Estimator` in PD, `Functional Estimator` in FD), estimated bootstrap standard error (`s.e.`) and  confidence intervals for diversity with q = 0, 1, and 2 (`LCL`, `UCL`). The estimated asymptotic and observed diversity can also be computed via the function `ObsAsy.link()`. The output are shown below:

Here only show first six rows:

```{r}
head(linkoutTD$AsyEst)
```


### GRAPHIC DISPLAYS: FUNCTION ggiNEXT.link()

The function `ggiNEXT.link()`, which extends `ggplot2` with default arguments, is described as follows: 

```{r eval=FALSE}
ggiNEXT.link(outcome, type = 1:3, facet.var = "Assemblage", color.var = "Order.q")  
```

Here `outcome` is the object of `iNEXT.link()`'s output. Three types of curves are allowed for different diversity dimensions:  

(1) Sample-size-based R/E curve (`type = 1`): This curve plots diversity estimates with confidence intervals as a function of sample size.  

(2) Sample completeness curve (`type = 2`): This curve plots the sample coverage with respect to sample size. 

(3) Coverage-based R/E curve (`type = 3`): This curve plots the diversity estimates with confidence intervals as a function of sample coverage. 

The argument `facet.var = "Order.q"` or `facet.var = "Assemblage"` is used to create a separate plot for each value of the specified variable. For example, the following code displays a separate plot of the diversity order q. The `ggiNEXT.link()` function is a wrapper with package `ggplot2` to create a R/E curve in a single line of code. The figure object is of class `"ggplot"`, so can be manipulated by using the `ggplot2` tools. 

When `facet.var = "Assemblage"` in `ggiNEXT.link` function, it creates a separate plot for each network and the different color lines represent each diversity order. Sample-size-based R/E curve (`type = 1`) as below:

```{r }
# Sample-size-based R/E curves, separating by "assemblage""
ggiNEXT.link(linkoutTD, type = 1, facet.var = "Assemblage")
```

When `facet.var = "Order.q"` in `ggiNEXT.link` function, it creates a separate plot for each diversity order and the different color lines represent each network. Sample-size-based R/E curve (`type = 1`) as below:

```{r }
# Sample-size-based R/E curves, separating by "Order.q"
ggiNEXT.link(linkoutTD, type = 1, facet.var = "Order.q")
```


The following command return the sample completeness (sample coverage) curve (`type = 2`) in which different colors are used for the three networks. 

```{r}
ggiNEXT.link(linkoutTD, type = 2, facet.var = "Order.q", color.var = "Assemblage")
```


The following commands return the coverage-based R/E sampling curves in which different colors are used for the three assemblages (`facet.var = "Assemblage"`) and for three diversity orders (`facet.var = "Order.q"`).

```{r}
ggiNEXT.link(linkoutTD, type = 3, facet.var = "Assemblage")
```

```{r}
ggiNEXT.link(linkoutTD, type = 3, facet.var = "Order.q")
```


### DATA INFORMATION FUNCTION: DataInfo.link()

```{r eval=FALSE}
DataInfo.link(data, diversity = "TD", row.tree = NULL, 
              col.tree = NULL, row.distM = NULL, col.distM = NULL) 
```

Here provide the function `DataInfo.link` to compute three diversity dimensions ('TD', 'PD', 'FD') data information, which including sample size, observed species richness, sample coverage estimate, and the first ten interaction frequency counts when `diversity = TD`. And so on for PD, FD.

```{r}
DataInfo.link(beetles, diversity = 'TD')
DataInfo.link(beetles, diversity = 'PD', col.tree = beetles_col_tree)
DataInfo.link(beetles, diversity = 'FD', col.distM = beetles_col_distM)
```


### POINT ESTIMATION FUNCTION: estimateD.link()

```{r eval=FALSE}
estimateD.link(data, diversity = "TD", q = c(0, 1, 2),  
               base = "coverage", level = NULL, nboot = 50, conf = 0.95, 
               PDtype = "meanPD", row.tree = NULL, col.tree = NULL, row.distM = NULL, 
               col.distM = NULL, FDtype = "AUC", FDtau = NULL) 
```

`estimateD.link` is used to compute three diversity dimensions (TD, PD, FD) estimates with q = 0, 1, 2 under any specified level of sample size (when `base = "size"`) or sample coverage (when `base = "coverage"`). If `level = NULL`, this function computes the diversity estimates for the minimum sample size among all samples extrapolated to double reference sizes (when `base = "size"`) or the minimum sample coverage among all samples extrapolated to double reference sizes (when `base = "coverage"`). 

For example, the following command returns the taxonomic diversity ('TD') with a specified level of sample coverage of 95% for the tree-beetles interaction data. For some networks, this coverage value corresponds to the rarefaction part whereas the others correspond to extrapolation, as indicated in the method of the output. 

```{r}
estimateD.link(beetles, diversity = 'TD', q = c(0, 1, 2), base = "coverage", level = 0.95)
```

### ASYMPTOTIC AND OBSERVED DIVERSITY FUNCTION: ObsAsy.link()

```{r,eval=FALSE}
ObsAsy.link(data, diversity = "TD", q = seq(0, 2, 0.2),  nboot = 30, 
            conf = 0.95, method = c("Asymptotic", "Observed"), 
            row.tree = NULL, col.tree = NULL, PDtype = "meanPD", row.distM = NULL, 
            col.distM = NULL, FDtype = "AUC", FDtau = NULL)
```

The function `ObsAsy.link()` compute three diversity dimensions (TD, PD, FD) for empirical (observed) diversity and estimated asymptotic diversity with any diversity order. For example, the following commands returns empirical and asymptotic taxonomic diversity ('TD') for dunes data, along with its confidence interval at diversity order q from 0 to 2. Here only show the first ten rows.

```{r}
out1 <- ObsAsy.link(beetles, diversity = 'TD', q = seq(0, 2, 0.2), method = c("Asymptotic", "Observed"),
                nboot = 5,conf = 0.95)

out1
```


### GRAPHIC DISPLAYS FUNCTION: ggObsAsy.link()

```{r,eval=FALSE}
ggObsAsy.link(outcome)
```

`ggObsAsy.link()` plots q-profile based on `ggplot2`. Here `outcome` is the object from the function `ObsAsy.link`.  

```{r }
# q profile curve
ggObsAsy.link(out1)
```


### SINGLE COMMUNITY FUNCTION: Completeness.link()

Function `Completeness.link()` provides a easy way to compute estimated sample completeness with order q. The arguments is below: 

```{r eval=FALSE}
Completeness.link(data, q = seq(0, 2, 0.2), nboot = 30, conf = 0.95) 
```


### GRAPHIC DISPLAYS FUNCTION:  ggCompleteness.link()

We also provides a realized function `ggCompleteness.link` to plot the output from `Completeness.link()`:

```{r eval=FALSE}
ggCompleteness.link(output)
```


Use data beetles to calculate sample completeness and plot it.

```{r fig.height=8, out.width="70%"}
out1 <- Completeness.link(data = beetles)
ggCompleteness.link(out1)
```



### SINGLE COMMUNITY FUNCTION: Spec.link()

The main function `Spec.link()` with default arguments: 

```{r eval=FALSE}
Spec.link(data, q = seq(0, 2, 0.2), method = "Estimated",
         nboot = 30, conf = 0.95, E.class = c(1:5), C = NULL) 
```


### GRAPHIC DISPLAYS FUNCTION: ggSpec.link()

The function `ggSpec.link()` is provided to plot the output from `Spec.link()`.

```{r eval=FALSE}
ggSpec.link(output)
```

There is an example for function `Spec.link` and function `ggSpec.link`. 

```{r}
out1 <- Spec.link(data = beetles)
ggSpec.link(out1)
```


### MULTI-COMMUNITIES FUNCTION: iNEXTbeta.link()

```{r eval=FALSE}
iNEXTbeta.link(data, diversity = "TD", level = seq(0.5, 1, 0.05), 
               q = c(0, 1, 2), nboot = 20, conf = 0.95, 
               PDtype = "meanPD", row.tree = NULL, col.tree = NULL, row.distM = NULL, 
               col.distM = NULL, FDtype = "AUC", FDtau = NULL, FDcut_number = 30) 
```

The arguments of this function are briefly described below, and will be explained in more details by illustrative examples in later text. This main function computes gamma, alpha and beta diversity estimates of order q at specified sample coverage and measures of dissimilarity. By default of <code>base = "coverage"</code> and <code>level = NULL</code>, then this function computes the gamma, alpha, beta diversity, and four dissimilarity-turnover indices estimates up to one (for q = 1, 2) or up to the coverage of double the reference sample size (for q = 0).  



```{r, echo=FALSE,warning=FALSE}
Des <- c("data"," data can be input as a matrix/data.frame (species by assemblages), or a list of matrices/data.frames, each matrix represents species-by-assemblages abundance matrix",
"diversity"," selection of diversity type: 'TD' = Taxonomic diversity, 'PD' = Phylogenetic diversity, and 'FD' = Functional diversity.",
"level","a sequence specifying the particular sample coverages (between 0 and 1). Default is seq(0.5, 1, 0.05).",
"q","	a numerical vector specifying the diversity orders. Default is c(0, 1, 2).",
"nboot"," a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter 0 to skip the bootstrap procedures. Default is 30.",
"conf","a positive number < 1 specifying the level of confidence interval. Default is 0.95.",
"PDtype"," (required only when diversity = 'PD'), select PD type: PDtype = 'PD' (effective total branch length) or PDtype = 'meanPD' (effective number of equally divergent lineages). Default is 'meanPD', where meanPD = PD/tree depth.",
"row.tree","(required only when diversity = 'PD' a phylogenetic tree of row assemblage in the pooled network row assemblage. ",
"col.tree","	(required only when diversity = 'PD') a phylogenetic tree of column assemblage in the pooled network column assemblage. ",
"row.distM"," (required only when </code>diversity = 'FD') a species pairwise distance matrix for all species of row assemblage in the pooled network row assemblage. ",
"col.distM"," (required only when diversity = 'FD') a species pairwise distance matrix for all species of column assemblage in the pooled network column assemblage.",
"FDtype"," (required only when diversity = 'FD'), select FD type: FDtype = 'tau_values' for FD under specified threshold values, or FDtype = 'AUC' (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is 'AUC'.",
"FDtau"," (required only when diversity = 'FD' and FDtype = 'tau_values'), a numerical vector between 0 and 1 specifying tau values (threshold levels). If NULL (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy)."
)
output <- 
  matrix(Des, 
         ncol=2, byrow = TRUE)

library(htmlTable)
htmlTable(output,
          header =  c("Argument","Description"),align = "l"
          )
```

the `iNEXTbeta.link()` function returns the `"iNEXTbeta.link"` object including seven data frames for each datasets: 

- gamma 
- alpha 
- beta
- C ( Sorensen-type non-overlap )
- U ( Jaccard-type  non-overlap )
- V ( Sorensen-type turnover )
- S ( Jaccard-type  turnover ) 

Here only show the first six rows in each table output:

```{r eval=FALSE}
# Taxonomic diversity
Abundance_TD = iNEXTbeta.link(data = beetles, diversity = 'TD', level = NULL, q = c(0, 1, 2))
Abundance_TD
```


```{r echo=FALSE}
Abundance_TD = iNEXTbeta.link(data = beetles, diversity = 'TD', level = NULL, q = c(0, 1, 2))

lapply(Abundance_TD$Dataset_1, function(x) {
    tmp = x[1:6,]
    tmp[,c(3,4,5,7,8,9)] = round(tmp[,c(3,4,5,7,8,9)], 3)
    tmp
})
```

The output contains seven data frames: `gamma`, `alpha`, `beta`, `C`, `U`, `V`, `S`. For each data frame, it includes the diversity estimate (`Estimate`), the diversity order (`Order.q`), `Method` (Rarefaction, Observed, or Extrapolation, depending on whether the size `m` is less than, equal to, or greater than the reference sample size), the sample coverage estimate (`SC`), the sample size (`Size`), the standard error from bootstrap replications (`s.e.`), the lower and upper confidence limits of diversity (`LCL`, `UCL`), and the name of dataset (`Dataset`). These diversity estimates with confidence intervals are used for plotting the diversity curve.


### GRAPHIC DISPLAYS FUNCTION: ggiNEXTbeta.link()

The function `ggiNEXTbeta.link()`, which extends `ggplot2` to the `"iNEXT.link"` object with default arguments, is described as follows: 


```{r, echo=FALSE,warning=FALSE}
Des <- c("output"," the output of <code>iNEXTbeta.link</code>.",
"type"," selection of plot type : type = 'B' for plotting the gamma, alpha, and beta diversity; type = 'D' for plotting 4 turnover dissimilarities.",
"scale"," Are scales shared across all facets (the default, 'fixed'), or do they vary across rows ('free_x'), columns ('free_y'), or both rows and columns ('free')? Default is 'free'."
)
output <- 
  matrix(Des, 
         ncol=2, byrow = TRUE)

library(htmlTable)
htmlTable(output,
          header =  c("Argument","Description"),align = "l"
          )
```


The `ggiNEXTbeta.link()` function is a wrapper around the `ggplot2` package to create a R/E curve using a single line of code. The resulting object is of class `"ggplot"`, so it can be manipulated using the `ggplot2` tools. Users can visualize the output of beta diversity or four dissimilarities by setting the parameter <code>**type**</code>:

```{r, fig.align='center', fig.height=6, fig.width=6}
ggiNEXTbeta.link(Abundance_TD, type = 'B')
```

```{r, fig.align='center', fig.height=8, fig.width=6}
ggiNEXTbeta.link(Abundance_TD, type = 'D')
```


### HOW TO CITE iNEXT.link
If you publish your work based on the results from the `iNEXT.link` package, you should make references to the following methodology paper:

- Chiu, C-H., Chao, A., Vogel, S., Kriegel, P. and Thorn, S. (2023). Quantifying and estimating ecological network diversity based on incomplete sampling data. Philosophical Transactions of the Royal Society B, 378: 20220183. https://doi.org/10.1098/rstb.2022.0183


### License
The iNEXT.link package is licensed under the GPLv3. To help refine `iNEXT.link`, your comments or feedback would be welcome (please send them to Anne Chao or report an issue on the iNEXT.link github [iNEXT.link_github](https://github.com/AnneChao/iNEXT.link). 

### References
- Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H., Dornelas, M and Magurran, A. E. (2021). Measuring temporal change in alpha diversity: a framework integrating taxonomic, phylogenetic and functional diversity and the iNEXT.3D standardization. Methods in Ecology and Evolution, 12, 1926-1940.

- Chao, A. & Jost, L. (2012) Coverage-based rarefaction and extrapolation: standardizing samples by completeness rather than size. Ecology, 93, 2533ÅM2547

- Chao, A., Y. Kubota, D. Zelen??, C.-H. Chiu, C.-F. Li, B. Kusumoto, M. Yasuhara, S. Thorn, C.-L. Wei, M. J. Costello, and R. K. Colwell (2020). Quantifying sample completeness and comparing diversities among assemblages. Ecological Research, 35, 292-314.

- Chiu, C-H., Chao, A., Vogel, S., Kriegel, P. and Thorn, S. (2023). Quantifying and estimating ecological network diversity based on incomplete sampling data. Philosophical Transactions of the Royal Society B, 378: 20220183. https://doi.org/10.1098/rstb.2022.0183


