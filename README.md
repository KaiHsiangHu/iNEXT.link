<!-- README.md is generated from README.Rmd. Please edit that file -->

# iNEXT.link (R package)

<h5 align="right">
Latest version: 2024-01-09
</h5>
<font color="394CAE">
<h3 color="394CAE" style="font-weight: bold">
Introduction to iNEXT.link (R package): Excerpt from iNEXT.link
UserGuide
</h3>
</font> <br>
<h5>
<b>Anne Chao, Kai-Hsiang Hu, K.W. Chen, C.G. Lo, S.Y. Wang</b> <br><br>
<i>Institute of Statistics, National Tsing Hua University, Hsin-Chu,
Taiwan 30043</i>
</h5>

<br> `iNEXT.link` (INterpolation and EXTrapolation in Network diversity)
is an R package, available in [Github](https://github.com/AnneChao).
Here we provide a quick introduction demonstrating how to run
iNEXT.link. An online version of [iNEXT.link
Online](https://chao.shinyapps.io/iNEXT_link/) is also available for
users without an R background. Detailed information about all functions
in iNEXT.link is provided in the iNEXT.link Manual in
[iNEXT.link_vignettes](http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/A%20Quick%20Introduction%20to%20iNEXT.link%20via%20Examples.html),
which is also available from [Anne Chao’s
website](http://chao.stat.nthu.edu.tw/wordpress/software_download/).

`iNEXT.link` is an R package that extends the concepts of `iNEXT.3D`
(Chao et al., 2021), `iNEXT.4step` (Chao et al., 2020) and
`iNEXT.beta3D` (Chao et al., 2023) to ecological networks (Chiu et al.,
2023). It is primarily designed to calculate and analyze various
measures of diversity in ecological networks. Specifically, the package
calculates three Hill numbers of order q (species richness, Shannon
diversity, and Simpson diversity) in taxonomic diversity level, as well
as phylogenetic and functional diversity levels.

For single ecological networks, `iNEXT.link` provides tools for
analyzing diversity. The package provides two types of rarefaction and
extrapolation (R/E) sampling curves to estimate diversity and confidence
intervals for single ecological networks. These include
sample-size-based (or size-based) R/E curves and coverage-based R/E
curves.

Moreover, `iNEXT.link` offers dissimilarity-turnover curves for the
coverage-based R/E curves for gamma, alpha, and beta diversity measures,
which can be used to compare diversity patterns across different
ecological networks.

### SOFTWARE NEEDED TO RUN INEXT.3D IN R

-   Required: [R](http://cran.rstudio.com/)
-   Suggested: [RStudio IDE](http://www.rstudio.com/ide/download/)

### HOW TO RUN iNEXT.link:

The `iNEXT.link` package can be downloaded from Anne Chao’s
[iNEXT.link_github](https://github.com/AnneChao/iNEXT.link) using the
following commands. For a first-time installation, additional
visualization extension packages (`ggplot2`) from CRAN and (`iNEXT.3D`),
(`iNEXT.4steps`), and (`iNEXT.beta3D`) from Anne Chao’s github must be
installed and loaded.

``` r
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

In this document, we provide a quick introduction demonstrating how to
run the package `iNEXT.link`(iNterpolation and EXTrapolation in Network
diversity). `iNEXT.link` has several main functions:

## Functions for Single community:

-   **iNEXT.link** : Computes rarefaction/extrapolation taxonomic,
    phylogenetic, and functional diversity estimates and sample coverage
    estimates.

-   **DataInfo.link** : exhibits basic data information

-   **estimateD.link** : computes species diversity with a particular
    user-specified level of sample size or sample coverage.

-   **ObsAsy.link**: compute asymptotic or empirical(observed) diversity
    of order q.

-   **Completeness.link** : Calculates estimated sample completeness
    with order q.

-   **Spec.link** : Computes standardized specialization estimation
    under specified sample coverage(or observed) with order q.

## Function for Multi-community:

-   **iNEXTbeta.link** : Computing standardized gamma, alpha, beta
    diversity, and four dissimilarity-turnover indices for three
    dimensions: taxonomic, phylogenetic and functional diversity at
    specified sample coverage.

## Functions for Visualizing Results:

-   **ggCompleteness.link** : Visualizing the output from the function
    `Completeness.link`

-   **ggSpec.link** : Visualizing the output from the function
    `Spec.link`

-   **ggObsAsy.link** : Visualizing the output from the function
    `ObsAsy.link`

-   **ggiNEXT.link** : Visualizing the output from the function
    `iNEXT.link`

-   **ggiNEXTbeta.link** : Visualizing the output from the function
    `iNEXTbeta.link`

First, we load data from `iNEXT.link`:

#### SINGLE COMMUNITY FUNCTION: iNEXT.link()

We first describe the main function `iNEXT.link()` with default
arguments:

``` r
iNEXT.link(data,diversity = "TD", q = c(0, 1, 2), size = NULL, nT = NULL,
           endpoint = NULL, knots = 40, conf = 0.95, nboot = 30, 
           row.tree = NULL, col.tree = NULL, PDtype = "meanPD", 
           row.distM = NULL, col.distM = NULL, FDtype = "AUC", FDtau = NULL)
```

The arguments of this function are briefly described below, and will be
explained in more details by illustrative examples in later text.This
main function computes diversity estimates of order q, the sample
coverage estimates and related statistics for K (if `knots = K`)
evenly-spaced knots (sample sizes) between size 1 and the `endpoint`,
where the endpoint is described below. Each knot represents a particular
sample size for which diversity estimates will be calculated. By
default, endpoint = double the reference sample size (total sample size
for interaction data (so calles abundance data in iNEXT.3D)). For
example, if `endpoint = 10`, `knot = 4`, diversity estimates will be
computed for a sequence of samples with sizes (1, 4, 7, 10).

<table class="gmisc_table" style="border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;">
<thead>
<tr>
<th style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">
Argument
</th>
<th style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">
Description
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">
data
</td>
<td style="text-align: left;">
a list of data.frames, each data.frames represents
col.species-by-row.species abundance matrix.
</td>
</tr>
<tr>
<td style="text-align: left;">
diversity
</td>
<td style="text-align: left;">
selection of diversity type: ‘TD’ = Taxonomic diversity, ‘PD’ =
Phylogenetic diversity, and ‘FD’ = Functional diversity.
</td>
</tr>
<tr>
<td style="text-align: left;">
q
</td>
<td style="text-align: left;">
a numerical vector specifying the diversity orders. Default is c(0, 1,
2).
</td>
</tr>
<tr>
<td style="text-align: left;">
size
</td>
<td style="text-align: left;">
an integer vector of sample sizes for which diversity estimates will be
computed. If NULL, then diversity estimates will be calculated for those
sample sizes determined by the specified/default endpoint and knots.
</td>
</tr>
<tr>
<td style="text-align: left;">
endpoint
</td>
<td style="text-align: left;">
an integer specifying the sample size that is the endpoint for R/E
calculation; If NULL, then endpoint=double the reference sample size;
</td>
</tr>
<tr>
<td style="text-align: left;">
knots
</td>
<td style="text-align: left;">
an integer specifying the number of equally-spaced knots between size 1
and the endpoint. Default is 40.
</td>
</tr>
<tr>
<td style="text-align: left;">
conf
</td>
<td style="text-align: left;">
a positive number \< 1 specifying the level of confidence interval.
Default is 0.95.
</td>
</tr>
<tr>
<td style="text-align: left;">
nboot
</td>
<td style="text-align: left;">
a positive integer specifying the number of bootstrap replications when
assessing sampling uncertainty and constructing confidence intervals.
Enter 0 to skip the bootstrap procedures. Default is 30.
</td>
</tr>
<tr>
<td style="text-align: left;">
row.tree
</td>
<td style="text-align: left;">
(required only when diversity = ‘PD’ a phylogenetic tree of row
assemblage in the pooled network row assemblage.
</td>
</tr>
<tr>
<td style="text-align: left;">
col.tree
</td>
<td style="text-align: left;">
(required only when diversity = ‘PD’) a phylogenetic tree of column
assemblage in the pooled network column assemblage.
</td>
</tr>
<tr>
<td style="text-align: left;">
PDtype
</td>
<td style="text-align: left;">
(required only when diversity = ‘PD’), select PD type: PDtype = ‘PD’
(effective total branch length) or PDtype = ‘meanPD’ (effective number
of equally divergent lineages). Default is ‘meanPD’, where meanPD =
PD/tree depth.
</td>
</tr>
<tr>
<td style="text-align: left;">
row.distM
</td>
<td style="text-align: left;">
(required only when </code>diversity = ‘FD’) a species pairwise distance
matrix for all species of row assemblage in the pooled network row
assemblage.
</td>
</tr>
<tr>
<td style="text-align: left;">
col.distM
</td>
<td style="text-align: left;">
(required only when diversity = ‘FD’) a species pairwise distance matrix
for all species of column assemblage in the pooled network column
assemblage.
</td>
</tr>
<tr>
<td style="text-align: left;">
FDtype
</td>
<td style="text-align: left;">
(required only when diversity = ‘FD’), select FD type: FDtype =
‘tau_values’ for FD under specified threshold values, or FDtype = ‘AUC’
(area under the curve of tau-profile) for an overall FD which integrates
all threshold values between zero and one. Default is ‘AUC’.
</td>
</tr>
<tr>
<td style="border-bottom: 2px solid grey; text-align: left;">
FDtau
</td>
<td style="border-bottom: 2px solid grey; text-align: left;">
(required only when diversity = ‘FD’ and FDtype = ‘tau_values’), a
numerical vector between 0 and 1 specifying tau values (threshold
levels). If NULL (default), then threshold is set to be the mean
distance between any two individuals randomly selected from the pooled
assemblage (i.e., quadratic entropy).
</td>
</tr>
</tbody>
</table>

## DATA FORMAT/INFORMATION

Supported Data Types:

Individual-based interaction data : Input data matrix for each
assemblage/site include samples species interactions in an empirical
sample of n total interactions (“reference sample”). When dealing with N
networks, the input data consists of N lists of species interaction
matrix.

## RAREFACTION/EXTRAPOLATION VIA EXAMPLES

The data set (tree-beetles interaction data) is included in iNEXT.link
package. The experiment took place in the Steigerwald forest in Germany,
where deadwood objects from six tree species were exposed in open, net,
and closed habitats. In each habitat, there are six plots (A, B, C, D,
E, F). Saproxilic beetles were sampled using stem emergence traps and
classified according to their functional traits. Data from four years
were pooled for each plot and habitat, and pairwise distances were
computed from the Gower distance. Here, the demonstration only uses data
from plot A in each habitat. For these data, the following commands
display the sample species interactions and run the `iNEXT.link()`
function for three types of diversty (`"TD"`, `"PD"`, `"FD"` with
specified threshold (default is dmean (quadratic entropy)), `"AUC"`
which integrates FD from threshold 0 to 1).

Under taxonomic diversity dimension, `iNEXT.link()` function returns
including: `$DataInfo` for summarizing data information; `$iNextEst` for
showing diversity estimates along with related statistics for a series
of rarefied and extrapolated samples; and `$AsyEst` for showing
asymptotic diversity estimates along with related statistics. Result
under phylogenetic diversity or functional diversity includes these
three parts, too.

`$DataInfo` in TD example, as shown below, returns basic data
information. It can also be presented using function `DataInfo.link()`
to get the same result.

Because the three kinds of diversity output are similar, the demo shows
TD only.

``` r
linkoutTD = iNEXT.link(data = beetles, diversity = 'TD', q = c(0,1,2), nboot = 30)
linkoutTD$DataInfo
NULL
```

Second part of output from function `iNEXT.link` is diversity estimates
and related statistics computed for these 40 knots by default, which
locates the reference sample size at the mid-point of the selected
knots. The diversity can be based on sample-size-based and sample
coverage-based. The first data frame of list `$iNextEst` (as shown below
for ‘size_based’) includes the sample size (`m`), the `Method`
(`Rarefaction`, `Observed`, or `Extrapolation`, depending on whether the
size `m` is less than, equal to, or greater than the reference sample
size), the diversity order (`Order.q`), the diversity estimate of order
q (`qD` in TD, `qPD` in PD, `qFD` in FD (under specified thresholds),
`qAUC` in FD (area under curve)), the lower and upper confidence limits
of diversity (`qD.LCL` and `qD.UCL` in TD, `qPD.LCL` and `qPD.UCL` in
PD, `qFD.LCL` and `qFD.UCL` in FD (under specified thresholds),
`qAUC.LCL` and `qAUC.UCL` in FD (area under curve)) conditioning on
sample size, and the sample coverage estimate (`SC`) along with the
lower and upper confidence limits of sample coverage (`SC.LCL`,
`SC.UCL`). These sample coverage estimates with confidence intervals are
used for plotting the sample completeness curve. It is time consuming
for `diversity = FD` and `FDtype = "AUC"`. If the argument `nboot` is
greater than zero, then the bootstrap method is applied to obtain the
confidence intervals for each diversity and sample coverage estimates.

Here only show first six rows:

``` r
head(linkoutTD$iNextEst$size_based)
NULL
```

The second data frame of list `$iNextEst` (as shown below for
‘coverage_based’) includes the sample coverage estimate (‘SC’), the
sample size (`m`), the `Method` (`Rarefaction`, `Observed`, or
`Extrapolation`, depending on whether the size `m` is less than, equal
to, or greater than the reference sample size), the diversity order
(`Order.q`), the diversity estimate of order q (`qD` in TD, `qPD` in PD,
`qFD` in FD (under specified thresholds), `qAUC` in FD (area under
curve)), the lower and upper confidence limits of diversity (`qD.LCL`
and `qD.UCL` in TD, `qPD.LCL` and `qPD.UCL` in PD, `qFD.LCL` and
`qFD.UCL` in FD (under specified thresholds), `qAUC.LCL` and `qAUC.UCL`
in FD (area under curve)) conditioning on sample coverage estimate.

Here only show first six rows:

``` r
head(linkoutTD$iNextEst$coverage_based)
NULL
```

The output `$AsyEst` lists the diversity labels (`Diversity` in TD,
`Phylogenetic Diversity` in PD, `Functional Diversity` in FD), the
observed diversity (`Observed` in TD, `Phylogenetic Observed` in PD,
`Functional Observed` in FD), asymptotic diversity estimates
(`Estimator` in TD, `Phylogenetic Estimator` in PD,
`Functional Estimator` in FD), estimated bootstrap standard error
(`s.e.`) and confidence intervals for diversity with q = 0, 1, and 2
(`LCL`, `UCL`). The estimated asymptotic and observed diversity can also
be computed via the function `ObsAsy.link()`. The output are shown
below:

Here only show first six rows:

``` r
head(linkoutTD$AsyEst)
NULL
```

### GRAPHIC DISPLAYS: FUNCTION ggiNEXT.link()

The function `ggiNEXT.link()`, which extends `ggplot2` with default
arguments, is described as follows:

``` r
ggiNEXT.link(outcome, type = 1:3, facet.var = "Assemblage", color.var = "Order.q")  
```

Here `outcome` is the object of `iNEXT.link()`’s output. Three types of
curves are allowed for different diversity dimensions:

1.  Sample-size-based R/E curve (`type = 1`): This curve plots diversity
    estimates with confidence intervals as a function of sample size.

2.  Sample completeness curve (`type = 2`): This curve plots the sample
    coverage with respect to sample size.

3.  Coverage-based R/E curve (`type = 3`): This curve plots the
    diversity estimates with confidence intervals as a function of
    sample coverage.

The argument `facet.var = "Order.q"` or `facet.var = "Assemblage"` is
used to create a separate plot for each value of the specified variable.
For example, the following code displays a separate plot of the
diversity order q. The `ggiNEXT.link()` function is a wrapper with
package `ggplot2` to create a R/E curve in a single line of code. The
figure object is of class `"ggplot"`, so can be manipulated by using the
`ggplot2` tools.

When `facet.var = "Assemblage"` in `ggiNEXT.link` function, it creates a
separate plot for each network and the different color lines represent
each diversity order. Sample-size-based R/E curve (`type = 1`) as below:

``` r
# Sample-size-based R/E curves, separating by "assemblage""
ggiNEXT.link(linkoutTD, type = 1, facet.var = "Assemblage")
[[1]]
```

<img src="README/README-unnamed-chunk-11-1.png" width="672" style="display: block; margin: auto;" />

When `facet.var = "Order.q"` in `ggiNEXT.link` function, it creates a
separate plot for each diversity order and the different color lines
represent each network. Sample-size-based R/E curve (`type = 1`) as
below:

``` r
# Sample-size-based R/E curves, separating by "Order.q"
ggiNEXT.link(linkoutTD, type = 1, facet.var = "Order.q")
[[1]]
```

<img src="README/README-unnamed-chunk-12-1.png" width="672" style="display: block; margin: auto;" />

The following command return the sample completeness (sample coverage)
curve (`type = 2`) in which different colors are used for the three
networks.

``` r
ggiNEXT.link(linkoutTD, type = 2, facet.var = "Order.q", color.var = "Assemblage")
[[1]]
```

<img src="README/README-unnamed-chunk-13-1.png" width="672" style="display: block; margin: auto;" />

The following commands return the coverage-based R/E sampling curves in
which different colors are used for the three assemblages
(`facet.var = "Assemblage"`) and for three diversity orders
(`facet.var = "Order.q"`).

``` r
ggiNEXT.link(linkoutTD, type = 3, facet.var = "Assemblage")
[[1]]
```

<img src="README/README-unnamed-chunk-14-1.png" width="672" style="display: block; margin: auto;" />

``` r
ggiNEXT.link(linkoutTD, type = 3, facet.var = "Order.q")
[[1]]
```

<img src="README/README-unnamed-chunk-15-1.png" width="672" style="display: block; margin: auto;" />

### DATA INFORMATION FUNCTION: DataInfo.link()

``` r
DataInfo.link(data, diversity = "TD", row.tree = NULL, 
              col.tree = NULL, row.distM = NULL, col.distM = NULL) 
```

Here provide the function `DataInfo.link` to compute three diversity
dimensions (‘TD’, ‘PD’, ‘FD’) data information, which including sample
size, observed species richness, sample coverage estimate, and the first
ten interaction frequency counts when `diversity = TD`. And so on for
PD, FD.

``` r
DataInfo.link(beetles, diversity = 'TD')
  Networks    n S.obs(row) S.obs(col) Links.obs Connectance  Coverage  f1 f2 f3 f4 f5 f6 f7 f8 f9 f10
1   Closed  816          6         83       178      0.3574 0.8799951  98 31 15  3  3  5  0  2  2   1
2     Open 1932          6         88       206      0.3902 0.9430819 110 33 15  8  7  5  3  2  3   2
DataInfo.link(beetles, diversity = 'PD', col.tree = beetles_col_tree)
  Networks    n S.obs(row) S.obs(col) Links.obs Connectance f1* f2*       g1       g2 PD.obs mean_T
1   Closed  816          6         83       178      0.3574 171  76 11508.99 3873.052    924    285
2     Open 1932          6         88       206      0.3902 190  76 13918.45 5172.410   1014    285
DataInfo.link(beetles, diversity = 'FD', col.distM = beetles_col_distM)
  Networks    n S.obs(row) S.obs(col) Links.obs Connectance  f1 f2 a1' a2' threshold
1   Closed  816          6         83       178      0.3574  98 31   0   0  8.930741
2     Open 1932          6         88       206      0.3902 110 33   0   0  8.021019
```

### POINT ESTIMATION FUNCTION: estimateD.link()

``` r
estimateD.link(data, diversity = "TD", q = c(0, 1, 2),  
               base = "coverage", level = NULL, nboot = 50, conf = 0.95, 
               PDtype = "meanPD", row.tree = NULL, col.tree = NULL, row.distM = NULL, 
               col.distM = NULL, FDtype = "AUC", FDtau = NULL) 
```

`estimateD.link` is used to compute three diversity dimensions (TD, PD,
FD) estimates with q = 0, 1, 2 under any specified level of sample size
(when `base = "size"`) or sample coverage (when `base = "coverage"`). If
`level = NULL`, this function computes the diversity estimates for the
minimum sample size among all samples extrapolated to double reference
sizes (when `base = "size"`) or the minimum sample coverage among all
samples extrapolated to double reference sizes (when
`base = "coverage"`).

For example, the following command returns the taxonomic diversity
(‘TD’) with a specified level of sample coverage of 95% for the
tree-beetles interaction data. For some networks, this coverage value
corresponds to the rarefaction part whereas the others correspond to
extrapolation, as indicated in the method of the output.

``` r
estimateD.link(beetles, diversity = 'TD', q = c(0, 1, 2), base = "coverage", level = 0.95)
  Assemblage Order.q   SC        m        Method        qTD       s.e.    qTD.LCL    qTD.UCL
1     Closed       0 0.95 1944.292 Extrapolation 268.252103 19.4789617 230.074040 306.430166
2     Closed       1 0.95 1944.292 Extrapolation  74.577021  3.6666959  67.390429  81.763613
3     Closed       2 0.95 1944.292 Extrapolation  33.378889  1.8831205  29.688040  37.069737
4       Open       0 0.95 2349.132 Extrapolation 228.271774 20.1609372 188.757063 267.786485
5       Open       1 0.95 2349.132 Extrapolation  23.704215  1.3613738  21.035971  26.372458
6       Open       2 0.95 2349.132 Extrapolation   8.065214  0.3508078   7.377643   8.752784
```

### ASYMPTOTIC AND OBSERVED DIVERSITY FUNCTION: ObsAsy.link()

``` r
ObsAsy.link(data, diversity = "TD", q = seq(0, 2, 0.2),  nboot = 30, 
            conf = 0.95, method = c("Asymptotic", "Observed"), 
            row.tree = NULL, col.tree = NULL, PDtype = "meanPD", row.distM = NULL, 
            col.distM = NULL, FDtype = "AUC", FDtau = NULL)
```

The function `ObsAsy.link()` compute three diversity dimensions (TD, PD,
FD) for empirical (observed) diversity and estimated asymptotic
diversity with any diversity order. For example, the following commands
returns empirical and asymptotic taxonomic diversity (‘TD’) for dunes
data, along with its confidence interval at diversity order q from 0 to
2. Here only show the first ten rows.

``` r
out1 <- ObsAsy.link(beetles, diversity = 'TD', q = seq(0, 2, 0.2), method = c("Asymptotic", "Observed"),
                nboot = 5,conf = 0.95)

out1
   Network Order.q        qTD       s.e.    qTD.LCL    qTD.UCL     Method
1   Closed     0.0 332.713393 60.7494794 248.807295 388.150139 Asymptotic
2   Closed     0.2 265.977501 40.0782717 210.483668 304.674929 Asymptotic
3   Closed     0.4 203.167245 23.0455368 171.395956 226.923193 Asymptotic
4   Closed     0.6 149.540048 11.1877823 134.610344 161.910172 Asymptotic
5   Closed     0.8 108.592546  4.7257148 103.413490 114.136074 Asymptotic
6   Closed     1.0  80.316088  2.5887035  77.135131  83.272186 Asymptotic
7   Closed     1.2  61.983531  2.3287095  59.172274  64.715248 Asymptotic
8   Closed     1.4  50.295815  2.1913699  47.830763  52.895399 Asymptotic
9   Closed     1.6  42.699384  1.9902762  40.563109  45.119820 Asymptotic
10  Closed     1.8  37.568079  1.7887645  35.732154  39.803041 Asymptotic
11  Closed     2.0  33.944467  1.6174707  32.375699  36.001151 Asymptotic
12    Open     0.0 389.238440 25.0436417 351.956315 411.302136 Asymptotic
13    Open     0.2 263.223968 16.6619512 239.923045 281.351119 Asymptotic
14    Open     0.4 159.174040 10.5518258 146.956595 171.822803 Asymptotic
15    Open     0.6  86.647814  6.4185253  79.429038  93.640386 Asymptotic
16    Open     0.8  45.571274  3.7189968  40.476046  49.255856 Asymptotic
17    Open     1.0  25.911799  2.1349317  22.779659  27.982319 Asymptotic
18    Open     1.2  16.963237  1.2985247  15.021985  18.161496 Asymptotic
19    Open     1.4  12.643810  0.8629161  11.349884  13.412179 Asymptotic
20    Open     1.6  10.333946  0.6238429   9.397577  10.866772 Asymptotic
21    Open     1.8   8.967960  0.4825937   8.242270   9.361798 Asymptotic
22    Open     2.0   8.089554  0.3928921   7.497158   8.395495 Asymptotic
23  Closed     0.0 178.000000 14.5705182 163.000000 195.000000   Observed
24  Closed     0.2 149.140436 12.9906767 135.699274 164.058644   Observed
25  Closed     0.4 122.284123 11.1339416 110.785042 134.787250   Observed
26  Closed     0.6  98.784711  9.1630339  89.408484 108.761439   Observed
27  Closed     0.8  79.539008  7.2927705  72.198747  87.535873   Observed
28  Closed     1.0  64.690053  5.6997735  59.081390  71.255242   Observed
29  Closed     1.2  53.711702  4.4581426  49.439587  59.094022   Observed
30  Closed     1.4  45.761569  3.5471212  42.454094  50.222956   Observed
31  Closed     1.6  40.008527  2.8997899  37.370939  43.775855   Observed
32  Closed     1.8  35.788959  2.4443288  33.608940  39.038698   Observed
33  Closed     2.0  32.627205  2.1223556  30.759029  35.489163   Observed
34    Open     0.0 206.000000  6.6858059 198.500000 214.600000   Observed
35    Open     0.2 147.911160  4.7333977 142.653699 153.717494   Observed
36    Open     0.4  98.276240  3.0741560  94.997952 101.634857   Observed
37    Open     0.6  60.906771  1.8905086  59.104208  63.299753   Observed
38    Open     0.8  36.832763  1.1698536  35.761650  38.478061   Observed
39    Open     1.0  23.311453  0.7590379  22.496702  24.386806   Observed
40    Open     1.2  16.195155  0.5259785  15.581739  16.917882   Observed
41    Open     1.4  12.394607  0.3925491  11.915707  12.918305   Observed
42    Open     1.6  10.237828  0.3139364   9.841360  10.644900   Observed
43    Open     1.8   8.920980  0.2654315   8.578456   9.257250   Observed
44    Open     2.0   8.059978  0.2338830   7.754363   8.352926   Observed
```

### GRAPHIC DISPLAYS FUNCTION: ggObsAsy.link()

``` r
ggObsAsy.link(outcome)
```

`ggObsAsy.link()` plots q-profile based on `ggplot2`. Here `outcome` is
the object from the function `ObsAsy.link`.

``` r
# q profile curve
ggObsAsy.link(out1)
```

<img src="README/README-unnamed-chunk-23-1.png" width="672" style="display: block; margin: auto;" />

### SINGLE COMMUNITY FUNCTION: Completeness.link()

Function `Completeness.link()` provides a easy way to compute estimated
sample completeness with order q. The arguments is below:

``` r
Completeness.link(data, q = seq(0, 2, 0.2), nboot = 30, conf = 0.95) 
```

### GRAPHIC DISPLAYS FUNCTION: ggCompleteness.link()

We also provides a realized function `ggCompleteness.link` to plot the
output from `Completeness.link()`:

``` r
ggCompleteness.link(output)
```

Use data beetles to calculate sample completeness and plot it.

``` r
out1 <- Completeness.link(data = beetles)
ggCompleteness.link(out1)
```

<img src="README/README-unnamed-chunk-26-1.png" width="70%" style="display: block; margin: auto;" />

### SINGLE COMMUNITY FUNCTION: Spec.link()

The main function `Spec.link()` with default arguments:

``` r
Spec.link(data, q = seq(0, 2, 0.2), method = "Estimated",
         nboot = 30, conf = 0.95, E.class = c(1:5), C = NULL) 
```

### GRAPHIC DISPLAYS FUNCTION: ggSpec.link()

The function `ggSpec.link()` is provided to plot the output from
`Spec.link()`.

``` r
ggSpec.link(output)
```

There is an example for function `Spec.link` and function `ggSpec.link`.

``` r
out1 <- Spec.link(data = beetles)
ggSpec.link(out1)
```

<img src="README/README-unnamed-chunk-29-1.png" width="672" style="display: block; margin: auto;" />

### MULTI-COMMUNITIES FUNCTION: iNEXTbeta.link()

``` r
iNEXTbeta.link(data, diversity = "TD", level = seq(0.5, 1, 0.05), 
               q = c(0, 1, 2), nboot = 20, conf = 0.95, 
               PDtype = "meanPD", row.tree = NULL, col.tree = NULL, row.distM = NULL, 
               col.distM = NULL, FDtype = "AUC", FDtau = NULL, FDcut_number = 30) 
```

The arguments of this function are briefly described below, and will be
explained in more details by illustrative examples in later text. This
main function computes gamma, alpha and beta diversity estimates of
order q at specified sample coverage and measures of dissimilarity. By
default of <code>base = “coverage”</code> and <code>level = NULL</code>,
then this function computes the gamma, alpha, beta diversity, and four
dissimilarity-turnover indices estimates up to one (for q = 1, 2) or up
to the coverage of double the reference sample size (for q = 0).

<table class="gmisc_table" style="border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;">
<thead>
<tr>
<th style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">
Argument
</th>
<th style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">
Description
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">
data
</td>
<td style="text-align: left;">
data can be input as a matrix/data.frame (species by assemblages), or a
list of matrices/data.frames, each matrix represents
species-by-assemblages abundance matrix
</td>
</tr>
<tr>
<td style="text-align: left;">
diversity
</td>
<td style="text-align: left;">
selection of diversity type: ‘TD’ = Taxonomic diversity, ‘PD’ =
Phylogenetic diversity, and ‘FD’ = Functional diversity.
</td>
</tr>
<tr>
<td style="text-align: left;">
level
</td>
<td style="text-align: left;">
a sequence specifying the particular sample coverages (between 0 and 1).
Default is seq(0.5, 1, 0.05).
</td>
</tr>
<tr>
<td style="text-align: left;">
q
</td>
<td style="text-align: left;">
a numerical vector specifying the diversity orders. Default is c(0, 1,
2).
</td>
</tr>
<tr>
<td style="text-align: left;">
nboot
</td>
<td style="text-align: left;">
a positive integer specifying the number of bootstrap replications when
assessing sampling uncertainty and constructing confidence intervals.
Bootstrap replications are generally time consuming. Enter 0 to skip the
bootstrap procedures. Default is 30.
</td>
</tr>
<tr>
<td style="text-align: left;">
conf
</td>
<td style="text-align: left;">
a positive number \< 1 specifying the level of confidence interval.
Default is 0.95.
</td>
</tr>
<tr>
<td style="text-align: left;">
PDtype
</td>
<td style="text-align: left;">
(required only when diversity = ‘PD’), select PD type: PDtype = ‘PD’
(effective total branch length) or PDtype = ‘meanPD’ (effective number
of equally divergent lineages). Default is ‘meanPD’, where meanPD =
PD/tree depth.
</td>
</tr>
<tr>
<td style="text-align: left;">
row.tree
</td>
<td style="text-align: left;">
(required only when diversity = ‘PD’ a phylogenetic tree of row
assemblage in the pooled network row assemblage.
</td>
</tr>
<tr>
<td style="text-align: left;">
col.tree
</td>
<td style="text-align: left;">
(required only when diversity = ‘PD’) a phylogenetic tree of column
assemblage in the pooled network column assemblage.
</td>
</tr>
<tr>
<td style="text-align: left;">
row.distM
</td>
<td style="text-align: left;">
(required only when </code>diversity = ‘FD’) a species pairwise distance
matrix for all species of row assemblage in the pooled network row
assemblage.
</td>
</tr>
<tr>
<td style="text-align: left;">
col.distM
</td>
<td style="text-align: left;">
(required only when diversity = ‘FD’) a species pairwise distance matrix
for all species of column assemblage in the pooled network column
assemblage.
</td>
</tr>
<tr>
<td style="text-align: left;">
FDtype
</td>
<td style="text-align: left;">
(required only when diversity = ‘FD’), select FD type: FDtype =
‘tau_values’ for FD under specified threshold values, or FDtype = ‘AUC’
(area under the curve of tau-profile) for an overall FD which integrates
all threshold values between zero and one. Default is ‘AUC’.
</td>
</tr>
<tr>
<td style="border-bottom: 2px solid grey; text-align: left;">
FDtau
</td>
<td style="border-bottom: 2px solid grey; text-align: left;">
(required only when diversity = ‘FD’ and FDtype = ‘tau_values’), a
numerical vector between 0 and 1 specifying tau values (threshold
levels). If NULL (default), then threshold is set to be the mean
distance between any two individuals randomly selected from the pooled
assemblage (i.e., quadratic entropy).
</td>
</tr>
</tbody>
</table>

the `iNEXTbeta.link()` function returns the `"iNEXTbeta.link"` object
including seven data frames for each datasets:

-   gamma
-   alpha
-   beta
-   C ( Sorensen-type non-overlap )
-   U ( Jaccard-type non-overlap )
-   V ( Sorensen-type turnover )
-   S ( Jaccard-type turnover )

Here only show the first six rows in each table output:

``` r
# Taxonomic diversity
Abundance_TD = iNEXTbeta.link(data = beetles, diversity = 'TD', level = NULL, q = c(0, 1, 2))
Abundance_TD
```

    $gamma
        Dataset Order.q    SC   Size  Gamma      Method  s.e.    LCL    UCL Diversity
    1 Dataset_1       0 0.500 20.941 13.978 Rarefaction 1.058 11.905 16.052        TD
    2 Dataset_1       0 0.525 25.102 16.017 Rarefaction 1.222 13.622 18.411        TD
    3 Dataset_1       0 0.550 30.426 18.489 Rarefaction 1.397 15.752 21.226        TD
    4 Dataset_1       0 0.575 37.171 21.450 Rarefaction 1.574 18.365 24.535        TD
    5 Dataset_1       0 0.600 45.591 24.932 Rarefaction 1.751 21.500 28.363        TD
    6 Dataset_1       0 0.625 56.005 28.975 Rarefaction 1.932 25.189 32.762        TD

    $alpha
        Dataset Order.q    SC   Size  Alpha      Method  s.e.    LCL    UCL Diversity
    1 Dataset_1       0 0.500 25.949  8.432 Rarefaction 0.411  7.626  9.238        TD
    2 Dataset_1       0 0.525 32.249  9.972 Rarefaction 0.477  9.037 10.907        TD
    3 Dataset_1       0 0.550 40.247 11.826 Rarefaction 0.541 10.765 12.887        TD
    4 Dataset_1       0 0.575 50.044 13.974 Rarefaction 0.606 12.786 15.161        TD
    5 Dataset_1       0 0.600 61.802 16.403 Rarefaction 0.677 15.076 17.730        TD
    6 Dataset_1       0 0.625 75.870 19.132 Rarefaction 0.760 17.643 20.621        TD

    $beta
        Dataset Order.q    SC   Size  Beta      Method  s.e.   LCL   UCL Diversity
    1 Dataset_1       0 0.500 25.949 1.658 Rarefaction 0.028 1.602 1.714        TD
    2 Dataset_1       0 0.525 32.249 1.606 Rarefaction 0.029 1.550 1.662        TD
    3 Dataset_1       0 0.550 40.247 1.563 Rarefaction 0.028 1.509 1.617        TD
    4 Dataset_1       0 0.575 50.044 1.535 Rarefaction 0.026 1.484 1.586        TD
    5 Dataset_1       0 0.600 61.802 1.520 Rarefaction 0.025 1.472 1.568        TD
    6 Dataset_1       0 0.625 75.870 1.514 Rarefaction 0.023 1.469 1.560        TD

    $`1-C`
        Dataset Order.q    SC   Size Dissimilarity      Method  s.e.   LCL   UCL Diversity
    1 Dataset_1       0 0.500 25.949         0.658 Rarefaction 0.028 0.602 0.714        TD
    2 Dataset_1       0 0.525 32.249         0.606 Rarefaction 0.029 0.550 0.662        TD
    3 Dataset_1       0 0.550 40.247         0.563 Rarefaction 0.028 0.509 0.617        TD
    4 Dataset_1       0 0.575 50.044         0.535 Rarefaction 0.026 0.484 0.586        TD
    5 Dataset_1       0 0.600 61.802         0.520 Rarefaction 0.025 0.472 0.568        TD
    6 Dataset_1       0 0.625 75.870         0.514 Rarefaction 0.023 0.469 0.560        TD

    $`1-U`
        Dataset Order.q    SC   Size Dissimilarity      Method  s.e.   LCL   UCL Diversity
    1 Dataset_1       0 0.500 25.949         0.794 Rarefaction 0.020 0.754 0.834        TD
    2 Dataset_1       0 0.525 32.249         0.755 Rarefaction 0.022 0.712 0.798        TD
    3 Dataset_1       0 0.550 40.247         0.721 Rarefaction 0.022 0.677 0.764        TD
    4 Dataset_1       0 0.575 50.044         0.697 Rarefaction 0.022 0.654 0.740        TD
    5 Dataset_1       0 0.600 61.802         0.684 Rarefaction 0.021 0.643 0.726        TD
    6 Dataset_1       0 0.625 75.870         0.679 Rarefaction 0.020 0.640 0.719        TD

    $`1-V`
        Dataset Order.q    SC   Size Dissimilarity      Method  s.e.   LCL   UCL Diversity
    1 Dataset_1       0 0.500 25.949         0.658 Rarefaction 0.028 0.602 0.714        TD
    2 Dataset_1       0 0.525 32.249         0.606 Rarefaction 0.029 0.550 0.662        TD
    3 Dataset_1       0 0.550 40.247         0.563 Rarefaction 0.028 0.509 0.617        TD
    4 Dataset_1       0 0.575 50.044         0.535 Rarefaction 0.026 0.484 0.586        TD
    5 Dataset_1       0 0.600 61.802         0.520 Rarefaction 0.025 0.472 0.568        TD
    6 Dataset_1       0 0.625 75.870         0.514 Rarefaction 0.023 0.469 0.560        TD

    $`1-S`
        Dataset Order.q    SC   Size Dissimilarity      Method  s.e.   LCL   UCL Diversity
    1 Dataset_1       0 0.500 25.949         0.794 Rarefaction 0.020 0.754 0.834        TD
    2 Dataset_1       0 0.525 32.249         0.755 Rarefaction 0.022 0.712 0.798        TD
    3 Dataset_1       0 0.550 40.247         0.721 Rarefaction 0.022 0.677 0.764        TD
    4 Dataset_1       0 0.575 50.044         0.697 Rarefaction 0.022 0.654 0.740        TD
    5 Dataset_1       0 0.600 61.802         0.684 Rarefaction 0.021 0.643 0.726        TD
    6 Dataset_1       0 0.625 75.870         0.679 Rarefaction 0.020 0.640 0.719        TD

The output contains seven data frames: `gamma`, `alpha`, `beta`, `C`,
`U`, `V`, `S`. For each data frame, it includes the diversity estimate
(`Estimate`), the diversity order (`Order.q`), `Method` (Rarefaction,
Observed, or Extrapolation, depending on whether the size `m` is less
than, equal to, or greater than the reference sample size), the sample
coverage estimate (`SC`), the sample size (`Size`), the standard error
from bootstrap replications (`s.e.`), the lower and upper confidence
limits of diversity (`LCL`, `UCL`), and the name of dataset (`Dataset`).
These diversity estimates with confidence intervals are used for
plotting the diversity curve.

### GRAPHIC DISPLAYS FUNCTION: ggiNEXTbeta.link()

The function `ggiNEXTbeta.link()`, which extends `ggplot2` to the
`"iNEXT.link"` object with default arguments, is described as follows:

<table class="gmisc_table" style="border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;">
<thead>
<tr>
<th style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">
Argument
</th>
<th style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">
Description
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">
output
</td>
<td style="text-align: left;">
the output of <code>iNEXTbeta.link</code>.
</td>
</tr>
<tr>
<td style="text-align: left;">
type
</td>
<td style="text-align: left;">
selection of plot type : type = ‘B’ for plotting the gamma, alpha, and
beta diversity; type = ‘D’ for plotting 4 turnover dissimilarities.
</td>
</tr>
<tr>
<td style="border-bottom: 2px solid grey; text-align: left;">
scale
</td>
<td style="border-bottom: 2px solid grey; text-align: left;">
Are scales shared across all facets (the default, ‘fixed’), or do they
vary across rows (‘free_x’), columns (‘free_y’), or both rows and
columns (‘free’)? Default is ‘free’.
</td>
</tr>
</tbody>
</table>

The `ggiNEXTbeta.link()` function is a wrapper around the `ggplot2`
package to create a R/E curve using a single line of code. The resulting
object is of class `"ggplot"`, so it can be manipulated using the
`ggplot2` tools. Users can visualize the output of beta diversity or
four dissimilarities by setting the parameter <code>**type**</code>:

``` r
ggiNEXTbeta.link(Abundance_TD, type = 'B')
```

<img src="README/README-unnamed-chunk-35-1.png" width="576" style="display: block; margin: auto;" />

``` r
ggiNEXTbeta.link(Abundance_TD, type = 'D')
```

<img src="README/README-unnamed-chunk-36-1.png" width="576" style="display: block; margin: auto;" />

### HOW TO CITE iNEXT.link

If you publish your work based on the results from the `iNEXT.link`
package, you should make references to the following methodology paper:

-   Chiu, C-H., Chao, A., Vogel, S., Kriegel, P. and Thorn, S. (2023).
    Quantifying and estimating ecological network diversity based on
    incomplete sampling data. Philosophical Transactions of the Royal
    Society B, 378: 20220183. <https://doi.org/10.1098/rstb.2022.0183>

### License

The iNEXT.link package is licensed under the GPLv3. To help refine
`iNEXT.link`, your comments or feedback would be welcome (please send
them to Anne Chao or report an issue on the iNEXT.link github
[iNEXT.link_github](https://github.com/AnneChao/iNEXT.link).

### References

-   Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H.,
    Dornelas, M and Magurran, A. E. (2021). Measuring temporal change in
    alpha diversity: a framework integrating taxonomic, phylogenetic and
    functional diversity and the iNEXT.3D standardization. Methods in
    Ecology and Evolution, 12, 1926-1940.

-   Chao, A. & Jost, L. (2012) Coverage-based rarefaction and
    extrapolation: standardizing samples by completeness rather than
    size. Ecology, 93, 2533M2547

-   Chao, A., Y. Kubota, D. Zelen??, C.-H. Chiu, C.-F. Li, B.
    Kusumoto, M. Yasuhara, S. Thorn, C.-L. Wei, M. J. Costello,
    and R. K. Colwell (2020). Quantifying sample completeness and
    comparing diversities among assemblages. Ecological Research, 35,
    292-314.

-   Chiu, C-H., Chao, A., Vogel, S., Kriegel, P. and Thorn, S. (2023).
    Quantifying and estimating ecological network diversity based on
    incomplete sampling data. Philosophical Transactions of the Royal
    Society B, 378: 20220183. <https://doi.org/10.1098/rstb.2022.0183>
