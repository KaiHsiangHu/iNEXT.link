<!-- README.md is generated from README.Rmd. Please edit that file -->

# iNEXT.link (R package)

<h5 align="right">
Latest version: 2023-04-10
</h5>
<font color="394CAE">
<h3 color="394CAE" style="font-weight: bold">
Introduction to iNEXT.link (R package): Excerpt from iNEXT.link
UserGuide
</h3>
</font> <br>
<h5>
<b>Anne Chao, K.S. Hu, K.W. Chen, C.G. Lo, S.Y. Wang</b> <br><br>
<i>Institute of Statistics, National Tsing Hua University, Hsin-Chu,
Taiwan 30043</i>
</h5>

<br> `iNEXT.link` is an R package that extends the concepts of iNEXT.3D,
iNEXT.4step and iNEXT.beta3D to ecological networks. In this document,
we provide a brief overview of `iNEXT.link` and its functionalities.
Detailed information about the `iNEXT.link` functions can be found in
the `iNEXT.link` manual, which is also available on
[Github](https://github.com/AnneChao). For users without an R
background, an online version(“<https://chao.shinyapps.io/iNEXT_link/>”)
of `iNEXT.link` is also available.

`iNEXT.link` is primarily designed to calculate and analyze various
measures of diversity in ecological networks. Specifically, the package
calculates three Hill numbers of order q (species richness, Shannon
diversity, and Simpson diversity), as well as phylogenetic and
functional diversity levels.

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

## HOW TO RUN iNEXT.link:

The `iNEXT.link` package can be downloaded from Anne Chao’s
[iNEXT.link_github](https://github.com/AnneChao/iNEXT.link) using the
following commands. For a first-time installation, additional
visualization extension packages (`ggplot2`) from CRAN and (`iNEXT.3D`),
(`iNEXT.4steps`), and (`iNEXT.beta3D`) from Anne Chao’s github must be
installed and loaded.

``` r
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

## For Single community:

-   **iNEXT.link** : Computes taxonomic, phylogenetic, and functional
    diversity estimates and sample coverage estimates.

-   **DataInfo.link** : exhibits basic data information

-   **estimateD.link** : computes species diversity with a particular
    user-specified level of sample size or sample coverage.

-   **AO.link**:compute asymptotic (or observed) diversity of order q.

-   **Completeness.link** : Calculates estimated sample completeness
    with order q.

-   **Spec.link** : Computes specialization estimation (observed) with
    order q.

## For Multi-community:

-   **iNEXTbeta.link** : Computing standardized gamma, alpha, beta
    diversity, and four dissimilarity-turnover indices for three
    dimensions: taxonomic, phylogenetic and functional diversity at
    specified sample coverage or sample size.

## Visualizing Results:

-   **ggCompleteness.link** : Visualizing the output from the function
    `Completeness.link`

-   **ggSpec.link** : Visualizing the output from the function
    `Spec.link`

-   **ggAO.link** : Visualizing the output from the function `AO.link`

-   **ggiNEXT.link** : Visualizing the output from the function
    `iNEXT.link`

-   **ggiNEXTbeta.link** : Visualizing the output from the function
    `iNEXTbeta.link`

First, we load data from `iNEXT.link`:

### SINGLE COMMUNITY FUNCTION: iNEXT.link()

We first describe the main function `iNEXT.link()` with default
arguments:

<br> iNEXT.link(data,diversity = “TD”, q = c(0, 1, 2), size = NULL, nT =
NULL, endpoint = NULL, knots = 40, conf = 0.95, nboot = 30, row.tree =
NULL, col.tree = NULL, PDtype = “meanPD”, row.distM = NULL, col.distM =
NULL, FDtype = “AUC”, FDtau = NULL) <br>

The arguments of this function are briefly described below, and will be
explained in more details by illustrative examples in later text.This
main function computes diversity estimates of order q, the sample
coverage estimates and related statistics for K (if `knots = K`)
evenly-spaced knots (sample sizes) between size 1 and the `endpoint`,
where the endpoint is described below. Each knot represents a particular
sample size for which diversity estimates will be calculated. By
default, endpoint = double the reference sample size (total sample size
for abundance data). For example, if `endpoint = 10`, `knot = 4`,
diversity estimates will be computed for a sequence of samples with
sizes (1, 4, 7, 10).

This function returns an “iNEXT.link” object which can be further used
to make plots using the function ggiNEXT.link() to be described below.

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

Individual-based abundance data : Input data for each assemblage/site
include samples species abundances in an empirical sample of n
individuals (“reference sample”). When dealing with N assemblages, the
input data consists of N lists of species abundances

## RAREFACTION/EXTRAPOLATION VIA EXAMPLES

The data sets (tree-beetles interaction data ) are included in
iNEXT.link package. The experiment took place in the Steigerwald forest
in Germany, where deadwood objects from six tree species were exposed in
open and closed habitats. Saproxilic beetles were sampled using stem
emergence traps and classified according to their functional traits.
Data from four years were pooled for each plot and habitat, and pairwise
distances were computed from the Gower distance. Here, the demonstration
only uses data from plot A. For these data, the following commands
display the sample species abundances and run the `iNEXT.link()`
function for three types of diversty (`"TD"`, `"PD"`, `"FD"` with
threshold dmean, `"AUC"` which integate FD from threshold 0-1) in
`q = 0`.

If one diversity class required, then `iNEXT.link()` function returns
including: `$Info` for summarizing data information; `$iNextEst` for
showing diversity estimates along with related statistics for a series
of rarefied and extrapolated samples; and `$AsyEst` for showing
asymptotic diversity estimates along with related statistics, otherwise,
returns lists which length equal to the number of diversity class
required and also named by diversity class. Among each list include
three data frames:`$Info`, `$iNextEst` and `$AsyEst`.

`$Info`, as shown below, returns basic data information. It can also be
presented using function `DataInfo()` (for “TD”), `PDInfo()` (for “PD”)
and `AUCInfo()` (for “FD”).

Because the three kinds of diversity output are similar therefore the
demo shows only TD.

``` r
linkoutTD = iNEXT.link(data = beetles,diversity = 'TD', q = c(0,1,2),nboot = 30)
linkoutTD$DataInfo
  Networks    n S.obs(row) S.obs(col) Links.obs Connectance Coverage  f1 f2 f3 f4 f5 f6 f7 f8 f9 f10
1   Closed  816          6         83       178      0.3574   0.8800  98 31 15  3  3  5  0  2  2   1
2     Open 1932          6         88       206      0.3902   0.9431 110 33 15  8  7  5  3  2  3   2
```

Diversity estimates and related statistics are computed for these 40
knots, which locates the reference sample at the midpoint of the
selected knots. If the argument se=TRUE, then the bootstrap method is
applied to obtain the 95% confidence intervals for each diversity and
sample coverage estimates.

For the sample size corresponding to each knot, the first data frame of
list `$iNextEst` (as shown below for “size_based”) under each diversity
class includes the sample size (`m`, i.e., each of the 40 knots), the
method (`interpolated`, `observed`, or `extrapolated`, depending on
whether the size `m` is less than, equal to, or greater than the
reference sample size), the diversity order, the diversity estimate of
order q (`qD`, `qPD`, `qFD` and `qAUC`), the 95% lower and upper
confidence limits of diversity (`qD.LCL` with`qD.UCL`, `qPD.LCL` with
`qPD.UCL`, `qFD.LCL` with `qFD.UCL` and `qAUC.LCL` with `qAUC.UCL`), and
the sample coverage estimate (`SC`) along with the 95% lower and upper
confidence limits of sample coverage (`SC.LCL`, `SC.UCL`). These sample
coverage estimates with confidence intervals are used for plotting the
sample completeness curve.

``` r
head(linkoutTD$iNextEst$size_based)
# A tibble: 6 x 10
  Assemblage     m Method      Order.q    qD qD.LCL qD.UCL     SC SC.LCL SC.UCL
  <chr>      <dbl> <chr>         <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
1 Closed         1 Rarefaction       0   1      1      1   0.0295 0.0252 0.0337
2 Closed        43 Rarefaction       0  28.0   26.8   29.2 0.543  0.509  0.577 
3 Closed        86 Rarefaction       0  44.8   42.3   47.3 0.661  0.631  0.691 
4 Closed       129 Rarefaction       0  58.2   54.4   61.9 0.713  0.686  0.740 
5 Closed       172 Rarefaction       0  69.8   65.0   74.6 0.745  0.720  0.770 
6 Closed       215 Rarefaction       0  80.2   74.5   86.0 0.767  0.744  0.790 
```

The second data frame of list `$iNextEst` (as shown below for
“coverage_based”) under each class includes real sample coverage (“SC”),
sample size (`m`, i.e., each of the 40 knots), the method
(`Rarefaction`, `Observed`, or `Extrapolation`, depending on whether the
size `m` is less than, equal to, or greater than the reference sample
size),the diversity estimate of order q (qD, qPD, qFD and qAUC)
conditioning on “SC”, the 95% lower and upper confidence limits of
diversity (`qD.LCL` with `qD.UCL`, `qPD.LCL` with `qPD.UCL`, `qFD.LCL`
with `qFD.UCL` and `qAUC.LCL` with `qAUC.UCL`). These sample coverage
estimates with confidence intervals are used for plotting the
coverage-based R/E curves.

``` r
head(linkoutTD$iNextEst$coverage_based)
# A tibble: 6 x 8
  Assemblage     SC      m Method      Order.q    qD qD.LCL qD.UCL
  <chr>       <dbl>  <dbl> <chr>         <dbl> <dbl>  <dbl>  <dbl>
1 Closed     0.0295   1.00 Rarefaction       0  1.00  0.937   1.06
2 Closed     0.543   43.0  Rarefaction       0 28.0  23.9    32.2 
3 Closed     0.661   86.0  Rarefaction       0 44.8  37.2    52.4 
4 Closed     0.713  129.   Rarefaction       0 58.2  47.6    68.7 
5 Closed     0.745  172.   Rarefaction       0 69.8  56.8    82.7 
6 Closed     0.767  215.   Rarefaction       0 80.2  65.3    95.2 
```

`$AsyEst` lists the observed diversity, asymptotic estimates, estimated
bootstrap s.e. and 95% confidence intervals for Hill numbers with q =
0(`Species richness`), 1(`Shannon diversity`), and
2(`Simpson diversity`). The estimated asymptotes and the observed
diversity are calculated via the functions `AO.link()` . The output for
the dunes data is shown below. All row and column variables are
self-explanatory.

``` r
head(linkoutTD$AsyEst)
  Assemblage         Diversity   Observed  Estimator       s.e.        LCL        UCL
1     Closed  Species richness 178.000000 332.713393 34.4345148 265.222985 400.203802
2     Closed Shannon diversity  64.690053  80.316088  4.5779439  71.343483  89.288693
3     Closed Simpson diversity  32.627205  33.944467  1.9604777  30.102001  37.786933
4       Open  Species richness 206.000000 389.238440 36.1090429 318.466017 460.010864
5       Open Shannon diversity  23.311453  25.911799  1.3610219  23.244245  28.579353
6       Open Simpson diversity   8.059978   8.089554  0.3307473   7.441301   8.737807
```

## GRAPHIC DISPLAYS: FUNCTION ggiNEXT.link()

The function `ggiNEXT.link()`, which extends `ggplot2` with default
arguments, is described as follows:

<br> ggiNEXT.link(outcome, type = 1:3,facet.var = “Assemblage”,
color.var = “Order.q”)  
<br>

Here `outcome` is the object of `iNEXT.link()`’s output. Three types of
curves are allowed for different diversity classes:

1.  Sample-size-based R/E curve (`type=1`): see Figs. 1a and 2a in the
    main text. This curve plots diversity estimates with confidence
    intervals as a function of sample size up to double the reference
    sample size, by default, or a user-specified `endpoint`.

2.  Sample completeness curve (`type=2`) with confidence intervals: see
    Figs. 1b and 2b in the main text. This curve plots the sample
    coverage with respect to sample size for the same range described in
    (1).

3.  Coverage-based R/E curve (`type=3`): see Figs. 1c and 2c in the main
    text. This curve plots the diversity estimates with confidence
    intervals as a function of sample coverage up to the maximum
    coverage obtained from the maximum size described in (1).

The argument `facet.var=("Order.q", "Assemblage")` is used to create a
separate plot for each value of the specified variable. For example, the
following code displays a separate plot (in Figs 1a and 1c) for each
value of the diversity order q. The `ggiNEXT.link()` function is a
wrapper around `ggplot2` package to create a R/E curve using a single
line of code. The resulting object is of class `"ggplot"`, so can be
manipulated using the `ggplot2` tools.

The argument `facet.var="Assemblage"` in `ggiNEXT.link` function creates
a separate plot for each assembalge, therefore the different Order.q
will seperated by different colours as shown below:

``` r
# Sample-size-based R/E curves, separating by "assemblage""
ggiNEXT.link(linkoutTD, type = 1, facet.var = "Assemblage")
[[1]]
```

<img src="README/README-unnamed-chunk-8-1.png" width="672" style="display: block; margin: auto;" />

The argument `facet.var="Order.q"` in `ggiNEXT.link` function creates a
separate plot for each order, therefore six assemblages will be
seperated by different colours as shown below:

``` r
# Sample-size-based R/E curves, separating by "Order.q"
ggiNEXT.link(linkoutTD, type = 1, facet.var = "Order.q")
[[1]]
```

<img src="README/README-unnamed-chunk-9-1.png" width="672" style="display: block; margin: auto;" />

The following commands return the sample completeness curve in which
different colors are used for the six assemblages. Since the sample
completeness curve are same for differnet class of diversity,
`ggiNEXT.link` returns only one plot:

``` r
ggiNEXT.link(linkoutTD, type = 2, facet.var = "Order.q", color.var="Assemblage")
[[1]]
```

<img src="README/README-unnamed-chunk-10-1.png" width="672" style="display: block; margin: auto;" />

The following commands return the coverage-based R/E sampling curves in
which different colors are used for the six assemblages
(`facet.var="Assemblage"`) and for three orders(`facet.var="Order.q"`)

``` r
ggiNEXT.link(linkoutTD, type = 3, facet.var="Assemblage")
[[1]]
```

<img src="README/README-unnamed-chunk-11-1.png" width="672" style="display: block; margin: auto;" />

``` r
ggiNEXT.link(linkoutTD, type = 3, facet.var="Order.q")
[[1]]
```

<img src="README/README-unnamed-chunk-12-1.png" width="672" style="display: block; margin: auto;" />

## DATA INFORMATION FUNCTION: DataInfo.link()

We can supply the function

<br> DataInfo.link(data, diversity = “TD”, row.tree = NULL, col.tree =
NULL, row.distM = NULL, col.distM = NULL) <br>

to compute three type diversity(‘TD’,‘PD’,‘FD’) data information, which
including sample size, observed species richness, sample coverage
estimate, and the first ten abundance frequency counts, and so on.

``` r
DataInfo.link(beetles, diversity = 'TD')
  Networks    n S.obs(row) S.obs(col) Links.obs Connectance Coverage  f1 f2 f3 f4 f5 f6 f7 f8 f9 f10
1   Closed  816          6         83       178      0.3574   0.8800  98 31 15  3  3  5  0  2  2   1
2     Open 1932          6         88       206      0.3902   0.9431 110 33 15  8  7  5  3  2  3   2
```

``` r
DataInfo.link(beetles, diversity = 'PD', col.tree = beetles_col_tree)
  Networks    n S.obs(row) S.obs(col) Links.obs Connectance f1* f2*       g1       g2 PD.obs mean_T
1   Closed  816          6         83       178      0.3574  98  31 11508.99 3873.052    924    285
2     Open 1932          6         88       206      0.3902 110  33 13918.45 5172.410   1014    285
```

``` r
DataInfo.link(beetles, diversity = 'FD', col.distM = beetles_col_distM)
  Networks    n S.obs(row) S.obs(col) Links.obs Connectance  f1 f2 a1' a2' threshold
1   Closed  816          6         83       178      0.3574  98 31   0   0  8.930741
2     Open 1932          6         88       206      0.3902 110 33   0   0  8.021019
```

## POINT ESTIMATION FUNCTION: estimateD.link()

We also supply the function

<br> estimateD.link(data, diversity = “TD”, q = c(0, 1, 2),  
base = “coverage”, level = NULL, nboot = 50, conf = 0.95, PDtype =
“meanPD”, row.tree = NULL, col.tree = NULL, row.distM = NULL, col.distM
= NULL, FDtype = “AUC”, FDtau = NULL) <br>

to compute three type diversity(‘TD’,‘PD’,‘FD’) estimates with q = 0, 1,
2 for any particular level of sample size (`base="size"`) or any
specified level of sample coverage (`base="coverage"`) for abundance
data . If `level=NULL`, this function computes the diversity estimates
for the minimum sample size/coverage among all assemblages.

For example, the following command returns the taxonomic diversity
(‘TD’) with a specified level of sample coverage of 70% for the
tree-beetles interaction data. For some assemblages, this coverage value
corresponds to the rarefaction part whereas the others correspond to
extrapolation, as indicated in the method of the output.

``` r
estimateD.link(beetles, diversity = 'TD', q = c(0,1,2), base = "coverage",level = 0.7)
  Assemblage  SC         m      Method Order.q        qD      s.e.    qD.LCL    qD.UCL
1     Closed 0.7 115.47429 Rarefaction       0 54.183847 3.7391230 46.855301 61.512394
2     Closed 0.7 115.47429 Rarefaction       1 37.823241 2.2323244 33.447965 42.198516
3     Closed 0.7 115.47429 Rarefaction       2 26.409739 1.3704514 23.723703 29.095774
4       Open 0.7  41.05942 Rarefaction       0 17.426867 1.4349989 14.614321 20.239413
5       Open 0.7  41.05942 Rarefaction       1 10.562168 0.6261249  9.334986 11.789351
6       Open 0.7  41.05942 Rarefaction       2  6.898403 0.2830484  6.343639  7.453168
```

## ASYMPTOTIC AND OBSERVED DIVERSITY FUNCTION: AO.link

<br> AO.link(data, diversity = “TD”, q = seq(0, 2, 0.2), nboot = 30,
conf = 0.95, method = c(“Asymptotic”, “Observed”), row.tree = NULL,
col.tree = NULL, PDtype = “meanPD”, row.distM = NULL, col.distM = NULL,
FDtype = “AUC”, FDtau = NULL) <br>

The function `AO.link()` can compute three type
diversity(‘TD’,‘PD’,‘FD’),which including empirical diversity and
asymptotic diversity. For any specified level of q can be compute.

For example, the following command returns an empirical taxonomic
diversity(‘TD’) and asymptotic taxonomic diversity(‘TD’) for dunes data,
along with its confidence interval, for a specified q level from 0 to 2.

``` r
out1 <- AO.link(beetles, diversity = 'TD', q = seq(0, 2, 0.2), method = c("Asymptotic", "Observed"),
                nboot = 5,conf = 0.95)

out1
   Order.q         qD       s.e.     qD.LCL     qD.UCL Network     Method
1      0.0 332.713393 27.6842710 308.568044 366.135164  Closed Asymptotic
2      0.2 265.977501 19.1639055 248.770603 286.512935  Closed Asymptotic
3      0.4 203.167245 12.2660226 192.055647 218.501068  Closed Asymptotic
4      0.6 149.540048  7.5456657 142.737347 160.265071  Closed Asymptotic
5      0.8 108.592546  4.9000398 104.688123 115.983548  Closed Asymptotic
6      1.0  80.316088  3.5691820  76.912147  85.381432  Closed Asymptotic
7      1.2  61.983531  2.8233555  58.450432  65.472350  Closed Asymptotic
8      1.4  50.295815  2.3375278  47.007803  52.776058  Closed Asymptotic
9      1.6  42.699384  2.0055407  39.721021  44.541638  Closed Asymptotic
10     1.8  37.568079  1.7797798  34.851500  38.999000  Closed Asymptotic
11     2.0  33.944467  1.6289016  31.418236  35.101545  Closed Asymptotic
12     0.0 389.238440 23.9500574 360.262312 418.294163    Open Asymptotic
13     0.2 263.223968 14.2916069 245.538231 278.399898    Open Asymptotic
14     0.4 159.174040  7.4655280 150.177567 165.700538    Open Asymptotic
15     0.6  86.647814  3.6648832  82.973284  91.157042    Open Asymptotic
16     0.8  45.571274  1.8440972  43.804805  48.145928    Open Asymptotic
17     1.0  25.911799  0.9564169  24.979600  27.222580    Open Asymptotic
18     1.2  16.963237  0.5277983  16.417248  17.635319    Open Asymptotic
19     1.4  12.643810  0.3233401  12.286290  13.013092    Open Asymptotic
20     1.6  10.333946  0.2228161  10.075356  10.557546    Open Asymptotic
21     1.8   8.967960  0.1713854   8.765575   9.125483    Open Asymptotic
22     2.0   8.089554  0.1444350   7.921235   8.249060    Open Asymptotic
23     0.0 178.000000  4.6690470 172.100000 182.600000  Closed  Empirical
24     0.2 149.140436  3.9496851 143.683988 153.083960  Closed  Empirical
25     0.4 122.284123  3.2189158 117.563672 125.444575  Closed  Empirical
26     0.6  98.784711  2.5527819  94.970468 101.184304  Closed  Empirical
27     0.8  79.539008  2.0204157  76.608498  81.313093  Closed  Empirical
28     1.0  64.690053  1.6433614  62.479137  66.444433  Closed  Empirical
29     1.2  53.711702  1.3937261  52.048203  55.458055  Closed  Empirical
30     1.4  45.761569  1.2287478  44.472832  47.440153  Closed  Empirical
31     1.6  40.008527  1.1169462  38.957960  41.604942  Closed  Empirical
32     1.8  35.788959  1.0407265  34.881725  37.311582  Closed  Empirical
33     2.0  32.627205  0.9903577  31.802314  34.092468  Closed  Empirical
34     0.0 206.000000  8.9888820 192.200000 213.300000    Open  Empirical
35     0.2 147.911160  6.7523332 137.630135 153.511141    Open  Empirical
36     0.4  98.276240  4.6859813  91.225503 102.391872    Open  Empirical
37     0.6  60.906771  3.0291282  56.587374  63.834306    Open  Empirical
38     0.8  36.832763  1.9021038  34.421597  38.881011    Open  Empirical
39     1.0  23.311453  1.2297123  22.002732  24.759636    Open  Empirical
40     1.2  16.195155  0.8511934  15.263719  17.261646    Open  Empirical
41     1.4  12.394607  0.6365222  11.665307  13.224986    Open  Empirical
42     1.6  10.237828  0.5093079   9.645089  10.919826    Open  Empirical
43     1.8   8.920980  0.4296588   8.419724   9.506050    Open  Empirical
44     2.0   8.059978  0.3771210   7.621307   8.578857    Open  Empirical
```

## GRAPHIC DISPLAYS FUNCTION: ggAO.link()

Plots q-profile, time-profile, and tau-profile based on the outcome of
AO.link using the ggplot2 package.

The function `ggAO.link()`, which extends `ggplot2` with default
arguments, is described as follows:

<br> ggAO.link(outcome) <br>

Here `outcome` is the object of `AO.link`’s output .

``` r
# q profile curve""
ggAO.link(out1)
```

<img src="README/README-unnamed-chunk-18-1.png" width="672" style="display: block; margin: auto;" />

## SINGLE COMMUNITY FUNCTION: Completeness.link()

Function `Completeness.link()` provides a easy way to compute estimated
sample completeness with order q. It has default arguments: <br>
Completeness.link(data, q = seq(0, 2, 0.2), nboot = 30, conf = 0.95)
<br>

## GRAPHIC DISPLAYS FUNCTION: ggCompleteness.link()

We also provides a realized function `ggCompleteness.link` to plot the
output from `Completeness.link()`:

<br>ggCompleteness.link(output) <br>

We use data beetles to calculate sample completeness and plot it.

``` r
out1 <- Completeness.link(data = beetles)
ggCompleteness.link(out1)
```

<img src="README/README-unnamed-chunk-19-1.png" width="70%" style="display: block; margin: auto;" />

## SINGLE COMMUNITY FUNCTION: Spec.link()

We describe the main function `Spec.link()` with default arguments:

<br> Spec.link(data, q = seq(0, 2, 0.2), method = “Estimated”, nboot =
30, conf = 0.95, E.class = c(1:5), C = NULL) <br>

## GRAPHIC DISPLAYS FUNCTION: ggSpec.link()

We provide a function `ggSpec.link()` to plot the output from
`Spec.link()`.

<br> ggSpec.link(output) <br>

There are an example for funciton `Spec.link` and function
`ggSpec.link`.

``` r
out1 <- Spec.link(data = beetles)
ggSpec.link(out1)
```

<img src="README/README-unnamed-chunk-20-1.png" width="672" style="display: block; margin: auto;" />

### Multi-community

## MULTI-COMMUNITY FUNCTION: iNEXTbeta.link()

``` r
iNEXTbeta.link(data, diversity = "TD", level = seq(0.5, 1, 0.05), 
               q = c(0, 1, 2), nboot = 20, conf = 0.95, 
               PDtype = "meanPD", row.tree = NULL, col.tree = NULL, row.distM = NULL, 
               col.distM = NULL, FDtype = "AUC", FDtau = NULL, FDcut_number = 30) 
```

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
including seven data frames for each regions:

-   gamma
-   alpha
-   beta
-   C ( Sorensen-type non-overlap )
-   U ( Jaccard-type non-overlap )
-   V ( Sorensen-type turnover )
-   S ( Jaccard-type turnover )

``` r
# Taxonomic diversity
Abundance_TD = iNEXTbeta.link(data = beetles, diversity = 'TD', q = c(0, 1, 2))
Abundance_TD
```

    $gamma
      Estimate Order.q      Method  SC   Size  s.e.    LCL    UCL   Region diversity
    1   13.978       0 Rarefaction 0.5 20.941 0.743 12.522 15.435 Region_1        TD
    2   11.552       1 Rarefaction 0.5 20.941 0.559 10.456 12.648 Region_1        TD
    3    9.072       2 Rarefaction 0.5 20.941 0.383  8.321  9.824 Region_1        TD
    4   24.932       0 Rarefaction 0.6 45.591 1.466 22.059 27.804 Region_1        TD
    5   17.460       1 Rarefaction 0.6 45.591 0.853 15.788 19.132 Region_1        TD
    6   11.614       2 Rarefaction 0.6 45.591 0.449 10.734 12.495 Region_1        TD

    $alpha
      Estimate Order.q      Method  SC   Size  s.e.    LCL    UCL   Region diversity
    1    8.432       0 Rarefaction 0.5 25.949 0.456  7.538  9.326 Region_1        TD
    2    6.706       1 Rarefaction 0.5 25.949 0.333  6.054  7.359 Region_1        TD
    3    5.011       2 Rarefaction 0.5 25.949 0.218  4.583  5.440 Region_1        TD
    4   16.403       0 Rarefaction 0.6 61.802 0.804 14.827 17.978 Region_1        TD
    5   10.562       1 Rarefaction 0.6 61.802 0.459  9.663 11.460 Region_1        TD
    6    6.342       2 Rarefaction 0.6 61.802 0.242  5.868  6.816 Region_1        TD

    $beta
      Estimate Order.q      Method  SC   Size  s.e.   LCL   UCL   Region diversity
    1    1.658       0 Rarefaction 0.5 25.949 0.025 1.608 1.707 Region_1        TD
    2    1.723       1 Rarefaction 0.5 25.949 0.021 1.681 1.764 Region_1        TD
    3    1.810       2 Rarefaction 0.5 25.949 0.015 1.781 1.839 Region_1        TD
    4    1.520       0 Rarefaction 0.6 61.802 0.020 1.481 1.559 Region_1        TD
    5    1.653       1 Rarefaction 0.6 61.802 0.016 1.622 1.684 Region_1        TD
    6    1.831       2 Rarefaction 0.6 61.802 0.011 1.810 1.852 Region_1        TD

    $C
      Estimate Order.q      Method  SC   Size  s.e.   LCL   UCL   Region diversity
    1    0.658       0 Rarefaction 0.5 25.949 0.025 0.608 0.707 Region_1        TD
    2    0.785       1 Rarefaction 0.5 25.949 0.018 0.750 0.819 Region_1        TD
    3    0.895       2 Rarefaction 0.5 25.949 0.009 0.877 0.913 Region_1        TD
    4    0.520       0 Rarefaction 0.6 61.802 0.020 0.481 0.559 Region_1        TD
    5    0.725       1 Rarefaction 0.6 61.802 0.014 0.698 0.753 Region_1        TD
    6    0.908       2 Rarefaction 0.6 61.802 0.006 0.895 0.921 Region_1        TD

    $U
      Estimate Order.q      Method  SC   Size  s.e.   LCL   UCL   Region diversity
    1    0.794       0 Rarefaction 0.5 25.949 0.019 0.757 0.830 Region_1        TD
    2    0.785       1 Rarefaction 0.5 25.949 0.018 0.750 0.819 Region_1        TD
    3    0.810       2 Rarefaction 0.5 25.949 0.015 0.781 0.839 Region_1        TD
    4    0.684       0 Rarefaction 0.6 61.802 0.017 0.650 0.718 Region_1        TD
    5    0.725       1 Rarefaction 0.6 61.802 0.014 0.698 0.753 Region_1        TD
    6    0.831       2 Rarefaction 0.6 61.802 0.011 0.810 0.852 Region_1        TD

    $V
      Estimate Order.q      Method  SC   Size  s.e.   LCL   UCL   Region diversity
    1    0.658       0 Rarefaction 0.5 25.949 0.025 0.608 0.707 Region_1        TD
    2    0.723       1 Rarefaction 0.5 25.949 0.021 0.681 0.764 Region_1        TD
    3    0.810       2 Rarefaction 0.5 25.949 0.015 0.781 0.839 Region_1        TD
    4    0.520       0 Rarefaction 0.6 61.802 0.020 0.481 0.559 Region_1        TD
    5    0.653       1 Rarefaction 0.6 61.802 0.016 0.622 0.684 Region_1        TD
    6    0.831       2 Rarefaction 0.6 61.802 0.011 0.810 0.852 Region_1        TD

    $S
      Estimate Order.q      Method  SC   Size  s.e.   LCL   UCL   Region diversity
    1    0.794       0 Rarefaction 0.5 25.949 0.019 0.757 0.830 Region_1        TD
    2    0.839       1 Rarefaction 0.5 25.949 0.014 0.811 0.867 Region_1        TD
    3    0.895       2 Rarefaction 0.5 25.949 0.009 0.877 0.913 Region_1        TD
    4    0.684       0 Rarefaction 0.6 61.802 0.017 0.650 0.718 Region_1        TD
    5    0.790       1 Rarefaction 0.6 61.802 0.012 0.767 0.813 Region_1        TD
    6    0.908       2 Rarefaction 0.6 61.802 0.006 0.895 0.921 Region_1        TD

The output contains seven data frames: `gamma`, `alpha`, `beta`, `C`,
`U`, `V`, `S`. For each data frame, it includes the diversity estimate
(`Estimate`), the diversity order (`Order.q`), `Method` (Rarefaction,
Observed, or Extrapolation, depending on whether the size `m` is less
than, equal to, or greater than the reference sample size), the sample
coverage estimate (`SC`), the sample size (`Size`), the standard error
from bootstrap replications (`s.e.`), the 95% lower and upper confidence
limits of diversity (`LCL`, `UCL`), and the name of region (`Region`).
These diversity estimates with confidence intervals are used for
plotting the diversity curve.

## GRAPHIC DISPLAYS FUNCTION: ggiNEXTbeta.link()

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

<img src="README/README-unnamed-chunk-26-1.png" width="576" style="display: block; margin: auto;" />

``` r
ggiNEXTbeta.link(Abundance_TD, type = 'D')
```

<img src="README/README-unnamed-chunk-27-1.png" width="576" style="display: block; margin: auto;" />

### How to cite

-   Chiu, C-H., Chao, A., Vogel, S., Kriegel, P. and Thorn, S. (2023).
    Network-diversity quantification and related statistical estimation:
    drawing on sampling models and methodologies from biodiversity
    research. To appear in Philosophical Transactions of the Royal
    Society B.

### Referance

Chiu, C-H., Chao, A., Vogel, S., Kriegel, P. and Thorn, S. (2023).
Network-diversity quantification and related statistical estimation:
drawing on sampling models and methodologies from biodiversity research.
To appear in Philosophical Transactions of the Royal Society B.
