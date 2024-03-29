% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MainFun.R
\name{ggObsAsy.link}
\alias{ggObsAsy.link}
\title{ggplot for Asymptotic Network diversity}
\usage{
ggObsAsy.link(output)
}
\arguments{
\item{output}{the output of the functions \code{ObsAsy.link} .\cr}
}
\value{
a figure of asymptotic or empirical (observed) diversity in q-profile.\cr\cr
}
\description{
\code{ggObsAsy.link} Plots q-profile based on the output of \code{ObsAsy.link} using the ggplot2 package.\cr
}
\examples{
## Taxonomic diversity
data(beetles)
output1 = ObsAsy.link(data = beetles, diversity = 'TD', q = seq(0, 2, 0.2))
ggObsAsy.link(output1)


## Phylogenetic diversity
data(beetles)
data(beetles_col_tree)
output2 = ObsAsy.link(data = beetles, diversity = 'PD', q = seq(0, 2, 0.2), col.tree = beetles_col_tree)
ggObsAsy.link(output2)


## Functional diversity under single threshold
data(beetles)
data(beetles_col_distM)
output3 = ObsAsy.link(data = beetles, diversity = 'FD', q = seq(0, 2, 0.2), col.distM = beetles_col_distM, FDtype = "tau_values")
ggObsAsy.link(output3)


## Functional diversity with thresholds integrating from 0 to 1
data(beetles)
data(beetles_col_distM)
output4 = ObsAsy.link(data = beetles, diversity = 'FD', q = seq(0, 2, 0.25), col.distM = beetles_col_distM, FDtype = "AUC", nboot = 0)
ggObsAsy.link(output4)
}
