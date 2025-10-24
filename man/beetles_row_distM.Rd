\name{beetles_row_distM}
\alias{beetles_row_distM}
\docType{data}
\title{
Species pairwise distance matrix for species in the dataset \code{beetles} and \code{beetles_plotA}.
}
\description{
This dataset includes species pairwise symmetric distance matrix for 238 species listed in the 
  \code{beetles} and \code{beetles_plotA} dataset. Each element in the matrix represents a Gower distance between two interaction computed from species traits.
}
\usage{data("beetles_row_distM")}
\format{
  \code{beetles_row_distM} is a 264 x 264 (species by species) symmetric data.frame. Each element of the data.frame is between zero (for species pairs with identical traits) and one. \cr
  $ Tachyta_nana          : num  0 0.203 0.243 0.257 0.262 ... \cr
  $ Plegaderus_vulneratus : num  0.203 0 0.222 0.236 0.24 ... \cr
  $ Plegaderus_caesus     : num  0.2435 0.2217 0 0.0536 0.0182 ... \cr
  ....... \cr
}
\references{
Vogel, S., et al. (2020). Optimizing enrichment of deadwood for biodiversity by varying sun   exposure and tree species: an experimental approach. Journal of Applied Ecology 57(10): 2075-2085.
}

\keyword{datasets}
