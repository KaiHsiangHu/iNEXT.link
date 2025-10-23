\name{iNEXT.link-package}
\alias{iNEXT.link-package}
\docType{package}
\title{
Network diversity iNterpolation/Extrapolation Estimation
}
\description{
\code{iNEXT.link} is an R package that extends the concepts of \code{iNEXT.3D} (Chao et al., 2021), \code{iNEXT.4step} (Chao et al., 2020) and \code{iNEXT.beta3D} (Chao et al., 2023) to ecological networks (Chiu et al., 2023). An online version "iNEXT.link Online"
(https://chao.shinyapps.io/iNEXT_link/) is also available for users without an R background.

\code{iNEXT.link} is designed to calculate and analyze various measures of diversity in ecological networks. Specifically, the package calculates three Hill numbers of order q (species richness, Shannon diversity, and Simpson diversity) in taxonomic diversity (TD) level, as well as their generalizations (Hill-Chao numbers) for phylogenetic diversity (PD) and functional diversity (FD) levels. 

For individual ecological networks, \code{iNEXT.link} provides tools for rarefaction and extrapolation (R/E), enabling diversity estimation and confidence intervals under both size-based and coverage-based standardization. Asymptotic 3D network diversity estimates are also available (see Chao et al., 2021). 

In addition, \code{iNEXT.link} offers dissimilarity-turnover curves based on coverage-based R/E for gamma, alpha, and beta diversity, allowing cross-network comparison of diversity patterns (see Chao et al., 2023). 

The package also assesses sample completeness for networks and supports coverage-based rarefaction and extrapolation for estimating specialization using five unevenness-based indices (see Chao et al., 2020).

This package contains thirteen main functions: \cr\cr
** Functions for single network **

1. \code{iNEXT.link} computes rarefaction/extrapolation taxonomic, phylogenetic, and functional diversity estimates and sample coverage estimates.

2. \code{DataInfo.link} exhibits basic data information
    
3. \code{estimateD.link} computes species diversity with a particular user-specified level of sample size or sample coverage.
    
4. \code{ObsAsy.link} computes asymptotic or empirical(observed) diversity of order q.

5. \code{Completeness.link} calculates estimated sample completeness with order q. 

6. \code{Spec.link.est} computes standardized specialization estimation under specified sample coverage with order q.

7. \code{Spec.link.ObsAsy} computes observed or asymptotic standardized specialization  with order q.

** Function for Multi-networks **

8. \code{iNEXTbeta.link} computing standardized gamma, alpha, beta diversity, and four dissimilarity-turnover indices for three dimensions: taxonomic, phylogenetic and functional diversity at specified sample coverage.


** Functions for Visualizing Results **

9. \code{ggCompleteness.link} visualizes the output from the function \code{Completeness.link}

10. \code{ggSpec.link} visualizes the output from the function \code{Spec.link.est} and \code{Spec.link.ObsAsy}

11. \code{ggObsAsy.link} visualizes the output from the function \code{ObsAsy.link}

12. \code{ggiNEXT.link} visualizes the output from the function \code{iNEXT.link}
    
13. \code{ggiNEXTbeta.link} visualizes the output from the function \code{iNEXTbeta.link}


An Online version of \code{iNEXT.link} is also available for users without an R background: \cr
https://chao.shinyapps.io/iNEXT_link/. \cr
In the detailed Online \code{iNEXT.link} User's Guide, we illustrate all the running procedures in an easily
accessible way through numerical examples with proper interpretations of portions of the output.
All the data of those illustrative examples are included in this package. \cr\cr
}
\author{
Anne Chao, Kai-Hsiang Hu \cr

Maintainer: Anne Chao <chao@stat.nthu.edu.tw>
}
\references{
Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H., Dornelas, M and Magurran, A. E. (2021). Measuring temporal change in alpha diversity: a framework integrating taxonomic, phylogenetic and functional diversity and the iNEXT.3D standardization. Methods in Ecology and Evolution, 12, 1926-1940. \cr\cr
Chao, A., Thorn, S., Chiu, C.-H., Moyes, F., Hu, K.-H., Chazdon, R. L., Wu, J., Dornelas, M., Zeleny, D., Colwell, R. K., and Magurran, A. E. (2023). Rarefaction and extrapolation with beta diversity under a framework of Hill numbers: the iNEXT.beta3D standardization. To appear in Ecological Monographs. \cr\cr
Chao, A., Y. Kubota, D. Zeleny, C.-H. Chiu, C.-F. Li, B. Kusumoto, M. Yasuhara, S. Thorn, C.-L. Wei, M. J. Costello, and R. K. Colwell (2020). Quantifying sample completeness and comparing diversities among assemblages. Ecological Research, 35, 292-314. \cr\cr
Chiu, C-H., Chao, A., Vogel, S., Kriegel, P. and Thorn, S. (2023). Quantifying and estimating ecological network diversity based on incomplete sampling data. Philosophical Transactions of the Royal Society B, 378: 20220183. https://doi.org/10.1098/rstb.2022.0183. \cr\cr
}
\keyword{ package }
