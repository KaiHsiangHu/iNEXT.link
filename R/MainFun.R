# DataInfo.link ----------------------
#' Exhibit basic data information
#'
#' \code{DataInfo.link}: exhibit basic data information for three dimensions of biodiversity.
#'
#' @param data  a \code{list} of \code{data.frames}, each \code{data.frames} represents col.species-by-row.species abundance matrix.
#' @param diversity selection of diversity type: \code{'TD'} = 'Taxonomic diversity', \code{'PD'} = 'Phylogenetic diversity', and \code{'FD'} = 'Functional diversity'.
#' @param row.tree (required only when \code{diversity = "PD"}), a phylogenetic tree of row assemblage in the pooled network row assemblage.
#' @param col.tree (required only when \code{diversity = "PD"}), a phylogenetic tree of column assemblage in the pooled network column assemblage.
#' @param row.distM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species of row assemblage in the pooled network row assemblage.
#' @param col.distM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species of column assemblage in the pooled network column assemblage.
#' 
#' @return a data.frame including basic data information.\cr\cr
#' Basic information shared by TD, mean-PD and FD includes Network name (\code{Network}),
#' number of observed interaction events in the reference sample (\code{n}), number of observed species in row assemblage (\code{S.obs(row)}), number of observed species in column assemblage (\code{S.obs(col)}), 
#' number of observed interactions in the reference sample (\code{Link.obs}), the proportion of links between species that are realized in the network matrix (\code{Connectance}), estimator of the sample coverage of the reference sample (\code{SC(n)}).\cr\cr
#' Other additional information is given below.\cr\cr
#' (1) TD: the first ten frequency counts in the reference sample (\code{f1}--\code{f10}).\cr\cr
#' (2) Mean-PD: the number of those singletons and doubletons in the node/branch set (\code{f1*},\code{f2*}) , 
#' the total branch length of those singletons (\code{g1}) and of those doubletons (\code{g2}) in the node/branch set,
#' the observed total branch length in the phylogenetic tree (\code{PD.obs}),  and the product of column tree depth and row tree depth (\code{T1*T2}).\cr\cr
#' (3) FD : the number of singletons and doubletons in the reference sample (\code{f1},\code{f2}), the number of singletons and doubletons in the functional group (\code{a1*},\code{a2*}),
#' and the mean weighted pairwise distance between any two interactions (\code{d_mean}).\cr\cr.\cr\cr
#' 
#' 
#' @examples
#' # Taxonomic network diversity for interaction data
#' data(beetles_plotA)
#' DataInfo.link(data = beetles_plotA, diversity = 'TD')
#'
#'
#' # Phylogenetic network diversity for interaction data
#' data(beetles_plotA)
#' data(beetles_row_tree)
#' DataInfo.link(data = beetles_plotA, diversity = 'PD', row.tree = beetles_row_tree)
#'
#'
#' # Functional network diversity for interaction data
#' data(beetles_plotA)
#' data(beetles_row_distM)
#' DataInfo.link(data = beetles_plotA, diversity = 'FD', row.distM = beetles_row_distM)
#' @export
#' 

# DataInfo.link(data = lapply(beetles_plotA,function(x){data.frame(t(x))}), diversity = 'PD', row.tree = beetles_row_tree)
# DataInfo.link(data = beetles_plotA , diversity = 'PD', row.tree = beetles_row_tree)
# DataInfo.link(data = beetles_plotA , diversity = 'FD', row.distM = beetles_row_distM)


DataInfo.link <- function(data, diversity = 'TD', row.tree = NULL, col.tree = NULL, row.distM = NULL, col.distM = NULL){

  datatype = "abundance"
  
  data_new = data
  
  # if(class(data == list)){
  #   if(names(data[[1]])[1] == names(data_new[[1]])[1]){
  #     row.tree = row.tree
  #     col.tree = col.tree
  #     row.distM = row.distM
  #     col.distM = col.distM
  #   }else{
  #     #change tree
  #     rowtree = row.tree
  #     coltree = col.tree
  #     row.tree = coltree
  #     col.tree = rowtree
  #     #change distM
  #     rowdistM = row.distM
  #     coldistM = col.distM
  #     row.distM = coldistM
  #     col.distM = rowdistM
  #   }
  # }else{
  #   if(colnames(data)[1] == colnames(data_new)[1]){
  #     row.tree = row.tree
  #     col.tree = col.tree
  #     row.distM = row.distM
  #     col.distM = col.distM
  #   }else{
  #     #change tree
  #     rowtree = row.tree
  #     coltree = col.tree
  #     row.tree = coltree
  #     col.tree = rowtree
  #     #change distM
  #     rowdistM = row.distM
  #     coldistM = col.distM
  #     row.distM = coldistM
  #     col.distM = rowdistM
  #   }
  # }
  
  
  if(diversity == 'PD'){

    if(!is.null(row.tree)){row.tree$tip.label = gsub('\\.', '_',row.tree$tip.label)}
    if(!is.null(col.tree)){col.tree$tip.label = gsub('\\.', '_',col.tree$tip.label)}
    
    if(inherits(data_new, "list")){
      table <- lapply(data_new, function(y){datainfphy(data = y, datatype = datatype,
                                                       row.tree = row.tree,col.tree = col.tree)})%>%
        do.call(rbind,.)
      rownames(table) <- names(data_new)
      table = tibble::rownames_to_column(table, var = "Networks")
    }else{
      table <- datainfphy(data = data_new, datatype = datatype,
                          row.tree = row.tree,col.tree = col.tree)
      rownames(table) <- "Network1"
      table = tibble::rownames_to_column(table, var = "Networks")
    }
    
  }else if(diversity == 'TD'){
    
    if(inherits(data_new, "list")){
      table <- lapply(data_new, function(y){datainf(data = y, datatype = datatype)})%>%do.call(rbind,.)
      rownames(table) <- names(data_new)
      table = tibble::rownames_to_column(table, var = "Networks")
    }else{
      table <- datainf(data = data_new, datatype = datatype)
      rownames(table) <- "Network1"
      table = tibble::rownames_to_column(table, var = "Networks")
    }
  }else if(diversity == 'FD'){
    
    if(inherits(data_new, "list")){
      table <- lapply(data_new, function(y){datainffun(data = y, datatype = datatype,
                                                       row.distM = row.distM,col.distM = col.distM)})%>%
        do.call(rbind,.)
      rownames(table) <- names(data_new)
      table = tibble::rownames_to_column(table, var = "Networks")
    }else{
      table <- datainffun(data = data_new, datatype = datatype,
                          row.distM = row.distM,col.distM = col.distM)
      rownames(table) <- "Network1"
      table = tibble::rownames_to_column(table, var = "Networks")
    }

    
  }
  # for(i in 1:length(data)){
  #   if(dim(data_new[[i]])[1] == dim(data[[i]])[1] & dim(data_new[[i]])[2] == dim(data[[i]])[2]){
  #     table[i,] <- table[i,]
  #   }else{
  #     temp <- c(table[i,3], table[i,4])
  #     table[i,3] <- temp[2]
  #     table[i,4] <- temp[1]
  #   }
  # }
  return(table)

}


# Completeness.link ----
#' Sample Completeness main function
#'
#' \code{Completeness.link}: Estimation of Sample Completeness with order q
#'
#' @param data a \code{list} of \code{data.frames}, each \code{data.frames} represents col.species-by-row.species abundance matrix.
#' @param q q a numerical vector specifying the diversity orders. Default is \code{seq(0, 2, 0.2)}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing
#' sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter 0 to skip the bootstrap procedures. Default is \code{30}.
#' @param conf a positive number < \code{1} specifying the level of confidence interval. Default is \code{0.95}.
#' @return a matrix of estimated sample completeness with order q: 
#'         \item{Order.q}{the network completeness order of q.}
#'         \item{Estimate.SC}{the estimated completeness of order q.}
#'         \item{s.e.}{standard error of the estimated network completeness.}
#'         \item{SC.LCL, SC.UCL}{the bootstrap lower and upper confidence limits for the expected network completeness of order q at the specified level (with a default value of \code{0.95}).}
#'         \item{Dataset}{the dataset name.}
#'
#' @examples
#' # Sample Completeness for interaction data
#' data(beetles_plotA)
#' output_com = Completeness.link(beetles_plotA)
#' output_com
#'
#' @references
#' Chao, A., Y. Kubota, D. Zeleny, C.-H. Chiu, C.-F. Li, B. Kusumoto, M. Yasuhara, S. Thorn, C.-L. Wei, M. J. Costello, and R. K. Colwell (2020). Quantifying sample completeness and comparing diversities among assemblages. Ecological Research, 35, 292-314.
#' @export
Completeness.link <- function(data, q = seq(0, 2, 0.2), nboot = 30, conf = 0.95){
  
  datatype = "abundance"
  
  # for(i in 1:length(data)){
  #   if(nrow(data[[i]]) > ncol(data[[i]])){
  #     data[[i]] <- as.data.frame(t(data[[i]]))
  #   }else{
  #     data[[i]] <- data[[i]]
  #   }
  # }
  
  if(inherits(data, c("data.frame", "matrix", "integer"))) data = list(Network1 = data)
  
  data_long <- lapply(data, function(tab){
    as.matrix(tab)%>%c()}
  )
  res = iNEXT.4steps::Completeness(data = data_long, q = q, datatype = datatype, nboot = nboot, conf = conf)
  names(res)[6] <- "Dataset"
  return(res)
}

# ggCompleteness.link -------------------------------------------------------------------
#' ggplot for Sample Completeness
#'
#' \code{ggCompleteness.link} The figure for estimation of Sample Completeness with order q
#'
#' @param output a table generated from Completeness.link function
#' @return a figure of estimated sample completeness with order q
#'
#' @examples
#' # Plot the completeness curve
#' data(beetles_plotA)
#' output_com = Completeness.link(beetles_plotA)
#' ggCompleteness.link(output_com)
#' 
#' @references
#' Chao, A., Y. Kubota, D. Zeleny, C.-H. Chiu, C.-F. Li, B. Kusumoto, M. Yasuhara, S. Thorn, C.-L. Wei, M. J. Costello, and R. K. Colwell (2020). Quantifying sample completeness and comparing diversities among assemblages. Ecological Research, 35, 292-314.
##' @export
ggCompleteness.link <- function(output){
  
  # Check if the number of unique 'Assemblage' is 8 or less
  if (length(unique(output$Dataset)) <= 8){
    cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                      "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  }else{
    # If there are more than 8 assemblages, start with the same predefined color palette
    # Then extend the palette by generating additional colors using the 'ggplotColors' function
    cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#330066", "#CC79A7", "#0072B2", "#D55E00"))
    #cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
    #                   "#330066", "#CC79A7", "red", "blue"))
    cbPalette <- c(cbPalette, ggplotColors(length(unique(output$Dataset))-8))
  }
  
  ggplot(output, aes(x = Order.q, y = Estimate.SC, colour = Dataset)) +
    geom_line(size = 1.2) + scale_colour_manual(values = cbPalette) +
    geom_ribbon(aes(ymin = SC.LCL, ymax = SC.UCL, fill = Dataset),
                alpha = 0.2, linetype = 0) + theme_bw() + scale_fill_manual(values = cbPalette) +
    labs(x = "Order q", y = "Sample completeness") + theme(text = element_text(size = 16)) +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(1.2, "cm"), legend.title = element_blank())
}


# iNEXT.link -------------------------------------------------------------------
#' Interpolation (rarefaction) and extrapolation of network diversity
#' @param data  a \code{list} of \code{data.frames}, each \code{data.frames} represents col.species-by-row.species abundance matrix.
#' @param diversity selection of diversity type: \code{'TD'} = 'Taxonomic diversity', \code{'PD'} = 'Phylogenetic diversity', and \code{'FD'} = 'Functional diversity'.
#' @param q a numerical vector specifying the diversity orders. Default is \code{c(0,1,2)}.
#' @param size an integer vector of sample sizes for which diversity estimates will be computed.
#' If \code{NULL}, then diversity estimates will be calculated for those sample sizes determined by the specified/default \code{endpoint}
#' and \code{knots}.
#' @param endpoint an integer vector specifying the sample size that is the \code{endpoint} for rarefaction/extrapolation.
#' If \code{NULL}, then \code{endpoint} = double reference sample size in each assemblage.
#' @param knots an integer specifying the number of equally-spaced \code{knots} between size 1 and the \code{endpoint}. Default is \code{40}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing
#' sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter 0 to skip the bootstrap procedures. Default is \code{30}.
#' @param conf a positive number < \code{1} specifying the level of confidence interval. Default is \code{0.95}.
#' @param row.tree (required only when \code{diversity = "PD"}), a phylogenetic tree of row assemblage in the pooled network row assemblage.
#' @param col.tree (required only when \code{diversity = "PD"}), a phylogenetic tree of column assemblage in the pooled network column assemblage.
#' @param PDtype (required only when \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"}(effective total branch length) or
#' \code{PDtype = "meanPD"}(effective number of equally divergent lineages).Default is \code{"meanPD"}.
#' @param row.distM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species of row assemblage in the pooled network row assemblage.
#' @param col.distM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species of column assemblage in the pooled network column assemblage.
#' @param FDtype (required only when \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_values"} for FD under specified threshold values, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{"AUC"}.
#' @param FDtau (required only when \code{diversity = "FD"} and \code{FDtype = "tau_values"}), a numerical vector between 0 and 1 specifying tau values (threshold levels). If \code{NULL} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy).
#' @import ape
#' @import tidyverse
#' @import magrittr
#' @import ggplot2
#' @import abind
#' @import phytools
#' @import phyclust
#' @import tidytree
#' @import colorRamps
#' @import iNEXT.3D
#' @import iNEXT.4steps
#' @import iNEXT.beta3D
#' @import future
#' @import future.apply
#' @import ade4
#' @import tidyr
#' @import tibble
#' @import reshape2
#' @import sets
#' @import dplyr
#' @importFrom utils head
#' @importFrom stats qnorm sd optimize quantile rbinom rmultinom
#' @importFrom grDevices hcl
#' @importFrom purrr map map_dfr
#' 
#' @return a list of three objects: \cr\cr
#' (1) \code{$TDInfo} (\code{$PDInfo}, or \code{$FDInfo}) for summarizing data information for q = 0, 1 and 2. Refer to the output of \code{DataInfo.link} for details.
#' (2) \code{$TDiNextEst}(\code{$PDiNextEst}, or \code{$FDiNextEst}) for showing network diversity estimates for rarefied and extrapolated samples along with related statistics. There are two data frames: \code{"$size_based"} and \code{"$coverage_based"}.
#'    In \code{"$size_based"}, the output includes:
#'    \item{Dataset}{the name of datasets.} 
#'    \item{Order.q}{the network diversity order of q.}
#'    \item{m}{the target standardized sample size.}
#'    \item{Method}{Rarefaction, Observed, or Extrapolation, depending on whether the size m is less than, equal to, or greater than the reference sample size.}
#'    \item{qiTD, qiPD, qiFD}{the estimated network diversity of order q for a sample of size m.}
#'    \item{qiTD.LCL, qiPD.LCL, qiFD.LCL and qiTD.UCL, qiPD.UCL, qiFD.UCL}{the bootstrap lower and upper confidence limits for the network diversity of order q at the level specified in the settings (with a default value of 0.95).}
#'    \item{SC}{the estimated network coverage for a sample of size m.}
#'    \item{SC.LCL, SC.UCL}{the bootstrap lower and upper confidence limits for the expected sample coverage at the level specified in the settings (with a default value of 0.95).}
#'    \item{Reftime}{the reference times for PD.}
#'    \item{Type}{\code{"PD"} (phylogenetic network diversity of effective total branch length),\code{"meanPD"} (effective number of equally divergent lineages) for PD.}
#'  Similar output is obtained for \code{"$coverage_based"}. \cr\cr
#' (3) \code{$TDAsyEst}(\code{$PDAsyEst}, or \code{$FDAsyEst}): for showing asymptotic diversity estimates along with related statistics:
#'    \item{Dataset}{the name of datasets.} 
#'    \item{qTD, qPD, qFD}{the network diversity order of q.}
#'    \item{TD_obs, PD_obs, FD_obs}{the observed network diversity.}
#'    \item{TD_asy, PD_asy, FD_asy}{the asymptotic network diversity estimate.}
#'    \item{s.e.}{standard error of asymptotic network diversity.}
#'    \item{qiTD.LCL, qiPD.LCL, qiFD.LCL and qiTD.UCL, qiPD.UCL, qiFD.UCL}{the bootstrap lower and upper confidence limits for asymptotic network diversity at the specified level (with a default value of 0.95).}
#'    \item{Reftime}{the reference times for PD.}
#'    \item{Type}{\code{"PD"} (phylogenetic network diversity of effective total branch length),\code{"meanPD"} (effective number of equally divergent lineages) for PD.}
#' 
#' 
#' @examples
#' # Compute standardized estimates of taxonomic network diversity for interaction data with
#' # order q = 0, 1, 2
#' data(beetles_plotA)
#' output_qiTD = iNEXT.link(beetles_plotA, diversity = 'TD', q = c(0,1,2), nboot = 30)
#' output_qiTD$TDInfo   # showing basic data information.
#' output_qiTD$TDiNextEst   # showing diversity estimates with rarefied and extrapolated.
#' output_qiTD$TDAsyEst   # showing asymptotic diversity estimates.
#'
#'
#' # Compute standardized estimates of phylogenetic network diversity for interaction data with
#' # order q = 0, 1, 2
#' data(beetles_plotA)
#' data(beetles_row_tree)
#' output_qiPD = iNEXT.link(beetles_plotA, diversity = 'PD', q = c(0, 1, 2), nboot = 20,
#'                          row.tree = beetles_row_tree)
#' output_qiPD
#'
#'
#' # Compute standardized estimates of functional network diversity for interaction data with
#' # order q = 0, 1, 2
#' data(beetles_plotA)
#' data(beetles_row_distM)
#' output_qiFD = iNEXT.link(data = beetles_plotA, diversity = 'FD', q = c(0,1,2), nboot = 20, 
#'                          row.distM = beetles_row_distM, FDtype = "AUC")
#' output_qiFD
#'
#' @references
#' Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H., Dornelas, M. and Magurran, A. E. (2021). Measuring temporal change in alpha diversity: a framework integrating taxonomic, phylogenetic and functional diversity and the iNEXT.3D standardization. 
#' \emph{Methods in Ecology and Evolution.}, 12, 1926-1940. \cr\cr
#' @export
#' 


iNEXT.link <- function(data, diversity = 'TD', q = c(0,1,2), size = NULL,
                       endpoint = NULL, knots = 40, conf = 0.95, nboot = 30,
                       row.tree = NULL, col.tree = NULL, PDtype = 'meanPD', row.distM = NULL, col.distM = NULL,
                       FDtype = "AUC", FDtau = NULL
                       ) {
  
  datatype = "abundance"
  nT = NULL
  # data_new <- list()
  # for(i in 1:length(data)){
  #   if(nrow(data[[i]]) > ncol(data[[i]])){
  #     data_new[[i]] <- as.data.frame(t(data[[i]]))
  #     names(data_new)[i] <- names(data)[i]
  #   }else{
  #     data_new[[i]] <- data[[i]]
  #     names(data_new)[i] <- names(data)[i]
  #   }
  # }
  
  if(inherits(data, c("data.frame", "matrix", "integer"))) data = list(Network1 = data)
  data_new = data
  
  
  
  # User interface
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  
  res = list()
  if(diversity == 'TD'){
      ## 1. datainfo
      datainfo = DataInfo.link(data = data, diversity = "TD")
      ## 2. iNterpolation/ Extrapolation
      data_long <- lapply(data_new, function(tab){
        as.matrix(tab)%>%c()}
      )
      if(0 %in% q){
        INEXT_est_0 <- iNEXT.3D::iNEXT3D(data_long, diversity = 'TD', q = q,conf = conf,
                                         nboot = nboot, knots = knots, endpoint = endpoint, size = size)
        INEXT_est_0[[2]]$coverage_based <- INEXT_est_0[[2]]$coverage_based[INEXT_est_0[[2]]$coverage_based$Order.q == 0,]
        
        asy_size <- sapply(1:length(data_new) ,function(i) coverage_to_size(data_new[[i]], 0.999, datatype = "abundance"))
        new_size <- lapply(1:length(data_new) ,function(i) round(seq(2*sum(data_new[[i]]),asy_size[i], length = 5)))
        size_check <-  check.size(data_new, datatype="abundance", size =size, endpoint=endpoint, knots=knots)
        size <- list()
        for(i in 1:length(data_new)){
          size[[i]] <- unique(c(size_check[[i]], new_size[[i]]))
        }
        INEXT_est_q <- iNEXT.3D::iNEXT3D(data_long, diversity = 'TD', q = q[q!=0],conf = conf,
                                         nboot = nboot, knots = NULL, endpoint = endpoint, size = size)
        INEXT_est <- INEXT_est_0
        INEXT_est[[2]]$coverage_based <- rbind(INEXT_est_0[[2]]$coverage_based, INEXT_est_q[[2]]$coverage_based)
        INEXT_est[[2]]$coverage_based <- INEXT_est[["TDiNextEst"]][["coverage_based"]][order(INEXT_est[["TDiNextEst"]][["coverage_based"]]$Assemblage),]
      }else{
        
        asy_size <- sapply(1:length(data_new) ,function(i) coverage_to_size(data_new[[i]], 0.999, datatype = "abundance"))
        new_size <- lapply(1:length(data_new) ,function(i) round(seq(2*sum(data_new[[i]]),asy_size[i], length = 5)))
        size_check <-  check.size(data_new, datatype="abundance", size =size, endpoint=endpoint, knots=knots)
        size <- list()
        for(i in 1:length(data_new)){
          size[[i]] <- unique(c(size_check[[i]], new_size[[i]]))
        }
        INEXT_est <- iNEXT.3D::iNEXT3D(data_long, diversity = 'TD', q = q,conf = conf,
                                       nboot = nboot, knots = knots, endpoint = endpoint, size = size)
      }
      
      res[[1]] = datainfo
      res[[2]] = INEXT_est$TDiNextEst
      res[[3]] = INEXT_est$TDAsyEst
      names(res) = c("TDInfo", "TDiNextEst", "TDAsyEst")
      res$TDiNextEst$size_based <- rename(res$TDiNextEst$size_based, c(qiTD = "qTD",qiTD.LCL = "qTD.LCL", qiTD.UCL = "qTD.UCL"))
      res$TDiNextEst$coverage_based <- rename(res$TDiNextEst$coverage_based, c(qiTD = "qTD",qiTD.LCL = "qTD.LCL", qiTD.UCL = "qTD.UCL"))
      
    }else if(diversity == 'PD'){
      datainfo = DataInfo.link(data = data, diversity = "PD", row.tree = row.tree, col.tree = col.tree)
      if(names(data[[1]])[1] == names(data_new[[1]])[1]){
        row.tree = row.tree
        col.tree = col.tree
        row.distM = row.distM
        col.distM = col.distM
      }else{
        #change tree
        rowtree = row.tree
        coltree = col.tree
        row.tree = coltree
        col.tree = rowtree
        #change distM
        rowdistM = row.distM
        coldistM = col.distM
        row.distM = coldistM
        col.distM = rowdistM
      }
      
      if(!is.null(row.tree)){row.tree$tip.label = gsub('\\.', '_',row.tree$tip.label)}
      if(!is.null(col.tree)){col.tree$tip.label = gsub('\\.', '_',col.tree$tip.label)}
      
      if(0 %in% q){
        INEXT_est_0 <- iNEXTPDlink(data_new, q = q, size = size,
                                   endpoint = endpoint, knots = knots, conf = conf,
                                   nboot = nboot,col.tree = col.tree,row.tree = row.tree, type = PDtype)
        INEXT_est_0[[2]]$coverage_based <- INEXT_est_0[[2]]$coverage_based[INEXT_est_0[[2]]$coverage_based$Order.q == 0,]
        
        asy_size <- sapply(1:length(data_new) ,function(i) coverage_to_size(data_new[[i]], 0.999, datatype = "abundance"))
        new_size <- lapply(1:length(data_new) ,function(i) round(seq(2*sum(data_new[[i]]),asy_size[i], length = 5)))
        size_check <-  check.size(data_new, datatype="abundance", size =size, endpoint=endpoint, knots=knots)
        size <- list()
        for(i in 1:length(data_new)){
          size[[i]] <- unique(c(size_check[[i]], new_size[[i]]))
        }
        INEXT_est_q <- iNEXTPDlink(data_new, q = q[q!=0], size = size,
                                   endpoint = endpoint, knots = knots, conf = conf,
                                   nboot = nboot,col.tree = col.tree,row.tree = row.tree, type = PDtype)
        INEXT_est <- INEXT_est_0
        INEXT_est[[2]]$coverage_based <- rbind(INEXT_est_0[[2]]$coverage_based, INEXT_est_q[[2]]$coverage_based)
        #INEXT_est[[2]]$coverage_based <- INEXT_est[["PDiNextEst"]][["coverage_based"]][order(INEXT_est[["PDiNextEst"]][["coverage_based"]]$Assemblage),]
      }else{
        
        asy_size <- sapply(1:length(data_new) ,function(i) coverage_to_size(data_new[[i]], 0.999, datatype = "abundance"))
        new_size <- lapply(1:length(data_new) ,function(i) round(seq(2*sum(data_new[[i]]),asy_size[i], length = 5)))
        size_check <-  check.size(data_new, datatype="abundance", size =size, endpoint=endpoint, knots=knots)
        size <- list()
        for(i in 1:length(data_new)){
          size[[i]] <- unique(c(size_check[[i]], new_size[[i]]))
        }
        INEXT_est <- iNEXTPDlink(data_new, q = q[q!=0], size = size,
                                 endpoint = endpoint, knots = knots, conf = conf,
                                 nboot = nboot,col.tree = col.tree,row.tree = row.tree, type = PDtype)
      }
      res = INEXT_est
      res[[1]] <- datainfo
      res$PDiNextEst$size_based <- rename(res$PDiNextEst$size_based, c(qiPD = "qPD",qiPD.LCL = "qPD.LCL", qiPD.UCL = "qPD.UCL"))
      res$PDiNextEst$coverage_based <- rename(res$PDiNextEst$coverage_based, c(qiPD = "qPD",qiPD.LCL = "qPD.LCL", qiPD.UCL = "qPD.UCL"))
      names(res[[3]])[c(2:4,6:7)] <- c("qiPD", "PD_obs", "PD_asy", "qiPD.LCL", "qiPD.UCL")
    }else if (diversity == "FD" & FDtype == "tau_values") {
      datainfo = DataInfo.link(data = data, diversity = "FD", row.distM = row.distM, col.distM = col.distM)
      if(names(data[[1]])[1] == names(data_new[[1]])[1]){
        row.tree = row.tree
        col.tree = col.tree
        row.distM = row.distM
        col.distM = col.distM
      }else{
        #change tree
        rowtree = row.tree
        coltree = col.tree
        row.tree = coltree
        col.tree = rowtree
        #change distM
        rowdistM = row.distM
        coldistM = col.distM
        row.distM = coldistM
        col.distM = rowdistM
      }
      if(0 %in% q){
        INEXT_est_0 <- iNEXTlinkFD(data_new, q = q, size = size,
                                   endpoint = endpoint, knots = knots, conf = conf,
                                   nboot = nboot, nT = nT, row.distM = row.distM, col.distM = col.distM, threshold = FDtau)
        INEXT_est_0[[2]]$coverage_based <- INEXT_est_0[[2]]$coverage_based[INEXT_est_0[[2]]$coverage_based$Order.q == 0,]
        
        asy_size <- sapply(1:length(data_new) ,function(i) coverage_to_size(data_new[[i]], 0.999, datatype = "abundance"))
        new_size <- lapply(1:length(data_new) ,function(i) round(seq(2*sum(data_new[[i]]),asy_size[i], length = 5)))
        size_check <-  check.size(data_new, datatype="abundance", size =size, endpoint=endpoint, knots=knots)
        size <- list()
        for(i in 1:length(data_new)){
          size[[i]] <- unique(c(size_check[[i]], new_size[[i]]))
        }
        INEXT_est_q <- iNEXTlinkFD(data_new, q = q[q!=0], size = size,
                                   endpoint = endpoint, knots = knots, conf = conf,
                                   nboot = nboot, nT = nT, row.distM = row.distM, col.distM = col.distM, threshold = FDtau)
        INEXT_est <- INEXT_est_0
        INEXT_est[[2]]$coverage_based <- rbind(INEXT_est_0[[2]]$coverage_based, INEXT_est_q[[2]]$coverage_based)
        INEXT_est[[2]]$coverage_based <- INEXT_est[["FDiNextEst"]][["coverage_based"]][order(INEXT_est[["FDiNextEst"]][["coverage_based"]]$Assemblage),]
      }else{
        
        asy_size <- sapply(1:length(data_new) ,function(i) coverage_to_size(data_new[[i]], 0.999, datatype = "abundance"))
        new_size <- lapply(1:length(data_new) ,function(i) round(seq(2*sum(data_new[[i]]),asy_size[i], length = 5)))
        size_check <-  check.size(data_new, datatype="abundance", size =size, endpoint=endpoint, knots=knots)
        size <- list()
        for(i in 1:length(data_new)){
          size[[i]] <- unique(c(size_check[[i]], new_size[[i]]))
        }
        INEXT_est <- iNEXTlinkAUC(data_new, q = q, size = size,
                                  endpoint = endpoint, knots = knots, conf = conf,
                                  nboot = nboot, nT = nT, row.distM = row.distM, col.distM = col.distM)
      }
      res = INEXT_est
      res[[1]] <- datainfo
      res$FDiNextEst$size_based <- rename(res$FDiNextEst$size_based, c(qiFD = "qFD",qiFD.LCL = "qFD.LCL", qiFD.UCL = "qFD.UCL"))
      res$FDiNextEst$coverage_based <- rename(res$FDiNextEst$coverage_based, c(qiFD = "qFD",qiFD.LCL = "qFD.LCL", qiFD.UCL = "qFD.UCL"))
    }else if (diversity == "FD" & FDtype == "AUC") {
      datainfo = DataInfo.link(data = data, diversity = "FD", row.distM = row.distM, col.distM = col.distM)
      if(names(data[[1]])[1] == names(data_new[[1]])[1]){
        row.tree = row.tree
        col.tree = col.tree
        row.distM = row.distM
        col.distM = col.distM
      }else{
        #change tree
        rowtree = row.tree
        coltree = col.tree
        row.tree = coltree
        col.tree = rowtree
        #change distM
        rowdistM = row.distM
        coldistM = col.distM
        row.distM = coldistM
        col.distM = rowdistM
      }
      if(0 %in% q){
        INEXT_est_0 <- iNEXTlinkAUC(data_new, q = q, size = size,
                                    endpoint = endpoint, knots = knots, conf = conf,
                                    nboot = nboot, nT = nT, row.distM = row.distM, col.distM = col.distM)
        INEXT_est_0[[2]]$coverage_based <- INEXT_est_0[[2]]$coverage_based[INEXT_est_0[[2]]$coverage_based$Order.q == 0,]
        
        asy_size <- sapply(1:length(data_new) ,function(i) coverage_to_size(data_new[[i]], 0.999, datatype = "abundance"))
        new_size <- lapply(1:length(data_new) ,function(i) round(seq(2*sum(data_new[[i]]),asy_size[i], length = 5)))
        size_check <-  check.size(data_new, datatype="abundance", size =size, endpoint=endpoint, knots=knots)
        size <- list()
        for(i in 1:length(data_new)){
          size[[i]] <- unique(c(size_check[[i]], new_size[[i]]))
        }
        INEXT_est_q <- iNEXTlinkAUC(data_new, q = q[q!=0], size = size,
                                    endpoint = endpoint, knots = knots, conf = conf,
                                    nboot = nboot, nT = nT, row.distM = row.distM, col.distM = col.distM)
        INEXT_est <- INEXT_est_0
        INEXT_est[[2]]$coverage_based <- rbind(INEXT_est_0[[2]]$coverage_based, INEXT_est_q[[2]]$coverage_based)
        INEXT_est[[2]]$coverage_based <- INEXT_est[["FDiNextEst"]][["coverage_based"]][order(INEXT_est[["FDiNextEst"]][["coverage_based"]]$Assemblage),]
      }else{
        
        asy_size <- sapply(1:length(data_new) ,function(i) coverage_to_size(data_new[[i]], 0.999, datatype = "abundance"))
        new_size <- lapply(1:length(data_new) ,function(i) round(seq(2*sum(data_new[[i]]),asy_size[i], length = 5)))
        size_check <-  check.size(data_new, datatype="abundance", size =size, endpoint=endpoint, knots=knots)
        size <- list()
        for(i in 1:length(data_new)){
          size[[i]] <- unique(c(size_check[[i]], new_size[[i]]))
        }
        INEXT_est <- iNEXTlinkAUC(data_new, q = q, size = size,
                                  endpoint = endpoint, knots = knots, conf = conf,
                                  nboot = nboot, nT = nT, row.distM = row.distM, col.distM = col.distM)
      }
      res = INEXT_est
      res[[1]] <- datainfo
      res$FDiNextEst$size_based <- rename(res$FDiNextEst$size_based, c(qiFD = "qFD",qiFD.LCL = "qFD.LCL", qiFD.UCL = "qFD.UCL"))
      res$FDiNextEst$coverage_based <- rename(res$FDiNextEst$coverage_based, c(qiFD = "qFD",qiFD.LCL = "qFD.LCL", qiFD.UCL = "qFD.UCL"))
      names(res[[3]])[c(2:4,6:7)] <- c("qiFD", "FD_obs", "FD_asy", "qiFD.LCL", "qiFD.UCL")
     }

  # class(res) = "iNEXT3D"
  names(res[[1]])[1] <- "Dataset"
  names(res[[2]]$size_based)[1] <- "Dataset"
  names(res[[2]]$coverage_based)[1] <- "Dataset"
  names(res[[3]])[1] <- "Dataset"
  return(res)
}





# ggiNEXT.link -------------------------------------------------------------------
#' ggplot2 extension for output from \code{iNEXT.link}
#'
#' \code{ggiNEXT.link} is a \code{ggplot} extension for an \code{iNEXT.link} object to plot sample-size- and coverage-based rarefaction/extrapolation sampling curves along with a bridging sample completeness curve.
#' @param output an \code{iNEXT.link} object computed by \code{iNEXT.link}.
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).            
#' @param facet.var create a separate plot for each value of a specified variable: 
#'  no separation (\code{facet.var = "None"}); 
#'  a separate plot for each diversity order (\code{facet.var = "Order.q"}); 
#'  a separate plot for each assemblage (\code{facet.var = "Assemblage"}); 
#'  a separate plot for each combination of diversity order and assemblage (\code{facet.var = "Both"}).              
#' @param color.var create curves in different colors for values of a specified variable:
#'  all curves are in the same color (\code{color.var = "None"}); 
#'  use different colors for diversity orders (\code{color.var = "Order.q"}); 
#'  use different colors for assemblages/sites (\code{color.var = "Assemblage"}); 
#'  use different colors for combinations of diversity orders and assemblage (\code{color.var = "Both"}).  
#' @return a \code{ggplot2} object for sample-size-based rarefaction/extrapolation curve (\code{type = 1}), sample completeness curve (\code{type = 2}), and coverage-based rarefaction/extrapolation curve (\code{type = 3}).
#' 
#' @examples
#' # Plot three types of curves of taxonomic network diversity with facet.var = "Assemblage"
#' # for interaction data with order q = 0, 1, 2
#' data(beetles_plotA)
#' output_qiTD = iNEXT.link(beetles_plotA, diversity = 'TD', q = c(0,1,2), nboot = 30)
#' ggiNEXT.link(output_qiTD, facet.var = "Assemblage")
#'
#'
#' # Plot two types (1 and 3) of curves of phylogenetic network diversity 
#' # for interaction data with order q = 0, 1, 2
#' data(beetles_plotA)
#' data(beetles_row_tree)
#' output_qiPD = iNEXT.link(beetles_plotA, diversity = 'PD', q = c(0, 1, 2), nboot = 20, 
#'                          row.tree = beetles_row_tree)
#' ggiNEXT.link(output_qiPD, type = c(1, 3))
#'
#'
#' # Plot three types of curves of functional network diversity
#' # for interaction data with order q = 0, 1, 2
#' data(beetles_plotA)
#' data(beetles_row_distM)
#' output_qiFD = iNEXT.link(data = beetles_plotA, diversity = 'FD', q = c(0,1,2), nboot = 0, 
#'                          row.distM = beetles_row_distM, FDtype = "AUC")
#' ggiNEXT.link(output_qiFD)
#' @export

ggiNEXT.link <- function(output, type = c(1,2,3), facet.var = "Assemblage", color.var = "Order.q"){
  
  names(output[[1]])[1] <- "Assemblage"
  names(output[[2]]$size_based)[1] <- "Assemblage"
  names(output[[2]]$coverage_based)[1] <- "Assemblage"
  names(output[[3]])[1] <- "Assemblage"
  
  if (names(output)[1] == 'TDInfo') {
    diversity = 'TD'
    output[[2]]$size_based <- rename(output[[2]]$size_based, c(qTD = "qiTD", qTD.LCL = "qiTD.LCL", qTD.UCL = "qiTD.UCL"))
    output[[2]]$coverage_based <- rename(output[[2]]$coverage_based, c(qTD = "qiTD", qTD.LCL = "qiTD.LCL", qTD.UCL = "qiTD.UCL"))
  } else if (names(output)[1] == "PDInfo") {
    diversity = 'PD'
    output[[2]]$size_based <- rename(output[[2]]$size_based, c(qPD = "qiPD", qPD.LCL = "qiPD.LCL", qPD.UCL = "qiPD.UCL"))
    output[[2]]$coverage_based <- rename(output[[2]]$coverage_based, c(qPD = "qiPD", qPD.LCL = "qiPD.LCL", qPD.UCL = "qiPD.UCL"))
  } else if (names(output)[1] %in% "FDInfo") {
    diversity = 'FD'
    output[[2]]$size_based <- rename(output[[2]]$size_based, c(qFD = "qiFD", qFD.LCL = "qiFD.LCL", qFD.UCL = "qiFD.UCL"))
    output[[2]]$coverage_based <- rename(output[[2]]$coverage_based, c(qFD = "qiFD", qFD.LCL = "qiFD.LCL", qFD.UCL = "qiFD.UCL"))
  }
  
  if(diversity == 'TD'){
    res = iNEXT.3D::ggiNEXT3D(output, type = c(1,2,3), facet.var = facet.var)
    res[[1]] = res[[1]]+ylab("Taxonomic network diversity")+xlab("Sample size")
    res[[2]] = res[[2]]+xlab("Sample size")
    res[[3]] = res[[3]]+ylab("Taxonomic network diversity")
    out = list()
    j = 1
    for(i in type){

      out[[j]] = res[[i]]
      j = j+1
    }
    out
  }else if(diversity == 'PD'){
    res = iNEXT.3D::ggiNEXT3D(output, type = c(1,2,3), facet.var = facet.var)
    res[[1]] = res[[1]]+ylab("Phylogenetic network diversity")+xlab("Sample size")
    res[[2]] = res[[2]]+xlab("Sample size")
    res[[3]] = res[[3]]+ylab("Phylogenetic network diversity")
    out = list()
    j = 1
    for(i in type){

      out[[j]] = res[[i]]
      j = j+1
    }
    out

  }else if(diversity == 'FD'){
    res = iNEXT.3D::ggiNEXT3D(output, type = c(1,2,3), facet.var = facet.var)
    res[[1]] = res[[1]]+ylab("Functional network diversity")+xlab("Sample size")
    res[[2]] = res[[2]]+xlab("Sample size")
    res[[3]] = res[[3]]+ylab("Functional network diversity")
    out = list()
    j = 1
    for(i in type){

      out[[j]] = res[[i]]
      j = j+1
    }
    out
  }
  
}


# ObsAsy.link -------------------------------------------------------------------
#' Asymptotic and Observed diversity q profile
#'
#' \code{ObsAsy.link}: The asymptotic (or observed) diversity of order q
#'
#' @param data a \code{list} of \code{data.frames}, each \code{data.frames} represents col.species-by-row.species abundance matrix.
#' @param diversity selection of diversity type: \code{'TD'} = Taxonomic network diversity, \code{'PD'} = Phylogenetic network diversity, and \code{'FD'} = Functional network diversity.
#' @param q a numerical vector specifying the diversity orders. Default is \code{seq(0, 2, by = 0.2)}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param method Select \code{'Asymptotic'} or \code{'Observed'}.
#' @param row.tree (required argument for \code{diversity = "PD"}), a phylogenetic tree in Newick format for row species in the row assemblage. 
#' @param col.tree (required argument for \code{diversity = "PD"}), a phylogenetic tree in Newick format for column species in the column assemblage. 
#' @param PDtype (argument only for \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"} (effective total branch length) or \code{PDtype = "meanPD"} (effective number of equally divergent lineages). Default is \code{"meanPD"}, where \code{meanPD = PD/tree depth}.
#' @param row.distM (required argument for \code{diversity = "FD"}), a species pairwise distance matrix for row species in the row assemblage. 
#' @param col.distM (required argument for \code{diversity = "FD"}), a species pairwise distance matrix for column species in the column assemblage. 
#' @param FDtype (argument only for \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_values"} for qiFD under specified threshold values, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall qiFD which integrates all threshold values between zero and one. Default is \code{"AUC"}.  
#' @param FDtau (argument only for \code{diversity = "FD"} and \code{FDtype = "tau_values"}), a numerical vector between 0 and 1 specifying tau values (threshold levels). If \code{NULL} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy).
#' @return a table of diversity table including the following arguments.
#' \item{Dataset}{the name of datasets.}
#' \item{Order.q}{the network diversity order of q.}
#' \item{qiTD, qiPD, qiFD}{the asymptotic/observed network network diversity estimates of order q.}
#' \item{s.e.}{standard error of the estimated network diversity.}
#' \item{qiTD.LCL, qiPD.LCL, qiFD.LCL and qiTD.UCL, qiPD.UCL, qiFD.UCL}{the bootstrap lower and upper confidence limits for the network diversity of order q at the specified level (with a default value of \code{0.95}).}
#' \item{Method}{\code{"Asymptotic"} means asymptotic network diversity and \code{"Observed"} means observed network diversity.}
#' \item{Reftime}{the reference times for qiPD.}
#' \item{Type}{\code{"PD"} (phylogenetic network diversity of effective total branch length),\code{"meanPD"} (effective number of equally divergent lineages) for qiPD.}
#'
#' @examples
#' # Compute the observed and asymptotic taxonomic network diversity 
#' # for interaction data with order q between 0 and 2 
#' # (in increments of 0.2 by default)
#' data(beetles_plotA)
#' output_ObsAsy_qiTD = ObsAsy.link(data = beetles_plotA, diversity = 'TD', q = seq(0, 2, 0.2))
#' output_ObsAsy_qiTD
#' 
#'
#' # Compute the observed and asymptotic phylogenetic network diversity
#' # for interaction data with order q between 0 and 2 
#' # (in increments of 0.2 by default)
#' data(beetles_plotA)
#' data(beetles_row_tree)
#' output_ObsAsy_qiPD = ObsAsy.link(data = beetles_plotA, diversity = 'PD', q = seq(0, 2, 0.2), 
#'                                  row.tree = beetles_row_tree, nboot = 10)
#' output_ObsAsy_qiPD
#'
#'
#' # Compute the observed and asymptotic functional network diversity
#' # for interaction data with order q between 0 and 2 and FDtype = "AUC"
#' # (in increments of 0.25 by default)
#' data(beetles_plotA)
#' data(beetles_row_distM)
#' output_ObsAsy_qiFD = ObsAsy.link(data = beetles_plotA, diversity = 'FD', q = seq(0, 2, 0.25), 
#'                                  row.distM = beetles_row_distM, FDtype = "AUC", nboot = 10)
#' output_ObsAsy_qiFD
#' 
#' @export
ObsAsy.link <- function(data, diversity = 'TD', q = seq(0, 2, 0.2), nboot = 30, conf = 0.95, method = c("Asymptotic", "Observed"),
                        row.tree = NULL, col.tree = NULL, PDtype = "meanPD", row.distM = NULL, col.distM = NULL, FDtype = "AUC", FDtau = NULL){
  
  datatype = "abundance"
  if(inherits(data, c("data.frame", "matrix", "integer"))) data = list(Network1 = data)
  
  if ( !(diversity %in% c('TD', 'PD', 'FD')) )
    stop("Please select one of below diversity: 'TD', 'PD', 'FD'", call. = FALSE)
  
  
  # for(i in 1:length(data)){
  #   if(nrow(data[[i]]) > ncol(data[[i]])){
  #     data[[i]] <- as.data.frame(t(data[[i]]))
  #   }else{
  #     data[[i]] <- data[[i]]
  #   }
  # }
  
  if (diversity == 'TD') {
    
    if (sum(method == 'Asymptotic') == length(method))
      NetDiv <- AsylinkTD(data, diversity = 'TD', q = q, datatype = datatype, nboot = nboot, conf = conf) else if (sum(method == 'Observed') == length(method))
        
        NetDiv <- ObslinkTD(data, diversity = 'TD', q = q, datatype = datatype, nboot = nboot, conf = conf) else if (sum(method == c('Asymptotic', 'Observed')) == length(method))
          
          NetDiv = rbind(AsylinkTD(data, diversity = 'TD', q = q, datatype = datatype, nboot = nboot, conf = conf),
                         ObslinkTD(data, diversity = 'TD', q = q, datatype = datatype, nboot = nboot, conf = conf))
        NetDiv <- rename(NetDiv, c(qiTD = "qTD", qiTD.LCL = "qTD.LCL", qiTD.UCL = "qTD.UCL"))
  }
  
  if (diversity == 'PD') {
    
    if (sum(method == 'Asymptotic') == length(method))
      NetDiv = AsylinkPD(data = data,q = q,B = nboot,row.tree = row.tree,
                         col.tree = col.tree,conf = conf,PDtype = PDtype) else if (sum(method == 'Observed') == length(method))
                           
                           NetDiv = ObslinkPD(data = data,q = q,B = nboot,row.tree = row.tree,
                                              col.tree = col.tree,conf = conf,PDtype = PDtype) else if (sum(method == c('Asymptotic', 'Observed')) == length(method))
                                                
                                                NetDiv = rbind(AsylinkPD(data = data,q = q,B = nboot,row.tree = row.tree,
                                                                         col.tree = col.tree,conf = conf,PDtype = PDtype),
                                                               ObslinkPD(data = data,q = q,B = nboot,row.tree = row.tree,
                                                                         col.tree = col.tree,conf = conf,PDtype = PDtype))
                                              
                                              NetDiv <- rename(NetDiv, c(qiPD = "qPD", qiPD.LCL = "qPD.LCL", qiPD.UCL = "qPD.UCL"))
                                              
  }
  
  if (diversity == 'FD' & FDtype == 'tau_values') {
    
    if (sum(method == 'Asymptotic') == length(method))
      NetDiv = AsylinkFD(data = data, q = q, datatype = datatype, nboot = nboot, conf = conf,
                         row.distM = row.distM, col.distM = col.distM, threshold = FDtau) else if (sum(method == 'Observed') == length(method))
                           
                           NetDiv = ObslinkFD(data = data, q = q, datatype = datatype, nboot = nboot, conf = conf,
                                              row.distM = row.distM, col.distM = col.distM, threshold = FDtau) else if (sum(method == c('Asymptotic', 'Observed')) == length(method))
                                                
                                                NetDiv = rbind(AsylinkFD(data = data, q = q, datatype = datatype, nboot = nboot, conf = conf,
                                                                         row.distM = row.distM, col.distM = col.distM, threshold = FDtau),
                                                               ObslinkFD(data = data, q = q, datatype = datatype, nboot = nboot, conf = conf,
                                                                         row.distM = row.distM, col.distM = col.distM, threshold = FDtau))
                                              
                                              NetDiv <- rename(NetDiv, c(qiFD = "qFD", qiFD.LCL = "qFD.LCL", qiFD.UCL = "qFD.UCL"))
                                              
  }
  
  if (diversity == 'FD' & FDtype == 'AUC') {
    
    if (sum(method == 'Asymptotic') == length(method))
      NetDiv = AsylinkAUC(data = data, q = q, datatype = datatype, nboot = nboot, conf = conf,
                          row.distM = row.distM, col.distM = col.distM) else if (sum(method == 'Observed') == length(method))
                            
                            NetDiv = ObslinkAUC(data = data, q = q, datatype = datatype, nboot = nboot, conf = conf,
                                                row.distM = row.distM, col.distM = col.distM) else if (sum(method == c('Asymptotic', 'Observed')) == length(method))
                                                  
                                                  NetDiv = rbind(AsylinkAUC(data = data, q = q, datatype = datatype, nboot = nboot, conf = conf,
                                                                            row.distM = row.distM, col.distM = col.distM),
                                                                 ObslinkAUC(data = data, q = q, datatype = datatype, nboot = nboot, conf = conf,
                                                                            row.distM = row.distM, col.distM = col.distM))
                                                
                                                NetDiv <- rename(NetDiv, c(qiFD = "qFD", qiFD.LCL = "qFD.LCL", qiFD.UCL = "qFD.UCL"))
                                                
  }


  return(NetDiv)

}


# ggObsAsy.link -------------------------------------------------------------------
#' ggplot for Asymptotic Network diversity
#'
#' \code{ggObsAsy.link} Plots q-profile based on the output of \code{ObsAsy.link} using the ggplot2 package.\cr
#'
#' @param output the output of the functions \code{ObsAsy.link} .\cr
#' @return a figure of asymptotic or empirical (observed) diversity in q-profile.\cr\cr
#'
#' @examples
#' # Plot q-profile of taxonomic network diversity for interaction data
#' # with order q between 0 and 2 (in increments of 0.2 by default).
#' data(beetles_plotA)
#' output_ObsAsy_qiTD = ObsAsy.link(data = beetles_plotA, diversity = 'TD', q = seq(0, 2, 0.2))
#' ggObsAsy.link(output_ObsAsy_qiTD)
#'
#'
#' # Plot q-profile of phylogenetic network diversity for interaction data
#' # with order q between 0 and 2 (in increments of 0.2 by default).
#' data(beetles_plotA)
#' data(beetles_row_tree)
#' output_ObsAsy_qiPD = ObsAsy.link(data = beetles_plotA, diversity = 'PD', q = seq(0, 2, 0.2), 
#'                                  row.tree = beetles_row_tree, nboot = 10)
#' ggObsAsy.link(output_ObsAsy_qiPD)
#'
#'
#' # Plot q-profile of functional network diversity for interaction data
#' # with order q between 0 and 2 (in increments of 0.2 by default)
#' # under tau values from 0 to 1
#' data(beetles_plotA)
#' data(beetles_row_distM)
#' output_ObsAsy_qiFD = ObsAsy.link(data = beetles_plotA, diversity = 'FD', q = seq(0, 2, 0.25), 
#'                                  row.distM = beetles_row_distM, FDtype = "AUC", nboot = 10)
#' ggObsAsy.link(output_ObsAsy_qiFD)
#'

#' @export
ggObsAsy.link <- function(output){
  
  if (colnames(output)[3] == 'qiTD') {
    diversity = 'TD'
    output <- rename(output, c(qTD = "qiTD", qTD.LCL = "qiTD.LCL", qTD.UCL = "qiTD.UCL"))
  } else if (colnames(output)[3] == 'qiPD') {
    diversity = 'PD'
    output <- rename(output, c(qPD = "qiPD", qPD.LCL = "qiPD.LCL", qPD.UCL = "qiPD.UCL"))
  } else if (colnames(output)[3] == 'qiFD') {
    diversity = 'FD'
    output <- rename(output, c(qFD = "qiFD", qFD.LCL = "qiFD.LCL", qFD.UCL = "qiFD.UCL"))
  }
  
  if(diversity == 'TD'){
    
    names(output)[names(output) == 'Network'] = 'Assemblage'
    
    iNEXT.3D::ggObsAsy3D(output) + ylab('Taxonomic network diversity')
    
  }else if(diversity == 'PD') {
    
    names(output)[names(output) == 'Network'] = 'Assemblage'
    
    iNEXT.3D::ggObsAsy3D(output) + 
      facet_grid(. ~ .) + 
      ylab('Phylogenetic network diversity')
    
  }else if(diversity == 'FD'){
    
    iNEXT.3D::ggObsAsy3D(output) + ylab("Functional network diversity")
    
  }

}



# estimateD.link  -------------------------------------------------------------------
#' Compute 3D with a particular of sample size/coverage
#'
#' \code{estimateD.link} computes species diversity (Hill numbers with q = 0, 1 and 2) with a particular user-specified level of sample size or sample coverage.
#'
#' @param data a \code{matrix}, \code{data.frame} (species by assemblages), or \code{list} of species abundance/incidence raw data.\cr
#' @param diversity a choice of three-level diversity: 'TD' = 'Taxonomic', 'PD' = 'Phylogenetic', and 'FD' = 'Functional' under certain threshold.
#' @param q a numerical vector of the order of Hill number. Default is \code{seq(0, 2, 0.2)}.
#' @param base comparison base: sample-size-based (\code{base = "size"}) or coverage-based \cr (\code{base = "coverage"}).
#' @param level a sequence specifying the particular sample sizes or sample coverages(between 0 and 1).
#' If \code{base="size"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample size among all sites extrapolated to double reference sizes.
#' If \code{base="coverage"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample coverage among all sites extrapolated to double reference sizes.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing
#' sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is \code{50}.
#' @param conf a positive number < 1 specifying the level of confidence interval, default is \code{0.95}.
#' @param row.tree phylogenetic tree of row assemblage in interaction matrix.
#' @param col.tree phylogenetic tree of column assemblage in interaction matrix.
#' @param PDtype (required only when \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"}(effective total branch length) or
#' \code{PDtype = "meanPD"}(effective number of equally divergent lineages).Default is \code{"meanPD"}.
#' @param row.distM (required only when \code{diversity = "FD"}), a row species pairwise distance matrix for all row species of row assemblage in interaction matrix.
#' @param col.distM (required only when \code{diversity = "FD"}), a column species pairwise distance matrix for all column species of column assemblage in interaction matrix.
#' @param FDtype (required only when \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_values"} for FD under specified threshold values, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{"AUC"}.
#' @param FDtau (required only when \code{diversity = "FD"} and \code{FDtype = "tau_values"}), a numerical vector between 0 and 1 specifying tau values (threshold levels). If \code{NULL} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy).
#' @return a \code{data.frame} of diversity table including the following arguments:
#' \item{Dataset}{the name of datasets.}
#' \item{Order.q}{the diversity order of q.}
#' \item{SC}{the target standardized coverage value.}
#' \item{m}{the corresponding sample size for the standardized completeness.}
#' \item{Method}{Rarefaction, Observed, or Extrapolation, depending on whether the target coverage is less than, equal to, or greater than the coverage of the reference sample.}

#' \item{qiTD, qiPD, qiFD}{the estimated network diversity of order q for the target coverage value The estimate for complete coverage (or \code{level = 1}) represents the estimated asymptotic diversity.} 
#' \item{s.e.}{standard error of network diversity estimate.}
#' \item{qiTD.LCL, qiPD.LCL, qiFD.LCL and qiTD.UCL, qiPD.UCL, qiFD.UCL}{the bootstrap lower and upper confidence limits for the network diversity of order q at the specified level (with a default value of \code{0.95}).}
#' \item{Reftime}{reference times for PD.}
#' \item{Type}{\code{"PD"} (phylogenetic network diversity of effective total branch length), \code{"meanPD"} (effective number of equally divergent lineages).}
#'
#' @examples
#' \donttest{
#' # Taxonomic network diversity for interaction data with two target coverages(93% and 97%)
#' data(beetles_plotA)
#' output_est_qiTD <- estimateD.link(beetles_plotA, diversity = 'TD', q = c(0,1,2),
#'                                   base = "coverage", level = c(0.93, 0.97))
#' output_est_qiTD
#' 
#' # Phylogenetic network diversity for interaction data with two target sizes (1500 and 3000)
#' data(beetles_plotA)
#' data(beetles_row_tree)
#' output_est_qiPD <- estimateD.link(beetles_plotA, diversity = 'PD', nboot = 10, 
#'                                   base = "size", level = c(1500, 3000), row.tree = beetles_row_tree)
#' output_est_qiPD
#' 
#' ## Functional network diversity for interaction data with two target coverages (93% and 97%)
#' data(beetles_plotA)
#' data(beetles_row_distM)
#' output_est_qiFD = estimateD.link(data = beetles_plotA, diversity = 'FD', q = c(0, 1, 2),
#'                                  base = "coverage", level = c(0.93, 0.97), nboot = 10,
#'                                  row.distM = beetles_row_distM, FDtype = "AUC")
#'output_est_qiFD
#' 
#' 
#' }
#' @export


estimateD.link = function(data, diversity = 'TD', q = c(0, 1, 2), base = "coverage",
                          level = NULL, nboot = 50, conf = 0.95, 
                          row.tree = NULL, col.tree = NULL, PDtype = 'meanPD', 
                          row.distM = NULL, col.distM = NULL, FDtype = "AUC", FDtau = NULL){
  
  datatype = "abundance"
  if(inherits(data, c("data.frame", "matrix", "integer"))) data = list(Network1 = data)
  
  # for(i in 1:length(data)){
  #   if(nrow(data[[i]]) > ncol(data[[i]])){
  #     data[[i]] <- as.data.frame(t(data[[i]]))
  #   }else{
  #     data[[i]] <- data[[i]]
  #   }
  # }
  
  if(diversity == 'TD'){

    div = lapply(1:length(data), function(i){
      x = data[[i]]
      assemblage = names(data)[[i]]
      long = as.matrix(x)%>%c()
      iNEXT.3D::estimate3D(long, q=q,datatype=datatype, base=base,
                           diversity = 'TD', nboot = nboot,conf=conf,level = level)%>%
        mutate(Assemblage = assemblage)
    })%>%do.call("rbind",.)
    
    div <- rename(div, c(qiTD = "qTD", qiTD.LCL = "qTD.LCL", qiTD.UCL = "qTD.UCL"))

    return(div)
  }else if(diversity == 'PD'){
    
    if(datatype=='abundance'){

      # if(inherits(data, c("data.frame", "matrix", "integer"))) data = list(Network1 = data)

      if(inherits(data, "list")){
        if(is.null(names(data))) region_names = paste0("Network_", 1:length(data)) else region_names = names(data)
        Ns = sapply(data, ncol)
        data_list = data
      }

    }
    if(is.null(conf)) conf = 0.95
    tmp = qnorm(1 - (1 - conf)/2)
    
    if (is.null(level) & base == "size") {
      
      level <- sapply(data, function(x) 2 * sum(x)) %>% min
    }else if (is.null(level) & base == "coverage") {
      
      level <- sapply(data, function(x) {
        ni <- sum(x)
        iNEXT.3D:::Coverage(data = x, datatype = datatype, m = 2 * ni)
      })
      
      level <- min(level)
    }
    
    for_each_region = function(data_2d, region_name, N){
      if (datatype=='abundance') {
        
        n = sum(data_2d)
        if(base == 'coverage'){
          size_m = sapply(level, function(i) coverage_to_size(data_2d, i, datatype='abundance'))
        }else{
          if(is.null(level)){
            size_m = n
          }else{
            size_m = level
          }
          level = iNEXT.3D:::Coverage(data_2d,m= size_m, datatype = 'abundance')

        }


        ref= iNEXT.3D:::Coverage(data_2d, m = n, datatype = 'abundance')
        #
        aL_table = create.aili(data_2d, row.tree = row.tree, col.tree = col.tree) %>%
          select(branch.abun, branch.length, tgroup)%>%
          filter(branch.abun>0)
        
        ## boot
        tbar <- sum(aL_table$branch.length*aL_table$branch.abun)/n
        
        qPDm <-iNEXT.3D:::PhD.m.est(ai = aL_table$branch.abun,
                                    Lis = aL_table$branch.length%>%as.matrix(),
                                    m = size_m,
                                    q = q,nt = n, reft = tbar,cal = PDtype) %>% as.vector()

        
        if(nboot >1 ){
          boot.sam <- sample.boot.phy(data_2d,nboot,row.tree = row.tree,col.tree = col.tree)
          PD.sd <- lapply(1:length(boot.sam), function(i){
            if(base == "size"){
              tmp = iNEXT.3D:::PhD.m.est(ai = boot.sam[[i]]$branch.abun,
                                         Lis = boot.sam[[i]]$branch.length%>%as.matrix(),
                                         m = size_m,
                                         q = q,nt = n, reft = tbar, cal = PDtype)%>%
                as.vector()%>%as.data.frame()
            }else{
              
              x_B = boot.sam[[i]] %>% filter(tgroup == "Tip") %>% .$branch.abun %>% as.matrix()
              ai_B <- boot.sam[[i]]$branch.abun %>% as.matrix()
              Li_b <- boot.sam[[i]]$branch.length %>% as.matrix()
              colnames(Li_b) = paste0("T",tbar)
              isn0 <- ai_B > 0
              
              tmp = iNEXT.3D:::invChatPD_abu(x = x_B, 
                                             ai = ai_B[isn0], 
                                             Lis = Li_b[isn0, , drop = F], 
                                             q = q, 
                                             Cs = level, 
                                             n = sum(x_B),
                                             reft = tbar, 
                                             cal = PDtype)$qPD%>%as.vector()%>%as.data.frame()
            }
            
            return(tmp)
          })%>%
            abind(along=3) %>% apply(1:2, sd)%>%as.vector()
        }else{
          PD.sd = rep(NA, length(qPDm))
        }
        
        
        ##
        len = length(q)
        res = data.frame(Assemblage = rep(region_name,len*length(size_m)),
                         Order.q = q,
                         SC = as.vector(sapply(1:length(size_m), function(i) rep(level[i],3))),
                         m = as.vector(sapply(1:length(size_m), function(i) rep(size_m[i],3))),
                         Method = rep(ifelse(level > ref, 'Extrapolation', ifelse(level == ref, 'Observed', 'Rarefaction')), each = len),
                         qPD = qPDm,
                         s.e. = PD.sd,
                         qPD.LCL = qPDm-tmp*PD.sd,
                         qPD.UCL = qPDm+tmp*PD.sd,
                         Reftime = tbar, 
                         Type = PDtype
                         )
        return(res)
      }
    }

    output = lapply(1:length(data), function(i) for_each_region(data_2d = data_list[[i]],
                                                                region_name = region_names[i], N = Ns[i]))%>%
      do.call('rbind',.)

    output <- rename(output, c(qiPD = "qPD", qiPD.LCL = "qPD.LCL", qiPD.UCL = "qPD.UCL"))
    
    return(output)
    
  }else if(diversity == 'FD'& FDtype == 'tau_values'){
    
    output = estimatelinkFD(data, row.distM = row.distM, col.distM = col.distM, datatype = datatype, q = q,
                            base = base, threshold = FDtau, level = level, nboot = nboot,
                            conf = conf)
    
    output <- rename(output, c(qiFD = "qFD", qiFD.LCL = "qFD.LCL", qiFD.UCL = "qFD.UCL"))
    
    return(output)
    
  }else if(diversity == 'FD'& FDtype == 'AUC'){
    
    output = estimatelinkAUC(data, row.distM = row.distM, col.distM = col.distM, datatype = datatype, q = q,
                             base = base, level = level, nboot = nboot,
                             conf = conf)
    
    output <- rename(output, c(qiFD = "qFD", qiFD.LCL = "qFD.LCL", qiFD.UCL = "qFD.UCL"))
    
    return(output)
  }
}




# iNEXTbeta.link ---------------------------
#' Interpolation (rarefaction) and extrapolation of network beta diversity
#'
#' \code{iNEXTbeta.link} Interpolation and extrapolation of beta diversity with order q.
#'
#' @param data data can be input as a \code{list} of \code{data.frame}, each \code{data.frame} represents col.species-by-row.species abundance matrix; see example 1 for an example.
#' @param diversity selection of diversity type: \code{'TD'} = 'Taxonomic diversity', \code{'PD'} = 'Phylogenetic diversity', and \code{'FD'} = 'Functional diversity'.
#' @param level a sequence specifying the particular sample coverages (between 0 and 1). Default is \code{seq(0.5, 1, 0.05)}.
#' @param q a numerical vector specifying the diversity orders. Default is \code{c(0,1,2)}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing
#' sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter 0 to skip the bootstrap procedures. Default is \code{30}.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is \code{0.95}.
#' @param col.tree (required only when \code{diversity = "PD"}), a phylogenetic tree of column assemblage in the pooled network column assemblage.
#' @param row.tree (required only when \code{diversity = "PD"}), a phylogenetic tree of row assemblage in the pooled network row assemblage.
#' @param PDtype (required only when \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"}(effective total branch length) or
#' \code{PDtype = "meanPD"}(effective number of equally divergent lineages).Default is \code{"meanPD"}.
#' @param col.distM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species of column assemblage in the pooled network column assemblage.
#' @param row.distM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species of row assemblage in the pooled network row assemblage.
#' @param FDtype (required only when \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_value"} for FD under specified threshold value, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{"AUC"}.
#' @param FDtau (required only when \code{diversity = "FD"} and \code{FDtype = "tau_value"}), a numerical vector between 0 and 1 specifying tau value (threshold levels). If \code{NULL} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy).
#' @param FDcut_number (required only when \code{diversity = "FD"} and \code{FDtype = "AUC"}), a numeric number to split zero to one into several equal-spaced length. Default is \code{30}.
#' @param by_pair a logical variable specifying whether to perform diversity decomposition for all pairs of assemblages or not. If \code{by_pair = TRUE}, alpha/beta/gamma diversity will be computed for all pairs of assemblages in the input data; if \code{by_pair = FALSE}, alpha/beta/gamma diversity will be computed for multiple assemblages (i.e, more than two assemblages) in the input data. Default is \code{FALSE}. 
#' @return A list of seven matrices with three diversity dimensions and four dissimilarity measures.
#' \item{Dataset}{the name of datasets.}
#' \item{Order.q}{the network diversity order of q.}
#' \item{SC}{the target standardized coverage value. The observed coverage and extrapolation limit for beta diversity are defined the same as those for alpha diversity. For \code{q = 0}, the extrapolation can be extended to a maximum coverage value \code{C(2n, alpha)} = coverage value of twice the alpha reference sample size; for \code{q = 1} and \code{2}, target coverage can be extended to \code{1} (complete coverage) if data are not sparse.}
#' \item{Size}{the corresponding sample size for the standardized coverage value.}
#' \item{Alpha/Beta/Gamma/Dissimilarity}{the estimated diversity or dissimilarity of order q for the target coverage value. The estimate for complete coverage (or \code{level = 1}) represents the estimated asymptotic network diversity.}
#' \item{Method}{Rarefaction, Observed, or Extrapolation, depending on whether the target coverage is less than, equal to, or greater than the coverage of the reference sample. (For beta diversity, observed coverage is defined as the coverage of the alpha reference sample).}
#' \item{s.e.}{standard error of network diversity estimate.}
#' \item{LCL, UCL}{the bootstrap lower and upper confidence limits for the network diversity of order q at the specified level (with a default value of \code{0.95}).}
#' \item{Diversity}{\code{"TD"} (taxonomic network diversity), \code{"PD"} (phylogenetic diversity of effective total branch length), \code{"meanPD"} (phylogenetic diversity of effective number of equally divergent lineages), \code{"FD_tau"} (functional diversity under a single tau), \code{"FD_AUC"} (functional diversity by integrating all threshold values between zero and one.}
#' 
#' 
#' @examples
#' ## Taxonomic network diversity for interaction data
#' # Coverage-based standardized alpha/beta/gamma/dissimilarity network diversity estimates and 
#' # related statistics
#' data(beetles_plotA)
#' output_beta_qiTD = iNEXTbeta.link(data = beetles_plotA, diversity = 'TD', level = NULL, 
#'                                   q = c(0, 1, 2))
#' output_beta_qiTD
#'
#' ## Phylogenetic network diversity for interaction data
#' # Coverage-based standardized alpha/beta/gamma/dissimilarity network diversity estimates and
#' # related statistics
#' data(beetles_plotA)
#' data(beetles_row_tree)
#' output_beta_qiPD = iNEXTbeta.link(data = beetles_plotA, diversity = 'PD', level = NULL,
#'                                   q = c(0, 1, 2),row.tree = beetles_row_tree, nboot = 10)
#' output_beta_qiPD
#'
#'
#' ## Functional network diversity for interaction data
#' # Coverage-based standardized alpha/beta/gamma/dissimilarity network diversity estimates and 
#' # related statistics
#' data(beetles_plotA)
#' data(beetles_row_distM)
#' output_beta_qiFD = iNEXTbeta.link(data = beetles_plotA, diversity = 'FD', level = NULL, 
#'                                   q = c(0, 1, 2), nboot = 10, 
#'                                   row.distM = beetles_row_distM, FDtype = "AUC")
#' output_beta_qiFD
#'
#' @references
#' Chao, A., Thorn, S., Chiu, C.-H., Moyes, F., Hu, K.-H., Chazdon, R. L., Wu, J., Dornelas, M., Zeleny, D., Colwell, R. K., and Magurran, A. E. (2023). Rarefaction and extrapolation with beta diversity under a framework of Hill numbers: the iNEXT.beta3D standardization. 
#' \emph{Ecological Monographs e1588.} \cr\cr
#' @export

iNEXTbeta.link = function(data, diversity = 'TD', level = NULL,
                          q = c(0, 1, 2), nboot = 30, conf = 0.95,
                          row.tree = NULL, col.tree = NULL, PDtype = 'meanPD', 
                          row.distM = NULL, col.distM = NULL, FDtype = "AUC", FDtau = NULL, FDcut_number = 30, by_pair = FALSE){
  
  datatype = 'abundance'
  # for(i in 1:length(data)){
  #   if(nrow(data[[i]]) > ncol(data[[i]])){
  #     data[[i]] <- as.data.frame(t(data[[i]]))
  #   }else{
  #     data[[i]] <- data[[i]]
  #   }
  # }
  
  if (inherits(data[[1]], c("data.frame", "matrix"))) {dat = list("Dataset" = data); } else {dat = data}

  combined_list = lapply(dat, function(y){

    long = ready4beta(y)%>%filter_all(any_vars(. != 0))
    rownames(long) = rownames(long)%>%gsub('\\.','_',.)
    colnames(long) = names(y)
    return(long)
  })

  if(by_pair == FALSE){
    if(diversity == 'TD'){
    
    dissimilarity <- iNEXTbeta3D(data = combined_list, diversity = 'TD',level = level, datatype = datatype,
                                 q = q ,nboot = nboot, conf = conf)
    
  }else if(diversity == 'PD'){
    
    if(!is.null(row.tree)){row.tree$tip.label = gsub('\\.', '_',row.tree$tip.label)}
    if(!is.null(col.tree)){col.tree$tip.label = gsub('\\.', '_',col.tree$tip.label)}
    
    dissimilarity = iNEXTbeta.PDlink(data = combined_list, level = level, datatype = datatype,
                                     q = q ,row.tree = row.tree,col.tree = col.tree, nboot = nboot, PDtype = PDtype)
    class(dissimilarity) = 'iNEXTbeta3D'
    
  }else if(diversity == 'FD' & FDtype == 'tau_value'){
    
    row_sp = c()
    col_sp = c()
    for(i in 1:length(dat)){
      for(j in 1:length(dat[[i]])){
        row_sp = c(row_sp,rownames(dat[[i]][[j]]))
        col_sp = c(col_sp,colnames(dat[[i]][[j]]))
      }
      
      
    }
    row_num = length(unique(row_sp))
    col_num = length(unique(col_sp))
    
    if(is.null(row.distM)){
      rdd = matrix(1,ncol = row_num,nrow = row_num)
      diag(rdd) = 0
      rownames(rdd) = unique(row_sp)
      colnames(rdd) = unique(row_sp)
      row.distM =  rdd}
    if(is.null(col.distM)){
      cdd = matrix(1,ncol = col_num,nrow = col_num)
      diag(cdd) = 0
      rownames(cdd) = unique(col_sp)
      colnames(cdd) = unique(col_sp)
      col.distM =  cdd}
    row.distM = as.matrix(row.distM)
    col.distM = as.matrix(col.distM)
    
    distM =  1-(1-row.distM[rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM))),rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM)))])*(1-col.distM[rep(1:nrow(col.distM),nrow(row.distM)),rep(1:nrow(col.distM),nrow(row.distM))])
    distM_name = paste0(rep(rownames(row.distM),rep(ncol(col.distM),nrow(row.distM))),"*",rep(colnames(col.distM),3))
    colnames(distM) = distM_name
    rownames(distM) = distM_name
    
    dissimilarity <- iNEXTbeta3D(data = combined_list, diversity = 'FD',level = level, datatype = datatype,
                                 q = q ,nboot = nboot, conf = conf, FDdistM = distM, FDtype = "tau_value", FDtau = FDtau)
    
  }else if(diversity == 'FD' & FDtype == 'AUC'){
    
    row_sp = c()
    col_sp = c()
    for(i in 1:length(dat)){
      for(j in 1:length(dat[[i]])){
        row_sp = c(row_sp,rownames(dat[[i]][[j]]))
        col_sp = c(col_sp,colnames(dat[[i]][[j]]))
      }
      
      
    }
    row_num = length(unique(row_sp))
    col_num = length(unique(col_sp))
    
    if(is.null(row.distM)){
      rdd = matrix(1,ncol = row_num,nrow = row_num)
      diag(rdd) = 0
      rownames(rdd) = unique(row_sp)
      colnames(rdd) = unique(row_sp)
      row.distM =  rdd}
    if(is.null(col.distM)){
      cdd = matrix(1,ncol = col_num,nrow = col_num)
      diag(cdd) = 0
      rownames(cdd) = unique(col_sp)
      colnames(cdd) = unique(col_sp)
      col.distM =  cdd}
    
    
    row.distM = as.matrix(row.distM)
    col.distM = as.matrix(col.distM)
    
    distM = 1-(1-row.distM[rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM))),rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM)))])*(1-col.distM[rep(1:nrow(col.distM),nrow(row.distM)),rep(1:nrow(col.distM),nrow(row.distM))])
    distM_name = paste0(rep(rownames(row.distM),rep(ncol(col.distM),nrow(row.distM))),"*",rep(colnames(col.distM),3))
    colnames(distM) = distM_name
    rownames(distM) = distM_name
    
    
    dissimilarity <- iNEXTbeta3D(data = combined_list, diversity = 'FD', level = level, datatype = datatype,
                                 q = q, nboot = nboot, conf = conf, FDdistM = distM, FDcut_number = FDcut_number)
  }
    dissimilarity[[1]]$gamma$Dataset <- paste(names(combined_list[[1]]),collapse = " vs. ")
    dissimilarity[[1]]$alpha$Dataset <- paste(names(combined_list[[1]]),collapse = " vs. ")
    dissimilarity[[1]]$beta$Dataset  <- paste(names(combined_list[[1]]),collapse = " vs. ")
    dissimilarity[[1]]$`1-C`$Dataset <- paste(names(combined_list[[1]]),collapse = " vs. ")
    dissimilarity[[1]]$`1-U`$Dataset <- paste(names(combined_list[[1]]),collapse = " vs. ")
    dissimilarity[[1]]$`1-V`$Dataset <- paste(names(combined_list[[1]]),collapse = " vs. ")
    dissimilarity[[1]]$`1-S`$Dataset <- paste(names(combined_list[[1]]),collapse = " vs. ")
    
  }else if(by_pair == TRUE){
    
    t <- length(names(combined_list[[1]]))
    dis <- list()
    if(diversity == 'TD'){
      for(i in 1:(t-1)){
        for(j in (i+1):t){
          com_list <- combined_list[[1]][,c(i,j)]
          dis_test <- iNEXTbeta3D(data = com_list, diversity = 'TD',level = level, datatype = datatype,
                                  q = q ,nboot = nboot, conf = conf)
          names(dis_test) <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$gamma$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$alpha$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$beta$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$`1-C`$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$`1-U`$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$`1-V`$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$`1-S`$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis <- c(dis, dis_test)
          
        }
      }
      dissimilarity <- list(Dataset = list())
      for(i in 1:7){
        for(j in 1:length(dis)){
          if(j == 1){
            dissimilarity[[1]][[i]] <-  dis[[j]][[i]]
          }else{
            dissimilarity[[1]][[i]] <- rbind(dissimilarity[[1]][[i]], dis[[j]][[i]]) 
          }
        }
        names(dissimilarity[[1]])[i] <- names(dis[[1]])[i]
      }
      
      
    }else if(diversity == 'PD'){
      
      if(!is.null(row.tree)){row.tree$tip.label = gsub('\\.', '_',row.tree$tip.label)}
      if(!is.null(col.tree)){col.tree$tip.label = gsub('\\.', '_',col.tree$tip.label)}
      
      for(i in 1:(t-1)){
        for(j in (i+1):t){
          com_list <- combined_list[[1]][,c(i,j)]
          dis_test <- iNEXTbeta.PDlink(data = com_list, level = level, datatype = datatype,
                                       q = q ,row.tree = row.tree,col.tree = col.tree, nboot = nboot, PDtype = PDtype)
          names(dis_test) <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$gamma$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$alpha$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$beta$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$`1-C`$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$`1-U`$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$`1-V`$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$`1-S`$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis <- c(dis, dis_test)
          
        }
      }
      dissimilarity <- list(Dataset = list())
      for(i in 1:7){
        for(j in 1:length(dis)){
          if(j == 1){
            dissimilarity[[1]][[i]] <-  dis[[j]][[i]]
          }else{
            dissimilarity[[1]][[i]] <- rbind(dissimilarity[[1]][[i]], dis[[j]][[i]]) 
          }
        }
        names(dissimilarity[[1]])[i] <- names(dis[[1]])[i]
      }
      
      class(dissimilarity) = 'iNEXTbeta3D'
      
    }else if(diversity == 'FD' & FDtype == 'tau_value'){
      row_sp = c()
      col_sp = c()
      for(i in 1:length(dat)){
        for(j in 1:length(dat[[i]])){
          row_sp = c(row_sp,rownames(dat[[i]][[j]]))
          col_sp = c(col_sp,colnames(dat[[i]][[j]]))
        }
        
        
      }
      row_num = length(unique(row_sp))
      col_num = length(unique(col_sp))
      
      if(is.null(row.distM)){
        rdd = matrix(1,ncol = row_num,nrow = row_num)
        diag(rdd) = 0
        rownames(rdd) = unique(row_sp)
        colnames(rdd) = unique(row_sp)
        row.distM =  rdd}
      if(is.null(col.distM)){
        cdd = matrix(1,ncol = col_num,nrow = col_num)
        diag(cdd) = 0
        rownames(cdd) = unique(col_sp)
        colnames(cdd) = unique(col_sp)
        col.distM =  cdd}
      row.distM = as.matrix(row.distM)
      col.distM = as.matrix(col.distM)
      
      distM =  1-(1-row.distM[rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM))),rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM)))])*(1-col.distM[rep(1:nrow(col.distM),nrow(row.distM)),rep(1:nrow(col.distM),nrow(row.distM))])
      distM_name = paste0(rep(rownames(row.distM),rep(ncol(col.distM),nrow(row.distM))),"*",rep(colnames(col.distM),3))
      colnames(distM) = distM_name
      rownames(distM) = distM_name
      
      for(i in 1:(t-1)){
        for(j in (i+1):t){
          com_list <- combined_list[[1]][,c(i,j)]
          dis_test <- iNEXTbeta3D(data = com_list, diversity = 'FD',level = level, datatype = datatype,
                                  q = q ,nboot = nboot, conf = conf, FDdistM = distM, FDtype = "tau_value", FDtau = FDtau)
          names(dis_test) <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$gamma$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$alpha$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$beta$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$`1-C`$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$`1-U`$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$`1-V`$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$`1-S`$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis <- c(dis, dis_test)
          
        }
      }
      dissimilarity <- list(Dataset = list())
      for(i in 1:7){
        for(j in 1:length(dis)){
          if(j == 1){
            dissimilarity[[1]][[i]] <-  dis[[j]][[i]]
          }else{
            dissimilarity[[1]][[i]] <- rbind(dissimilarity[[1]][[i]], dis[[j]][[i]]) 
          }
        }
        names(dissimilarity[[1]])[i] <- names(dis[[1]])[i]
      }
      
      
    }else if(diversity == 'FD' & FDtype == 'AUC'){
      
      row_sp = c()
      col_sp = c()
      for(i in 1:length(dat)){
        for(j in 1:length(dat[[i]])){
          row_sp = c(row_sp,rownames(dat[[i]][[j]]))
          col_sp = c(col_sp,colnames(dat[[i]][[j]]))
        }
        
        
      }
      row_num = length(unique(row_sp))
      col_num = length(unique(col_sp))
      
      if(is.null(row.distM)){
        rdd = matrix(1,ncol = row_num,nrow = row_num)
        diag(rdd) = 0
        rownames(rdd) = unique(row_sp)
        colnames(rdd) = unique(row_sp)
        row.distM =  rdd}
      if(is.null(col.distM)){
        cdd = matrix(1,ncol = col_num,nrow = col_num)
        diag(cdd) = 0
        rownames(cdd) = unique(col_sp)
        colnames(cdd) = unique(col_sp)
        col.distM =  cdd}
      
      
      row.distM = as.matrix(row.distM)
      col.distM = as.matrix(col.distM)
      
      distM = 1-(1-row.distM[rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM))),rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM)))])*(1-col.distM[rep(1:nrow(col.distM),nrow(row.distM)),rep(1:nrow(col.distM),nrow(row.distM))])
      distM_name = paste0(rep(rownames(row.distM),rep(ncol(col.distM),nrow(row.distM))),"*",rep(colnames(col.distM),3))
      colnames(distM) = distM_name
      rownames(distM) = distM_name
      
      for(i in 1:(t-1)){
        for(j in (i+1):t){
          com_list <- combined_list[[1]][,c(i,j)]
          dis_test <- iNEXTbeta3D(data = com_list, diversity = 'FD', level = level, datatype = datatype,
                                  q = q ,nboot = nboot, conf = conf, FDdistM = distM, FDcut_number = FDcut_number)
          names(dis_test) <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$gamma$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$alpha$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$beta$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$`1-C`$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$`1-U`$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$`1-V`$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis_test[[1]]$`1-S`$Dataset <- paste(names(com_list),collapse = " vs. ")
          dis <- c(dis, dis_test)
          
        }
      }
      dissimilarity <- list(Dataset = list())
      for(i in 1:7){
        for(j in 1:length(dis)){
          if(j == 1){
            dissimilarity[[1]][[i]] <-  dis[[j]][[i]]
          }else{
            dissimilarity[[1]][[i]] <- rbind(dissimilarity[[1]][[i]], dis[[j]][[i]]) 
          }
        }
        names(dissimilarity[[1]])[i] <- names(dis[[1]])[i]
      }
      
    }
    }
  
  return(dissimilarity)
}


# ggiNEXTbeta.link -------------------------------------------------------------------
#' ggplot2 extension for output from \code{iNEXTbeta.link}
#'
#' \code{ggiNEXTbeta.link}: ggplot for Interpolation and extrapolation of beta diversity with order q
#'
#' @param output the output from \code{"iNEXTbeta.link"}
#' @param type selection of plot type : \code{type = 'B'} for plotting the gamma, alpha, and beta diversity ;
#' \code{type = 'D'} for plotting 4 turnover dissimilarities.
# @param scale Are scales shared across all facets (\code{"fixed"}), or do they vary across rows (\code{"free_x"}), columns (\code{"free_y"}), or both rows and columns (\code{"free"})? Default is \code{"free"}.
#'
#' @return a figure for gamma, alpha, and beta diversity or four dissimilarity measures.
#' 
#' @examples
#' ## Taxonomic network diversity for interaction data
#' # Plot coverage-based standardized alpha/beta/gamma/dissimilarity network diversity estimates and 
#' # related statistics
#' data(beetles_plotA)
#' output_beta_qiTD = iNEXTbeta.link(data = beetles_plotA, diversity = 'TD', level = NULL, 
#'                                   q = c(0, 1, 2))
#' ggiNEXTbeta.link(output_beta_qiTD, type = 'B')
#' ggiNEXTbeta.link(output_beta_qiTD, type = 'D')
#'
#' ## Phylogenetic network diversity for interaction data
#' # Plot coverage-based standardized alpha/beta/gamma/dissimilarity network diversity estimates and 
#' # related statistics
#' data(beetles_plotA)
#' data(beetles_row_tree)
#' output_beta_qiPD = iNEXTbeta.link(data = beetles_plotA, diversity = 'PD', level = NULL, 
#'                                   q = c(0, 1, 2), row.tree = beetles_row_tree, nboot = 10)
#' ggiNEXTbeta.link(output_beta_qiPD, type = 'B')
#' ggiNEXTbeta.link(output_beta_qiPD, type = 'D')
#'
#'
#' ## Functional network diversity for interaction data
#' # Plot coverage-based standardized alpha/beta/gamma/dissimilarity network diversity estimates and 
#' # related statistics
#' data(beetles_plotA)
#' data(beetles_row_distM)
#' output_beta_qiFD = iNEXTbeta.link(data = beetles_plotA, diversity = 'FD', level = NULL, 
#'                                   q = c(0, 1, 2), row.distM = beetles_row_distM, FDtype = "AUC", 
#'                                   nboot = 10)
#' ggiNEXTbeta.link(output_beta_qiFD, type = 'B')
#' ggiNEXTbeta.link(output_beta_qiFD, type = 'D')
#'
#'
#' @export

ggiNEXTbeta.link <- function(output, type = c('B', 'D')){

  if (type == 'B'){

    gamma = lapply(output, function(y) y[["gamma"]]) %>% do.call(rbind,.) %>% rename("Estimate" = "Gamma") %>% mutate(div_type = "Gamma") %>% as_tibble()
    alpha = lapply(output, function(y) y[["alpha"]]) %>% do.call(rbind,.) %>% rename("Estimate" = "Alpha") %>% mutate(div_type = "Alpha") %>% as_tibble()
    beta =  lapply(output, function(y) y[["beta"]])  %>% do.call(rbind,.) %>% rename("Estimate" = "Beta")  %>% mutate(div_type = "Beta")  %>% as_tibble()
    # beta = beta %>% filter(Method != 'Observed')
    
    alpha[alpha == 'Observed_SC(n, alpha)'] = 'Observed'
    alpha[alpha == 'Extrap_SC(2n, alpha)'] = 'Extrapolation'
    gamma[gamma == 'Observed_SC(n, gamma)'] = 'Observed'
    gamma[gamma == 'Extrap_SC(2n, gamma)'] = 'Extrapolation'
    beta[beta == 'Observed_SC(n, alpha)'] = 'Observed'
    beta[beta == 'Extrap_SC(2n, alpha)'] = 'Extrapolation'
    
    df = rbind(gamma, alpha, beta)
    for (i in unique(gamma$Order.q)) df$Order.q[df$Order.q == i] = paste0('q = ', i)
    df$div_type <- factor(df$div_type, levels = c("Gamma","Alpha","Beta"))
    
    id_obs = which(df$Method == 'Observed')
    
    if (length(id_obs) > 0) {
      for (i in 1:length(id_obs)) {
        
        new = df[id_obs[i],]
        new$SC = new$SC - 0.0001
        new$Method = 'Rarefaction'
        
        newe = df[id_obs[i],]
        newe$SC = newe$SC + 0.0001
        newe$Method = 'Extrapolation'
        
        df = rbind(df, new, newe)
        
      }
    }
    
    if (unique(output[[1]]$gamma$Diversity) == 'TD') { ylab = "Taxonomic diversity" }
    if (unique(output[[1]]$gamma$Diversity) %in% c('PD','meanPD')) { ylab = "Phylogenetic diversity" }
    if (unique(output[[1]]$gamma$Diversity) == 'FD_tau') { ylab = "Functional diversity" }
    if (unique(output[[1]]$gamma$Diversity) == 'FD_AUC') { ylab = "Functional diversity (AUC)" }

  }
  if (type == 'D'){

    C = lapply(output, function(y) y[["1-C"]]) %>% do.call(rbind,.) %>% rename("Estimate" = "Dissimilarity") %>% mutate(div_type = "1-CqN") %>% as_tibble()
    U = lapply(output, function(y) y[["1-U"]]) %>% do.call(rbind,.) %>% rename("Estimate" = "Dissimilarity") %>% mutate(div_type = "1-UqN") %>% as_tibble()
    V = lapply(output, function(y) y[["1-V"]]) %>% do.call(rbind,.) %>% rename("Estimate" = "Dissimilarity") %>% mutate(div_type = "1-VqN") %>% as_tibble()
    S = lapply(output, function(y) y[["1-S"]]) %>% do.call(rbind,.) %>% rename("Estimate" = "Dissimilarity") %>% mutate(div_type = "1-SqN") %>% as_tibble()
    # C = C %>% filter(Method != 'Observed')
    # U = U %>% filter(Method != 'Observed')
    # V = V %>% filter(Method != 'Observed')
    # S = S %>% filter(Method != 'Observed')
    C[C == 'Observed_SC(n, alpha)'] = U[U == 'Observed_SC(n, alpha)'] = V[V == 'Observed_SC(n, alpha)'] = S[S == 'Observed_SC(n, alpha)'] = 'Observed'
    C[C == 'Extrap_SC(2n, alpha)'] = U[U == 'Extrap_SC(2n, alpha)'] = V[V == 'Extrap_SC(2n, alpha)'] = S[S == 'Extrap_SC(2n, alpha)'] = 'Extrapolation'
    
    df = rbind(C, U, V, S)
    for (i in unique(C$Order.q)) df$Order.q[df$Order.q == i] = paste0('q = ', i)
    df$div_type <- factor(df$div_type, levels = c("1-CqN", "1-UqN", "1-VqN", "1-SqN"))
    
    id_obs = which(df$Method == 'Observed')
    
    if (length(id_obs) > 0) {
      for (i in 1:length(id_obs)) {
        
        new = df[id_obs[i],]
        new$SC = new$SC - 0.0001
        new$Method = 'Rarefaction'
        
        newe = df[id_obs[i],]
        newe$SC = newe$SC + 0.0001
        newe$Method = 'Extrapolation'
        
        df = rbind(df, new, newe)
        
      }
    }
    
    if (unique(output[[1]]$gamma$Diversity) == 'TD') { ylab = "Taxonomic dissimilarity" }
    if (unique(output[[1]]$gamma$Diversity) %in% c('PD','meanPD')) { ylab = "Phylogenetic dissimilarity" }
    if (unique(output[[1]]$gamma$Diversity) == 'FD_tau') { ylab = "Functional dissimilarity" }
    if (unique(output[[1]]$gamma$Diversity) == 'FD_AUC') { ylab = "Functional dissimilarity (AUC)" }

  }
  
  lty = c(Rarefaction = "solid", Extrapolation = "dashed")
  # lty = c(Rarefaction = "solid", Extrapolation = "twodash")

  df$Method = factor(df$Method, levels = c('Rarefaction', 'Extrapolation', 'Observed'))

  double_size = unique(df[df$Method=="Observed",]$Size)*2
  double_extrapolation = df %>% filter(Method=="Extrapolation" & round(Size) %in% double_size)
  
  # Check if the number of unique 'Assemblage' is 8 or less
  if (length(unique(df$Dataset)) <= 8){
    cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  }else{
    # If there are more than 8 assemblages, start with the same predefined color palette
    # Then extend the palette by generating additional colors using the 'ggplotColors' function
    cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#330066", "#CC79A7", "#0072B2", "#D55E00"))
    cbPalette <- c(cbPalette, ggplotColors(length(unique(df$Dataset))-8))
  }
  point_size = 2

  ggplot(data = df, aes(x = SC, y = Estimate, col = Dataset)) +
    geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Dataset, col = NULL), alpha=0.4) +
    geom_line(data = subset(df, Method!='Observed'), aes(linetype=Method), size=1.1) +
    scale_linetype_manual(values = lty) +
    scale_fill_manual(values = cbPalette) + scale_colour_manual(values = cbPalette)+
    # geom_line(lty=2) +
    geom_point(data = subset(df, Method=='Observed' & div_type=="Gamma"),shape=19, size=point_size) +
    geom_point(data = subset(df, Method=='Observed' & div_type!="Gamma"),shape=1, size=point_size,stroke=1.5)+
    geom_point(data = subset(double_extrapolation, div_type == "Gamma"),shape=17, size=point_size) +
    geom_point(data = subset(double_extrapolation, div_type!="Gamma"),shape=2, size=point_size,stroke=1.5) +
    facet_grid(div_type~Order.q, scales = 'free') +
    # facet_wrap(div_type~Order.q, scales = scale, switch="both") +
    theme_bw() + 
    theme(legend.position = "bottom", 
          legend.title = element_blank(),
          strip.text = element_text(size = 15, face = 'bold'),
          axis.title = element_text(hjust = 0.5, size = 15, face = 'bold'),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.box = "vertical",
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-10, -10, -5, -10),
          legend.text = element_text(size = 13),
          plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
    labs(x='Sample coverage', y=ylab) +
    guides(linetype = guide_legend(keywidth = 2.5))
}


# Spec.link.ObsAsy  -------------------------------------------------------------------
#' Asymptotic estimation (or observed) of specialization with order q
#' @param data a \code{list} of \code{data.frames}, each \code{data.frames} represents col.species-by-row.species abundance matrix.
#' @param q a numerical vector specifying the diversity orders. Default is \code{seq(0, 2, 0.2)}.
#' @param method a binary calculation method with \code{"Asymptotic"} or \code{"Observed"}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing
#' sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter 0 to skip the bootstrap procedures. Default is \code{30}.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is \code{0.95}.
#' @param E.class an integer vector between 1 to 5.
#' @param level the type of specialization measure: Choose "weighted" for weighted species-level specialization (default), or "network" for network-level specialization.
#' 
#' @return A list of several tables containing estimated (or observed) evenness with order q.\cr
#'         Each tables represents a class of specialization.
#'         \item{Order.q}{the network diversity order of q.}
#'         \item{Specialization}{the specialization of order q.}
#'         \item{s.e.}{standard error of evenness.}
#'         \item{Spec.LCL, Spec.UCL}{the bootstrap lower and upper confidence limits for the evenness of order q at the specified level (with a default value of \code{0.95}).}
#'         \item{Method}{\code{"Asymptotic"} or \code{"Observed"}.}
#'         \item{Dataset}{the dataset name.}
#'         \item{Measure}{specialization class.}
#'         \item{level}{type of specialization measure.}
#'         
#'
#' @examples
#' data(beetles)
#' output_spec_ObsAsy = Spec.link.ObsAsy(beetles)
#' output_spec_ObsAsy
#' @export
#' 

Spec.link.ObsAsy <- function(data, q = seq(0, 2, 0.2),
                          method = "Asymptotic",
                          nboot = 30,
                          conf = 0.95,
                          E.class = c(1:5), level = "weighted"){
  
  datatype = "abundance"
  diversity = 'TD'
  if(inherits(data, c("data.frame", "matrix", "integer"))) data = list(Network1 = data)
  
  long = lapply(data, function(da){da%>%as.data.frame()%>%gather(key = "col_sp", value = "abundance")%>%.[,2]})
  if ( !(method %in% c("Asymptotic", "Observed")) )
    stop("Please select one of below method: 'Asymptotic', 'Observed'", call. = FALSE)
  if ( !(level %in% c("weighted", "network")) )
    stop("Please select one of below level: 'weighted', 'network'", call. = FALSE)
  
  if (level == "weighted"){
    index = lapply(data,function(x){
      sub = apply(x,2,function(i){sum(i!=0)})
      return(index = sum(sub == 1))
    })
    if(sum(unlist(index)) != 0){
      stop("Under the current setting level = 'weighted', cases where a column species interacts with only a single row species cannot be properly computed. Please inspect and preprocess the input list accordingly, or consider using level = 'network' instead.", call. = FALSE)
    }
  }
  
  if(level == "weighted"){
    
    Spec <- lapply(E.class, function(e){
      each_class = lapply(1:length(data), function(i){
        
        assemblage = names(data)[[i]]
        if(is.null(assemblage)){
          assemblage = paste0("Network",i)
        }
        sub = data[[i]]
        if(method == "Observed"){
          res = iNEXT.4steps::Evenness(sub, q = q,datatype = datatype,
                                       method = method, nboot=0, E.class = e)
        }else if(method == "Asymptotic"){
          res = Evenness_asym(sub, q = q,datatype = datatype, nboot=0, E.class = e)
        }
        
        wk = colSums(sub)/sum(sub)
        even = sapply(q, function(r){
          tmp = res[[1]] |> filter(Order.q == r)
          return(sum(tmp$Evenness*wk))
        })
        
        if(nboot > 1){
          

          # plan(multisession)
          se = do.call(rbind,future_lapply(1:nboot, function(i){
            
            bootstrap_population = iNEXT.beta3D:::bootstrap_population_multiple_assemblage(data = sub, data_gamma = rowSums(sub), datatype = 'abundance')
            bootstrap_sample = sapply(1:ncol(sub), function(k) rmultinom(n = 1, size = sum(sub[,k]), prob = bootstrap_population[,k]))
            
            if(method == "Observed"){
              res_boost = iNEXT.4steps::Evenness(bootstrap_sample, q = q,datatype = datatype,
                                                 method = method, nboot=0, E.class = e)
            }else if(method == "Asymptotic"){
              res_boost = Evenness_asym(sub, q = q,datatype = datatype, nboot=0, E.class = e)
            }
            

            wk = colSums(bootstrap_sample)/sum(bootstrap_sample)
            even_boost = sapply(q, function(r){
              tmp = res_boost[[1]] |> filter(Order.q == r)
              return(sum(tmp$Evenness*wk))
            }) 
            
            even_boost
            return(even_boost)
          }, future.seed = TRUE)) %>% apply(2,sd) %>% as.data.frame()
          plan(sequential)
          
          colnames(se) = c("s.e.")
          
        }else{
          
          se = matrix(NA, ncol = 1, nrow = length(even))  %>%  as.data.frame()
          colnames(se) = c("s.e.")
          
        }
        
        tmp = qnorm(1 - (1 - conf)/2)
        
        res = cbind('Specialization' = 1-even,se)  %>% 
          mutate(Network = assemblage, Method = method ,Order.q = q,
                 "Spec.LCL" = Specialization - tmp*`s.e.`, "Spec.UCL"= Specialization + tmp*`s.e.`, class = paste0("1 - E",e), level = "Weighted species-level") %>%
          select(c("Order.q", 'Specialization',"s.e.", "Spec.UCL", "Spec.LCL","Method","Network","class","level"))
        
        return(res)
      })%>%do.call("rbind",.)
    })
    
    names(Spec) = paste0("1 - E",E.class)
    
    for(i in 1:length(E.class)){
      Spec[[i]][,4:5] <- Spec[[i]][,c(5,4)]
      names(Spec[[i]])[7:8] <- c("Dataset", "Measure")
    }
    
  }else if(level == "network"){

    
    long = lapply(data, function(da){da%>%as.data.frame()%>%gather(key = "col_sp", value = "abundance")%>%.[,2]})
    
    Spec <- lapply(E.class, function(e){
      each_class = lapply(seq_along(long), function(i){
        if(method == "Observed"){
          res = iNEXT.4steps::Evenness(long[[i]], q = q,datatype = datatype,
                                       method = method, nboot=nboot, E.class = e)
        }else if(method == "Asymptotic"){
          res = Evenness_asym(long[[i]], q = q,datatype = datatype, nboot=nboot, E.class = e)
        }
        res['Coverage'] = NULL
        res = lapply(res, function(each_class){
          each_class%>%
            mutate(Evenness = 1-Evenness, Even.LCL = 1-Even.LCL, Even.UCL = 1-Even.UCL) %>% 
            select(-Assemblage)%>%
            rename('Specialization'='Evenness', 'Spec.UCL' ='Even.LCL', 'Spec.LCL' ='Even.UCL')%>%
            mutate(Network = names(long)[[i]],class = paste0("1 - E",e), level = "Network-level")
        })
        # if(method == "Observed") index = 1
        # if(method == "Estimated") index = 2
        # return(res[[index]]%>%mutate(Assemblage = names(long)[[i]]))
        return(res[[1]])
      }) %>% do.call("rbind",.)
      
      # each_class %>% mutate(class = paste0("1 - E",e))
    })
    names(Spec) = paste0("1 - E",E.class)
    
    for(i in 1:length(E.class)){
      Spec[[i]][,4:5] <- Spec[[i]][,c(5,4)]
      names(Spec[[i]])[4:5] <- c("Spec.LCL", "Spec.UCL")
      names(Spec[[i]])[7:8] <- c("Dataset", "Measure")
    }
  }
  
  return(Spec)
  
}

# Spec.link.est -------------------------------------------------------------------
#' Standardized estimation of specialization with order q
#' @param data a \code{list} of \code{data.frames}, each \code{data.frames} represents col.species-by-row.species abundance matrix.
#' @param q a numerical vector specifying the diversity orders. Default is \code{seq(0, 2, 0.2)}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing
#' sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter 0 to skip the bootstrap procedures. Default is \code{30}.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is \code{0.95}.
#' @param E.class an integer vector between 1 to 5.
#' @param SC a standardized coverage for calculating specialization index. If \code{NULL}, then this function computes the diversity estimates for the minimum sample coverage among all samples extrapolated to double reference sizes (\code{SC = Cmax}).
#' @param level the type of specialization measure: Choose "weighted" for weighted species-level specialization (default), or "network" for network-level specialization.
#'
#' 
#' @return A list of several tables containing    estimated (or observed) evenness with order q.\cr
#'         Each tables represents a class of specialization.
#'         \item{Order.q}{the network diversity order of q.}
#'         \item{Specialization}{the specialization of order q.}
#'         \item{s.e.}{standard error of evenness.}
#'         \item{Spec.LCL, Spec.UCL}{the bootstrap lower and upper confidence limits for the evenness of order q at the specified level (with a default value of \code{0.95}).}
#'         \item{SC}{the target standardized coverage value.}
#'         \item{Dataset}{the dataset name.}
#'         \item{Measure}{specialization class.}
#'         \item{level}{type of specialization measure.}
#'         
#'
#' @examples
#' data(beetles)
#' output_spec = Spec.link.est(beetles)
#' output_spec
#' @export
#' 


Spec.link.est <- function(data, q = seq(0, 2, 0.2),
                          nboot = 30,                      
                          conf = 0.95,
                          E.class = c(1:5),
                          SC = NULL, 
                          level = "weighted"){
  
  method = "Estimated"
  datatype = "abundance"
  diversity = 'TD'
  if(inherits(data, c("data.frame", "matrix", "integer"))) data = list(Network1 = data)
  
  if ( !(level %in% c("weighted", "network")) )
    stop("Please select one of below level: 'weighted', 'network'", call. = FALSE)
  
  if (level == "weighted"){
    index = lapply(data,function(x){
      sub = apply(x,2,function(i){sum(i!=0)})
      return(index = sum(sub == 1))
    })
    if(sum(unlist(index)) != 0){
      stop("Under the current setting level = 'weighted', cases where a column species interacts with only a single row species cannot be properly computed. Please inspect and preprocess the input list accordingly, or consider using level = 'network' instead.", call. = FALSE)
    }
  }
  
  long = lapply(data, function(da){da%>%as.data.frame()%>%gather(key = "col_sp", value = "abundance")%>%.[,2]})
  
  if(method == "Estimated"){
    if (is.null(SC)) SC = sapply(long, function(x) iNEXT.3D:::Coverage(x, datatype = 'abundance', 2 * sum(x))) %>% min
  }else{
    if(is.null(SC)){
      SC = DataInfo.link(data, diversity="TD")$Coverage
    }else{
      SC = SC 
    }
  }
  
  if(level == "weighted"){
    
    Spec <- lapply(E.class, function(e){
      each_class = lapply(1:length(data), function(i){
        
        assemblage = names(data)[[i]]
        if(is.null(assemblage)){
          assemblage = paste0("Network",i)
        }
        sub = data[[i]]
        
        res = iNEXT.4steps::Evenness(sub, q = q,datatype = datatype,
                                     method = method, nboot=0, E.class = e, SC = SC)

        wk = colSums(sub)/sum(sub)
        even = sapply(q, function(r){
          tmp = res[[1]] |> filter(Order.q == r)
          return(sum(tmp$Evenness*wk))
        })
        
        if(nboot > 1){
          
          # plan(multisession)
          se = do.call(rbind,future_lapply(1:nboot, function(i){
            
            bootstrap_population = iNEXT.beta3D:::bootstrap_population_multiple_assemblage(data = sub, data_gamma = rowSums(sub), datatype = 'abundance')
            bootstrap_sample = sapply(1:ncol(sub), function(k) rmultinom(n = 1, size = sum(sub[,k]), prob = bootstrap_population[,k]))
            
            res_boost = iNEXT.4steps::Evenness(bootstrap_sample, q = q,datatype = datatype,
                                         method = method, nboot=0, E.class = e, SC = SC)
            wk = colSums(bootstrap_sample)/sum(bootstrap_sample)
            even_boost = sapply(q, function(r){
              tmp = res_boost[[1]] |> filter(Order.q == r)
              return(sum(tmp$Evenness*wk))
            }) 
            
            even_boost
            return(even_boost)
          }, future.seed = TRUE)) %>% apply(2,sd) %>% as.data.frame()
          plan(sequential)
          
          colnames(se) = c("s.e.")
          
        }else{
          
          se = matrix(NA, ncol = 1, nrow = length(even)) %>% as.data.frame()
          colnames(se) = c("s.e.")
          
        }
        
        tmp = qnorm(1 - (1 - conf)/2)
        
        res = cbind('Specialization' = 1-even,se)  %>% 
          mutate(Network = assemblage, Method = method ,Order.q = q, SC = SC,
                 "Spec.LCL" = Specialization - tmp*`s.e.`, "Spec.UCL"= Specialization + tmp*`s.e.`, class = paste0("1 - E",e), level = "Weighted species-level") %>%
          select(c("Order.q", 'Specialization',"s.e.", "Spec.UCL", "Spec.LCL","Method","SC","Network","class","level"))
        
        return(res)
      })%>%do.call("rbind",.)
    })
    
    names(Spec) = paste0("1 - E",E.class)
    
    Spec <- lapply(Spec, function(x) x %>% mutate('SC' = SC))
    for(i in 1:length(E.class)){
      Spec[[i]][,4:5] <- Spec[[i]][,c(5,4)]
      names(Spec[[i]])[4:5] <- c("Spec.LCL", "Spec.UCL")
      names(Spec[[i]])[8:9] <- c("Dataset", "Measure")
      Spec[[i]]$Method <- NULL
    }
    
  }else if(level == "network"){
    
    long = lapply(data, function(da){da%>%as.data.frame()%>%gather(key = "col_sp", value = "abundance")%>%.[,2]})
    
    Spec <- lapply(E.class, function(e){
        each_class = lapply(seq_along(long), function(i){
          res = iNEXT.4steps::Evenness(long[[i]], q = q,datatype = datatype,
                                       method = method, nboot=nboot, E.class = e, SC = SC)
          res['Coverage'] = NULL
          res = lapply(res, function(each_class){
            each_class%>%
              mutate(Evenness = 1-Evenness, Even.LCL = 1-Even.LCL, Even.UCL = 1-Even.UCL) %>% 
              select(-Assemblage)%>%
              rename('Specialization'='Evenness', 'Spec.UCL' ='Even.LCL', 'Spec.LCL' ='Even.UCL')%>%
              mutate(Network = names(long)[[i]],class = paste0("1 - E",e), level = "Network-level")
          })
          # if(method == "Observed") index = 1
          # if(method == "Estimated") index = 2
          # return(res[[index]]%>%mutate(Assemblage = names(long)[[i]]))
          return(res[[1]])
        }) %>% do.call("rbind",.)
        
        # each_class %>% mutate(class = paste0("1 - E",e))
      })
      names(Spec) = paste0("1 - E",E.class)
    
      for(i in 1:length(E.class)){
        Spec[[i]][,4:5] <- Spec[[i]][,c(5,4)]
        names(Spec[[i]])[4:5] <- c("Spec.LCL", "Spec.UCL")
        names(Spec[[i]])[8:9] <- c("Dataset", "Measure")
        Spec[[i]]$Method <- NULL
      }
  }
  
  return(Spec)
    
}


# ggSpec.link -------------------------------------------------------------------
#' ggplot2 extension for output from \code{Spec.link.est} and \code{Spec.link.ObsAsy}.
#' 
#'  \code{ggSpec.link}: Plots q-profile based on the output of \code{Spec.link.est} and \code{Spec.link.ObsAsy} using the ggplot2 package
#' 
#' @param output a table generated from Specialization function
#' @return a figure of estimated sample completeness with order q
#'
#' @examples
#' data(beetles_plotA)
#' output_spec = Spec.link.est(beetles_plotA, nboot = 5)
#' ggSpec.link(output_spec)
#' @export

ggSpec.link = function (output)
{

  classdata = cbind(do.call(rbind, output))
  
  # Check if the number of unique 'Network' is 8 or less
  if (length(unique(classdata$Dataset)) <= 8){
    cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  }else{
    # If there are more than 8 networks, start with the same predefined color palette
    # Then extend the palette by generating additional colors using the 'ggplotColors' function
    cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#330066", "#CC79A7", "#0072B2", "#D55E00"))
    cbPalette <- c(cbPalette, ggplotColors(length(unique(classdata$Dataset))-8))
  }
  
  fig = ggplot(classdata, aes(x = Order.q, y = Specialization, colour = Dataset)) +
    geom_line(size = 1.2) + 
    geom_ribbon(aes(ymin = Spec.LCL, ymax = Spec.UCL, fill = Dataset), alpha = 0.2, linetype = 0) +
    scale_colour_manual(values = cbPalette) + 
    scale_fill_manual(values = cbPalette) +
    labs(x = "Order q", y = "Specialization") + 
    theme_bw() + 
    theme(legend.position = "bottom", legend.box = "vertical", 
          legend.key.width = unit(1.2, "cm"), 
          legend.title = element_blank(), 
          legend.margin = margin(0, 0, 0, 0), 
          legend.box.margin = margin(-10, -10, -5, -10), 
          text = element_text(size = 16), plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) + 
    guides(linetype = guide_legend(keywidth = 2.5))
  if (length(output) != 1)
    fig = fig + facet_wrap(~Measure) + 
    theme(strip.text.x = element_text(size = 12, colour = "purple", face = "bold"))
  
  return(fig)
}


