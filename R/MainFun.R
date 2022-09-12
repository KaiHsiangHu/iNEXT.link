# DataInfo.link ----------------------
#' Exhibit basic data information
#'
#' \code{DataInfo.link}: exhibits basic data information
#'
#' @param data  a \code{list} of \code{data.frames}, each \code{data.frames} represents col.species-by-row.species abundance matrix.
#' @param diversity selection of diversity type: \code{'TD'} = 'Taxonomic diversity', \code{'PD'} = 'Phylogenetic diversity', and \code{'FD'} = 'Functional diversity'.
#' @return
#' a data.frame of basic data information incliuding sample size, observed species richness, sample coverage estimate, and the first ten abundance frequency counts.
#' @examples
#' #' ## Taxonomic diversity
#' data(beetles)
#' DataInfo.link(data = beetles, diversity = 'TD')
#'
#'
#' ## Phylogenetic diversity
#' data(beetles)
#' data(beetles_col_tree)
#' DataInfo.link(data = beetles, diversity = 'PD', col.tree = beetles_col_tree)
#'
#'
#' ## Functional diversity
#' data(beetles)
#' data(beetles_col_distM)
#' DataInfo.link(data = beetles, diversity = 'FD', col.distM = beetles_col_distM)
#' @export
DataInfo.link <- function(data, diversity = 'TD', datatype = "abundance", row.tree = NULL, col.tree = NULL, row.distM = NULL, col.distM = NULL){

  if(diversity == 'PD'){

    if(!is.null(row.tree)){row.tree$tip.label = gsub('\\.', '_',row.tree$tip.label)}
    if(!is.null(col.tree)){col.tree$tip.label = gsub('\\.', '_',col.tree$tip.label)}

    table <- lapply(data, function(y){datainfphy(data = y, datatype = datatype,
                                                 row.tree = row.tree,col.tree = col.tree)})%>%
      do.call(rbind,.)
    rownames(table) <- names(data)
    table = tibble::rownames_to_column(table, var = "Networks")
  }else if(diversity == 'TD'){
    table <- lapply(data, function(y){datainf(data = y, datatype = datatype)})%>%do.call(rbind,.)
    rownames(table) <- names(data)
    table = tibble::rownames_to_column(table, var = "Networks")
  }else if(diversity == 'FD'){


    table <- lapply(data, function(y){datainffun(data = y, datatype = datatype,
                                                 row.distM = row.distM,col.distM = col.distM)})%>%
      do.call(rbind,.)
    rownames(table) <- names(data)
    table = tibble::rownames_to_column(table, var = "Networks")
  }
  return(table)

}

# Completeness.link ----
#' Sample Completeness main function
#'
#' \code{Completeness.link} Estimation of Sample Completeness with order q
#'
#' @param data a \code{list} of \code{data.frames}, each \code{data.frames} represents col.species-by-row.species abundance matrix.
#' @param q q a numerical vector specifying the diversity orders. Default is \code{seq(0, 2, 0.2)}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing
#' sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter 0 to skip the bootstrap procedures. Default is \code{30}.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is \code{0.95}.
#' @return a matrix of estimated sample completeness with order q: \cr\cr
#'
#' @examples
#' data(beetles)
#' output = Completeness.link(beetles)
#' output
#'
#' @references
#' Chao,A.,Y.Kubota,D.Zelen??,C.-H.Chiu.
#' Quantifying sample completeness and comparing diversities among assemblages.
#' @export
Completeness.link <- function(data, q = seq(0, 2, 0.2), datatype = "abundance", nboot = 30, conf = 0.95){
  data_long <- lapply(data, function(tab){
    as.matrix(tab)%>%c()}
  )
  res = iNEXT.4steps::Completeness(data = data_long, q = q, datatype = datatype, nboot = nboot, conf = conf)
  return(res)
}

# ggCompleteness.link -------------------------------------------------------------------
#' ggplot for Sample Completeness
#'
#' \code{ggCompleteness.link} The figure for estimation of Sample Completeness with order q
#'
#' @param outcome a table generated from Completeness.link function
#' @return a figure of estimated sample completeness with order q
#'
#' @examples
#' data(beetles)
#' output = Completeness.link(beetles)
#' ggCompleteness.link(output)
#' @export
ggCompleteness.link <- function(outcome){
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  ggplot(outcome, aes(x = Order.q, y = Estimate.SC, colour = Assemblage)) +
    geom_line(size = 1.2) + scale_colour_manual(values = cbPalette) +
    geom_ribbon(aes(ymin = SC.LCL, ymax = SC.UCL, fill = Assemblage),
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
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is \code{0.95}.
#' @param col.tree (required only when \code{diversity = "PD"}), a phylogenetic tree of column assemblage in the pooled network column assemblage.
#' @param row.tree (required only when \code{diversity = "PD"}), a phylogenetic tree of row assemblage in the pooled network row assemblage.
#' @param PDtype (required only when \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"}(effective total branch length) or
#' \code{PDtype = "meanPD"}(effective number of equally divergent lineages).Default is \code{"meanPD"}.
#' @param col.distM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species of column assemblage in the pooled network column assemblage.
#' @param row.distM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species of row assemblage in the pooled network row assemblage.
#' @param FDtype (required only when \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_values"} for FD under specified threshold values, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{"AUC"}.
#' @param FDtau (required only when \code{diversity = "FD"} and \code{FDtype = "tau_values"}), a numerical vector between 0 and 1 specifying tau values (threshold levels). If \code{NULL} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy).
#' @import ape
#' @import dplyr
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
#' @import iNEXT.Beta3D
#' @import future
#' @import future.apply
#' @import ade4
#' @import tidyr
#' @import tibble
#' @import reshape2
#' @import sets
#' @return

#' \itemize{
#'  \item{\code{$DataInfo}: A dataframe summarizing data information}
#'  \item{\code{$iNextEst}: coverage-based diversity estimates along with confidence intervals}
#'  for showing diversity estimates for rarefied and extrapolated samples along with related statistics;
#'  \item{\code{$AsyEst}: for
#' showing asymptotic diversity estimates along with related statistics.}
#' }
#'
#' @examples
#' ## Taxonomic diversity
#' data(beetles)
#' output1 = iNEXT.link(data = beetles, diversity = 'TD', q = c(0,1,2))
#' output1$DataInfo # showing basic data information.
#' output1$iNextEst # showing diversity estimates with rarefied and extrapolated.
#' output1$AsyEst # showing asymptotic diversity estimates.
#'
#'
#' ## Phylogenetic diversity
#' data(beetles)
#' data(beetles_col_tree)
#' output2 = iNEXT.link(data = beetles, diversity = 'PD', q = c(0,1,2), col.tree = beetles_col_tree)
#' output2
#'
#'
#' ## Functional diversity under single threshold
#' data(beetles)
#' data(beetles_col_distM)
#' output3 = iNEXT.link(data = beetles, diversity = 'FD', q = c(0,1,2), col.distM = beetles_col_distM, FDtype = "tau_values")
#' output3
#'
#'
#' ## Functional diversity with thresholds integrating from 0 to 1
#' data(beetles)
#' data(beetles_col_distM)
#' output4 = iNEXT.link(data = beetles, diversity = 'FD', q = c(0,1,2), col.distM = beetles_col_distM, FDtype = "AUC")
#' output4
#'
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. \emph{Philosophical Transactions of the Royal Society B.}, 365, 3599-3609. \cr\cr
#' Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015). Rarefaction and extrapolation of phylogenetic diversity. \emph{Methods in Ecology and Evolution}, 6, 380-388.\cr\cr
#' Chao, A., Chiu C.-H. and Jost L. (2016). Phylogenetic diversity measures and their decomposition: a framework based on Hill numbers. pp. 141-172 in Pellens R. and Grandcolas P. (eds)
#' \emph{Biodiversity Conservation and Phylogenetic Systematics: Preserving our Evolutionary Heritage in an Extinction Crisis}, Springer. \cr\cr
#' Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. \emph{Systematic Biology}, 66, 100-111.
#' @export

iNEXT.link <- function(data, diversity = 'TD', q = c(0,1,2), datatype = "abundance", size = NULL, nT = NULL,
                       endpoint = NULL, knots = 40, conf = 0.95, nboot = 30,
                       row.tree = NULL, col.tree = NULL, PDtype = 'meanPD', row.distM = NULL, col.distM = NULL,
                       FDtype = "AUC", FDtau = NULL
){
  # User interface
  TYPE <- c("abundance", "incidence", "incidence_freq", "incidence_raw")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  class_x <- class(data)[1]

  if(datatype == "incidence"){
    stop('datatype="incidence" was no longer supported after v2.0.8,
         please try datatype="incidence_freq".')
  }
  if(datatype=="incidence_freq") datatype <- "incidence"

  if(datatype=="incidence_raw"){
    if(class_x=="list"){
      data <- lapply(data, as.incfreq)
    }else{
      data <- as.incfreq(data)
    }
    datatype <- "incidence"
  }
  if ( sum(!(diversity %in% c('TD', 'PD', 'FD', 'AUC')))>0 ){stop("Please select one of below diversity: 'TD', 'PD', 'FD'",
                                                                  call. = FALSE)}

  res = list()
  if(diversity == 'TD'){
    ## 1. datainfo
    datainfo = DataInfo.link(data = data, diversity = diversity, datatype = datatype)
    ## 2. iNterpolation/ Extrapolation
    data_long <- lapply(data, function(tab){
      as.matrix(tab)%>%c()}
    )
    INEXT_est <- iNEXT.3D::iNEXT3D(data_long, diversity = 'TD', q = q,conf = conf,
                                   nboot = nboot, knots = knots, endpoint = endpoint, size = size)

    res[[1]] = datainfo
    res[[2]] = INEXT_est$iNextEst
    res[[3]] = INEXT_est$AsyEst
    names(res) = c("DataInfo", "iNextEst", "AsyEst")

  }else if(diversity == 'PD'){

    if(!is.null(row.tree)){row.tree$tip.label = gsub('\\.', '_',row.tree$tip.label)}
    if(!is.null(col.tree)){col.tree$tip.label = gsub('\\.', '_',col.tree$tip.label)}
    res = iNEXTPDlink(data, q = q, datatype = datatype, size = size,
                      endpoint = endpoint, knots = knots, conf = conf,
                      nboot = nboot,col.tree = col.tree,row.tree = row.tree, type = PDtype)

  }else if (diversity == "FD" & FDtype == "tau_values") {

    res = iNEXTlinkFD(data, q = q, datatype = datatype, size = size,
                      endpoint = endpoint, knots = knots, conf = conf,
                      nboot = nboot, nT = nT, row.distM = row.distM, col.distM = col.distM, threshold = FDtau)
  }
  else if (diversity == "FD" & FDtype == "AUC") {
    res = iNEXTlinkAUC(data, q = q, datatype = datatype, size = size,
                       endpoint = endpoint, knots = knots, conf = conf,
                       nboot = nboot, nT = nT, row.distM = row.distM, col.distM = col.distM)
  }

  return(res)
}





# ggiNEXT.link -------------------------------------------------------------------
#' ggplot2 extension for outcome from \code{iNEXT.link}
#'
#' \code{ggiNEXT.link}: the \code{\link[ggplot2]{ggplot}} extension for \code{\link{iNEXT.link}} Object to plot sample-size-based and coverage-based rarefaction/extrapolation curves along with a bridging sample completeness curve
#' @param outcome a list object computed by \code{\link{iNEXT.link}}.
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1});
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).
#' @param se a logical variable to display confidence interval around the estimated sampling curve.
#' @param facet.var create a separate plot for each value of a specified variable:
#'  no separation \cr (\code{facet.var="None"});
#'  a separate plot for each diversity order (\code{facet.var="Order.q"});
#'  a separate plot for each assemblage (\code{facet.var="Assemblage"});
#'  a separate plot for each combination of order x assemblage (\code{facet.var="Both"}).
#' @param color.var create curves in different colors for values of a specified variable:
#'  all curves are in the same color (\code{color.var="None"});
#'  use different colors for diversity orders (\code{color.var="Order.q"});
#'  use different colors for sites (\code{color.var="Assemblage"});
#'  use different colors for combinations of order x assemblage (\code{color.var="Both"}).
#' @param grey a logical variable to display grey and white ggplot2 theme.
#' @param ... other arguments passed on to methods. Not currently used.
#' @return a ggplot2 object
#' @examples
#' ## Taxonomic diversity
#' data(beetles)
#' output1 = iNEXT.link(data = beetles, diversity = 'TD', q = c(0,1,2))
#' ggiNEXT.link(output1, diversity = 'TD')
#'
#'
#' ## Phylogenetic diversity
#' data(beetles)
#' data(beetles_col_tree)
#' output2 = iNEXT.link(data = beetles, diversity = 'PD', q = c(0,1,2), col.tree = beetles_col_tree)
#' ggiNEXT.link(output2, diversity = 'PD')
#'
#'
#' ## Functional diversity under single threshold
#' data(beetles)
#' data(beetles_col_distM)
#' output3 = iNEXT.link(data = beetles, diversity = 'FD', q = c(0,1,2), col.distM = beetles_col_distM, FDtype = "tau_values")
#' ggiNEXT.link(output3, diversity = 'FD')
#'
#'
#' ## Functional diversity with thresholds integrating from 0 to 1
#' data(beetles)
#' data(beetles_col_distM)
#' output4 = iNEXT.link(data = beetles, diversity = 'FD', q = c(0,1,2), col.distM = beetles_col_distM, FDtype = "AUC")
#' ggiNEXT.link(output4, diversity = 'FD')
#' @export
ggiNEXT.link <- function(outcome, diversity = 'TD', type = c(1,2,3), se = TRUE, facet.var = "Assemblage",
                         color.var = "Order.q"){
  if(diversity == 'TD'){
    res = iNEXT.3D::ggiNEXT3D(outcome, type = c(1,2,3),facet.var = facet.var)
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
    res = iNEXT.3D::ggiNEXT3D(outcome, type = type,facet.var = facet.var)
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
    # # output = outcome
    # # output$iNextEst$size_based = output$iNextEst$size_based%>%
    # #   rename('qD'="PD", 'qD.UCL'="PD.UCL",'qD.LCL'="PD.LCL")
    # # iNEXT.3D::ggiNEXT3D(output, type = 1)
    # iNE <- outcome$iNextEst
    # iNE.sub <- iNE[iNE$method == "observed",]
    # iNE[iNE$method == "observed",]$method <-  "Rarefaction"
    # ex <- iNE.sub
    # ex$method <- "Extrapolation"
    # iNE <- rbind(iNE,ex)
    # iNE$method <- factor(iNE$method,levels = c("Rarefaction","Extrapolation"))
    # iNE$Order.q = paste0("q = ", iNE$Order.q)
    # iNE.sub$Order.q = paste0("q = ", iNE.sub$Order.q)
    #
    # if(type == 1){
    #   # size-based
    #   plot <- ggplot(iNE, aes(x = m,y = PD)) + geom_line(aes(color = Region,linetype = method),size = 1.2) +
    #     facet_wrap(~Order.q, scales = "free") +
    #     # facet_grid()
    #     geom_ribbon(aes(x = m,ymax = PD.UCL ,ymin = PD.LCL,fill = Region),alpha = 0.25) + theme_bw()+
    #     geom_point(aes(x = m,y = PD ,color = Region,shape = Region),size = 5,data = iNE.sub) +
    #     xlab("Number of individuals")
    # }else if(type == 2){
    #   # output <- outcome$size_based
    #   # if (length(unique(output$Order.q)) > 1) output <- subset(output, Order.q == unique(output$Order.q)[1])
    #   # output$y.lwr <- output$SC.LCL
    #   # output$y.upr <- output$SC.UCL
    #   # id <- match(c(x_name, "Method", "SC", "SC.LCL", "SC.UCL", "Assemblage", "Order.q", "qD", "qD.LCL", "qD.UCL"), names(output), nomatch = 0)
    #   # output[,1:10] <- output[, id]
    #   #
    #   # xlab_name <- paste0("Number of ", xlab_name)
    #   # ylab_name <- "Sample Coverage"
    #
    # }else if(type == 3){
    #   # coverage-based
    #   plot <- ggplot(iNE) + geom_line(aes(x = SC,y = PD,color = Region,linetype = method),size = 1.2) +
    #     facet_wrap(~Order.q, scales = "free") +
    #     geom_ribbon(aes(x = SC,ymax = PD.UCL ,ymin = PD.LCL,fill = Region),alpha = 0.25) +
    #     geom_point(aes(x = SC,y = PD ,color = Region,shape = Region),size = 5,data = iNE.sub) + theme_bw() +
    #     xlab("Sample coverage")
    # }
    #
    # plot +
    #   theme(legend.position = "bottom",
    #         legend.title=element_blank(), strip.text = element_text(size = stript.size),
    #         text=element_text(size=text.size),
    #         legend.key.width = unit(0.8,"cm"))  +
    #   labs(y = "Phylogenetic network diversity", lty = "Method")
    # # if(grey){
    # #   g <- g +
    # #     scale_fill_grey(start = 0, end = .4) +
    # #     scale_colour_grey(start = .2, end = .2)
    # # }


  }else if(diversity == 'FD'){
    res = iNEXT.3D::ggiNEXT3D(outcome, type = type,facet.var = facet.var)
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


# AO.link -------------------------------------------------------------------
#' Asymptotic and Observed diversity q profile
#'
#' \code{AO.link}: The asymptotic (or observed) diversity of order q
#'
#' @param data a \code{list} of \code{data.frames}, each \code{data.frames} represents col.species-by-row.species abundance matrix.
#' @param diversity selection of diversity type: \code{'TD'} = Taxonomic diversity, \code{'PD'} = Phylogenetic diversity, and \code{'FD'} = Functional diversity.
#' @param q a numerical vector specifying the diversity orders. Default is \code{c(0,1,2)}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter 0 to skip the bootstrap procedures. Default is \code{30}.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is \code{0.95}.
#' @param method asymptotic or Observed
#' @param col.tree (required only when \code{diversity = "PD"}), a phylogenetic tree of column assemblage in the pooled network column assemblage.
#' @param row.tree (required only when \code{diversity = "PD"}), a phylogenetic tree of row assemblage in the pooled network row assemblage.
#' @param PDtype (required only when \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"}(effective total branch length) or \code{PDtype = "meanPD"}(effective number of equally divergent lineages).Default is \code{"meanPD"}.
#' @param col.distM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species of column assemblage in the pooled network column assemblage.
#' @param row.distM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species of row assemblage in the pooled network row assemblage.
#' @param FDtype (required only when \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_values"} for FD under specified threshold values, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{"AUC"}.
#' @param FDtau (required only when \code{diversity = "FD"} and \code{FDtype = "tau_values"}), a numerical vector between 0 and 1 specifying tau values (threshold levels). If \code{NULL} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy).
#' @return a table of Asymptotic network diversity q profile.
#'
#' @examples
#' ## Taxonomic diversity
#' data(beetles)
#' output1 = AO.link(data = beetles, diversity = 'TD', q = seq(0, 2, 0.2))
#' output1
#'
#'
#' ## Phylogenetic diversity
#' data(beetles)
#' data(beetles_col_tree)
#' output2 = AO.link(data = beetles, diversity = 'PD', q = seq(0, 2, 0.2), col.tree = beetles_col_tree)
#' output2
#'
#'
#' ## Functional diversity under single threshold
#' data(beetles)
#' data(beetles_col_distM)
#' output3 = AO.link(data = beetles, diversity = 'FD', q = seq(0, 2, 0.2), col.distM = beetles_col_distM, FDtype = "tau_values")
#' output3
#'
#'
#' ## Functional diversity with thresholds integrating from 0 to 1
#' data(beetles)
#' data(beetles_col_distM)
#' output4 = AO.link(data = beetles, diversity = 'FD', q = seq(0, 2, 0.25),
#'                   col.distM = beetles_col_distM, FDtype = "AUC", nboot = 0)
#' output4
#' @export
AO.link <- function(data, diversity = 'TD', q = seq(0, 2, 0.2), datatype = "abundance", nboot = 30, conf = 0.95, method = c("Asymptotic", "Observed"),
                    row.tree = NULL, col.tree = NULL, PDtype = "meanPD", row.distM = NULL, col.distM = NULL, FDtype = "AUC", FDtau = NULL){

  if ( !(diversity %in% c('TD', 'PD', 'FD')) )
    stop("Please select one of below diversity: 'TD', 'PD', 'FD'", call. = FALSE)

  if (diversity == 'TD') {

    if (sum(method == 'Asymptotic') == length(method))
      NetDiv <- AsylinkTD(data, diversity = 'TD', q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95) else if (sum(method == 'Observed') == length(method))

        NetDiv <- ObslinkTD(data, diversity = 'TD', q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95) else if (sum(method == c('Asymptotic', 'Observed')) == length(method))

          NetDiv = rbind(AsylinkTD(data, diversity = 'TD', q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95),
                         ObslinkTD(data, diversity = 'TD', q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95))

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
    NetDiv = NetDiv %>%
      dplyr::select(Order.q, Estimate, s.e., LCL, UCL, Method, Region, Reftime, Type ) %>%
      set_colnames(c('Order.q', 'qD', "s.e.",'qD.LCL','qD.UCL', 'Method',"Network", "Reftime","Type"))
    
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

  }

  return(NetDiv)

}


# ggAO.link -------------------------------------------------------------------
#' ggplot for Asymptotic Network diversity
#'
#' \code{ggAO.link} Plots q-profile based on the outcome of \code{AO.link} using the ggplot2 package.\cr
#'
#' @param outcome the outcome of the functions \code{AO.link} .\cr
#' @param diversity diversity type
#' @param text.size control the text size of the output plot.
#' @return a figure of asymptotic (or observed) diversity with order q\cr\cr
#'
#' @examples
#' ## Taxonomic diversity
#' data(beetles)
#' output1 = AO.link(data = beetles, diversity = 'TD', q = seq(0, 2, 0.2))
#' ggAO.link(output1, diversity = 'TD')
#'
#'
#' ## Phylogenetic diversity
#' data(beetles)
#' data(beetles_col_tree)
#' output2 = AO.link(data = beetles, diversity = 'PD', q = seq(0, 2, 0.2), col.tree = beetles_col_tree)
#' ggAO.link(output2, diversity = 'PD')
#'
#'
#' ## Functional diversity under single threshold
#' data(beetles)
#' data(beetles_col_distM)
#' output3 = AO.link(data = beetles, diversity = 'FD', q = seq(0, 2, 0.2), col.distM = beetles_col_distM, FDtype = "tau_values")
#' ggAO.link(output3, diversity = 'FD')
#'
#'
#' ## Functional diversity with thresholds integrating from 0 to 1
#' data(beetles)
#' data(beetles_col_distM)
#' output4 = AO.link(data = beetles, diversity = 'FD', q = seq(0, 2, 0.25), col.distM = beetles_col_distM, FDtype = "AUC", nboot = 0)
#' ggAO.link(output4, diversity = 'FD')
#' @export
ggAO.link <- function(outcome, diversity = 'TD', text.size = 14){
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  # if (sum(unique(outcome$method) %in% c("Estimate", "Observed")) == 0)
  #   stop("Please use the outcome from specified function 'AsyD'")
  if(diversity %in% c('TD','PD')){
    if(diversity == 'TD'){
      ylab = 'Taxonomic network diversity'
    }else if(diversity == 'PD') {
      ylab = 'Phylogenetic network diversity'
    }
    ggplot(outcome, aes(x = Order.q, y = qD, colour = Network, lty = Method)) +
      geom_line(size = 1.2) + scale_colour_manual(values = cbPalette) +
      geom_ribbon(data = outcome[outcome$Method == "Asymptotic", ],
                  aes(ymin = qD.LCL, ymax = qD.UCL, fill = Network), alpha = 0.2, linetype = 0) +
      scale_fill_manual(values = cbPalette) +
      scale_linetype_manual(values = c(Asymptotic = 1, Empirical = 2)) +
      labs(x = "Order q", y = ylab) + theme_bw() + theme(text = element_text(size = 10)) +
      theme(legend.position = "bottom", legend.box = "vertical",
            legend.key.width = unit(1.1, "cm"), legend.title = element_blank(),
            legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-10,-10, -5, -10),
            text = element_text(size = text.size)
      )

    # +labs(x = "Order q", y = "Network phylogenetic diversity", lty = "Method") + scale_linetype_manual(values=c("dashed","solid"))

  }else if(diversity == 'FD'){
    iNEXT.3D:::ggAO3D(outcome)+ylab("Functional network diversity")
  }

}



# estimateD.link  -------------------------------------------------------------------
#' Compute species diversity with a particular of sample size/coverage
#'
#' \code{estimateD.link} computes species diversity (Hill numbers with q = 0, 1 and 2) with a particular user-specified level of sample size or sample coverage.
#'
#' @param data a \code{matrix}, \code{data.frame} (species by assemblages), or \code{list} of species abundance/incidence raw data.\cr
#' @param diversity a choice of three-level diversity: 'TD' = 'Taxonomic', 'PD' = 'Phylogenetic', and 'FD' = 'Functional' under certain threshold. Besides,'AUC' is the fourth choice which
#' integrates several threshold functional diversity to get diversity.
#' @param q a numerical vector of the order of Hill number.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}).
#' @param base comparison base: sample-size-based (\code{base="size"}) or coverage-based \cr (\code{base="coverage"}).
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing
#' sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 50
#' @param level a sequence specifying the particular sample sizes or sample coverages(between 0 and 1).
#' If \code{base="size"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample size among all sites extrapolated to double reference sizes.
#' If \code{base="coverage"} and \code{level=NULL}, then this function computes the diversity estimates for the minimum sample coverage among all sites extrapolated to double reference sizes.
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
#' @param PDtype (required only when \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"}(effective total branch length) or
#' \code{PDtype = "meanPD"}(effective number of equally divergent lineages).Default is \code{"meanPD"}.
#' @param col.tree phylogenetic tree of column assemblage in interaction matrix.
#' @param row.tree phylogenetic tree of row assemblage in interaction matrix.
#' @param col.distM (required only when diversity = "FD"), a column species pairwise distance matrix for all column species of column assemblage in interaction matrix.
#' @param row.distM (required only when diversity = "FD"), a row species pairwise distance matrix for all row species of row assemblage in interaction matrix.
#' @param FDtype (required only when diversity = "FD"), select FD type: FDtype = "tau_values" for FD under specified threshold values, or FDtype = "AUC" (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is "AUC".
#' @param FDtau (required only when diversity = "FD" and FDtype = "tau_values"), a numerical vector between 0 and 1 specifying tau values (threshold levels). If NULL (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy).
#' @return a data.frame of species diversity table including the sample size, sample coverage, method (rarefaction or extrapolation), and diversity estimates with q = 0, 1, and 2 for the user-specified sample size or sample coverage.
#'
#' @examples
#' \dontrun{
#' ## Taxonomic diversity
#' data(beetles)
#' output1 <- estimateD.link(beetles, diversity = 'TD', datatype = "abundance",
#'                           base = "coverage", level = 0.7, nboot = 30)
#' output1
#' 
#' ## Phylogenetic diversity
#' output2 <- estimateD.link(beetles, diversity = 'PD', datatype = "abundance",
#'                           base = "size", level = NULL, nboot = 30, col.tree = beetles_col_tree)
#' output2
#' 
#' ## Functional diversity under single threshold
#' data(beetles)
#' data(beetles_col_distM)
#' output3 = estimateD.link(data = beetles, diversity = 'FD', col.distM = beetles_col_distM, FDtype = "tau_values")
#' output3
#'
#'
#' ## Functional diversity with thresholds integrating from 0 to 1
#' data(beetles)
#' data(beetles_col_distM)
#' output4 = estimateD.link(data = beetles, diversity = 'FD',
#'                          col.distM = beetles_col_distM, FDtype = "AUC", nboot = 0)
#' output4
#' 
#' }
#' @export
estimateD.link = function(data, diversity = 'TD', q = c(0, 1, 2), datatype = "abundance", base = "coverage",
                          level = NULL, nboot = 50, conf = 0.95, PDtype = 'meanPD',
                          row.tree = NULL, col.tree = NULL, row.distM = NULL, col.distM = NULL, FDtype = "AUC", FDtau = NULL){
  if(diversity == 'TD'){

    div = lapply(1:length(data), function(i){
      x = data[[i]]
      assemblage = names(data)[[i]]
      long = as.matrix(x)%>%c()
      iNEXT.3D::estimate3D(long, q=q,datatype=datatype, base=base,
                           diversity = 'TD', nboot = nboot,conf=conf,level = level)%>%
        mutate(Assemblage = assemblage)
    })%>%do.call("rbind",.)

    return(div)
  }else if(diversity == 'PD'){

    if(datatype=='abundance'){

      if(class(data)=="data.frame" | class(data)=="matrix"| class(data)=="integer" ) data = list(Region_1 = data)

      if(class(data)== "list"){
        if(is.null(names(data))) region_names = paste0("Region_", 1:length(data)) else region_names = names(data)
        Ns = sapply(data, ncol)
        data_list = data
      }

    }
    if(is.null(conf)) conf = 0.95
    tmp = qnorm(1 - (1 - conf)/2)

    for_each_region = function(data_2d, region_name, N){
      if (datatype=='abundance') {
        
        n = sum(data_2d)
        if(base == 'coverage'){
          size_m = sapply(level, function(i) coverage_to_size(data_2d, i, datatype='abundance'))
        }else if(base == 'size'){
          if(is.null(level)){
            size_m = n
          }else{
            size_m = level
          }
          level = iNEXT.3D:::Coverage(data_2d,m= n, datatype = 'abundance')

        }


        ref= iNEXT.3D:::Coverage(data_2d,m= n, datatype = 'abundance')
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
          PD.sd <- lapply(boot.sam, function(aL_boot){

            tmp = iNEXT.3D:::PhD.m.est(ai = aL_boot$branch.abun,
                                       Lis = aL_boot$branch.length%>%as.matrix(),
                                       m = size_m,
                                       q = q,nt = n,cal = PDtype)%>%
              as.vector()%>%as.data.frame()
            return(tmp)
          })%>%
            abind(along=3) %>% apply(1:2, sd)%>%as.vector()
        }else{
          PD.sd = rep(NA, length(qPDm))
        }
        
        
        ##
        len = length(q)
        res = data.frame(Assemblage = rep(region_name,len),
                         SC = rep(level, rep(len,length(size_m))),
                         m = rep(size_m,rep(len,length(size_m))),
                         Method = ifelse(level > ref, 'Extrapolation', 'Rarefaction'),
                         Order.q = q,
                         qPD = qPDm,
                         qPD.LCL = qPDm-1.96*PD.sd,
                         qPD.UCL = qPDm+1.96*PD.sd,
                         Reftime = tbar, 
                         Type = PDtype
        )
        return(res)
      }
    }

    output = lapply(1:length(data), function(i) for_each_region(data_2d = data_list[[i]],
                                                                region_name = region_names[i], N = Ns[i]))%>%
      do.call('rbind',.)

    return(output)
  }else if(diversity == 'FD'& FDtype == 'tau_values'){
    output = estimatelinkFD(data, row.distM = row.distM, col.distM = col.distM, datatype = datatype, q = q,
                            base = base, threshold = FDtau, level = level, nboot = nboot,
                            conf = conf)
    return(output)
  }else if(diversity == 'FD'& FDtype == 'AUC'){
    output = estimatelinkAUC(data, row.distM = row.distM, col.distM = col.distM, datatype = datatype, q = q,
                             base = base, level = level, nboot = nboot,
                             conf = conf)
    return(output)
  }
}




# iNEXTBeta.link ---------------------------
#' Interpolation (rarefaction) and extrapolation of network Beta diversity

#' Function \code{iNEXTBeta.link} Interpolation and extrapolation of Beta diversity with order q
#'
#' @param data data can be input as a \code{lisst} of \code{data.frame}, each \code{data.frame} represents col.species-by-row.species abundance matrix; see Note 1 for an example.
#' @param diversity selection of diversity type: \code{'TD'} = 'Taxonomic diversity', \code{'PD'} = 'Phylogenetic diversity', and \code{'FD'} = 'Functional diversity'.
#' @param q a numerical vector specifying the diversity orders. Default is \code{c(0,1,2)}.
#' @param level a sequence specifying the particular sample coverages (between 0 and 1). Default is \code{seq(0.5, 1, 0.05)}.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing
#' sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter 0 to skip the bootstrap procedures. Default is \code{30}.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is \code{0.95}.
#' @param col.tree (required only when \code{diversity = "PD"}), a phylogenetic tree of column assemblage in the pooled network column assemblage.
#' @param row.tree (required only when \code{diversity = "PD"}), a phylogenetic tree of row assemblage in the pooled network row assemblage.
#' @param PDtype (required only when \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"}(effective total branch length) or
#' \code{PDtype = "meanPD"}(effective number of equally divergent lineages).Default is \code{"meanPD"}.
#' @param col.distM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species of column assemblage in the pooled network column assemblage.
#' @param row.distM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species of row assemblage in the pooled network row assemblage.
#' @param FDtype (required only when \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_values"} for FD under specified threshold values, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{"AUC"}.
#' @param FDtau (required only when \code{diversity = "FD"} and \code{FDtype = "tau_values"}), a numerical vector between 0 and 1 specifying tau values (threshold levels). If \code{NULL} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy).
#' @return A list of seven lists with three-diversity and four-dissimilarity.
#' @examples
#' ## Taxonomic diversity
#' data(beetles)
#' output1 = iNEXTBeta.link(data = beetles, diversity = 'TD', level = seq(0.5, 0.9, 0.4), q = c(0, 1, 2))
#' output1
#'
#' ## Phylogenetic diversity
#' data(beetles)
#' data(beetles_col_tree)
#' output2 = iNEXTBeta.link(data = beetles, diversity = 'PD', level = seq(0.5, 0.9, 0.4), q = c(0, 1, 2), col.tree = beetles_col_tree)
#' output2
#'
#'
#' ## Functional diversity under single threshold
#' data(beetles)
#' data(beetles_col_distM)
#' output3 = iNEXTBeta.link(data = beetles, diversity = 'FD', level = seq(0.5, 0.9, 0.4), q = c(0, 1, 2), col.distM = beetles_col_distM, FDtype = "tau_values")
#' output3
#'
#'
#' ## Functional diversity with thresholds integrating from 0 to 1
#' data(beetles)
#' data(beetles_col_distM)
#' output4 = iNEXTBeta.link(data = beetles, diversity = 'FD', level = seq(0.5, 0.9, 0.4), q = c(0, 1, 2), col.distM = beetles_col_distM, FDtype = "AUC")
#' output4
#' @references
#' Chao, A., Chazdon, R. L., Colwell, R. K. and Shen, T.-J.(2005). A new statistical approach for assessing similarity of species composition with incidence and abundance data. Ecology Letters 8, 148-159. (pdf file) Spanish translation in pp. 85-96 of Halffter, G. Soberon, J., Koleff, P. and Melic, A. (eds) 2005 Sobre Diversidad Biologica: el Sognificado de las Diversidades Alfa, Beta y Gamma. m3m-Monografias 3ercer Milenio, vol. 4, SEA, CONABIO, Grupo DIVERSITAS & CONACYT, Zaragoza. IV +242 pp.
#' Chiu, C.-H., Jost, L. and Chao*, A. (2014). Phylogenetic beta diversity, similarity, and differentiation measures based on Hill numbers. Ecological Monographs 84, 21-44.
#' @export

iNEXTBeta.link = function(data, diversity = 'TD', level = seq(0.5, 1, 0.05), datatype = 'abundance',
                          q = c(0, 1, 2), nboot = 20, conf = 0.95, PDtype = 'meanPD',
                          row.tree = NULL, col.tree = NULL, row.distM = NULL, col.distM = NULL,
                          FDtype = "AUC", FDtau = NULL, FDcut_number = 30){

  if(class(data[[1]]) == 'data.frame' ){dat = list(data); }else{dat = data}

  combined_list = lapply(dat, function(y){

    long = ready4beta(y)%>%filter_all(any_vars(. != 0))
    rownames(long) = rownames(long)%>%gsub('\\.','_',.)
    colnames(long) = names(y)
    return(long)
  })

  if(diversity == 'TD'){
    dissimilarity <- iNEXTBeta3D(data = combined_list, diversity = 'TD',level = level, datatype = datatype,
                                 q = q ,nboot = nboot, conf = conf)



  }
  else if(diversity == 'PD'){

    if(!is.null(row.tree)){row.tree$tip.label = gsub('\\.', '_',row.tree$tip.label)}
    if(!is.null(col.tree)){col.tree$tip.label = gsub('\\.', '_',col.tree$tip.label)}

    dissimilarity = iNEXTbeta.PDlink(data = combined_list, level = level, datatype = datatype,
                                     q =q ,row.tree = row.tree,col.tree = col.tree, nboot = nboot, PDtype = PDtype)
  }else if(diversity == 'FD' & FDtype == 'tau_values'){
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

    dissimilarity <- iNEXTBeta3D(data = combined_list, diversity = 'FD',level = level, datatype = datatype,
                                 q =q ,nboot = nboot, conf = conf, FDdistM = distM, FDtype = "tau_value", FDtau = FDtau)
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


    dissimilarity <- iNEXTBeta3D(data = combined_list, diversity = 'FD', level = level, datatype = datatype,
                                 q = q ,nboot = nboot, conf = conf, FDdistM = distM, FDcut_number = FDcut_number)
  }

  return(dissimilarity)
}


# ggiNEXTBeta.link -------------------------------------------------------------------
#' ggplot2 extension for outcome from \code{iNEXTBeta.link}
#'
#' \code{ggiNEXTBeta.link}: ggplot for Interpolation and extrapolation of Beta diversity with order q
#'
#' @param outcome the outcome from \code{"iNEXTBeta.link"}
#' @param type selection of plot type : \code{type = 'B'} for plotting the gamma, alpha, and beta diversity ;
#' \code{type = 'D'} for plotting 4 turnover dissimilarities.
#' @param scale Are scales shared across all facets (\code{"fixed"}), or do they vary across rows (\code{"free_x"}), columns (\code{"free_y"}), or both rows and columns (\code{"free"})? Default is \code{"free"}.
#'
#' @return a figure for Beta diversity or dissimilarity diversity.
#' @examples
#' ## Taxonomic diversity
#' data(beetles)
#' output1 = iNEXTBeta.link(data = beetles, diversity = 'TD', level = seq(0.5, 0.9, 0.4), q = c(0, 1, 2))
#' ggiNEXTBeta.link(output1, type = 'B')
#' ggiNEXTBeta.link(output1, type = 'D')
#'
#' ## Phylogenetic diversity
#' data(beetles)
#' data(beetles_col_tree)
#' output2 = iNEXTBeta.link(data = beetles, diversity = 'PD', level = seq(0.5, 0.9, 0.4), q = c(0, 1, 2), col.tree = beetles_col_tree)
#' ggiNEXTBeta.link(output2, type = 'B')
#' ggiNEXTBeta.link(output2, type = 'D')
#'
#'
#' ## Functional diversity under single threshold
#' data(beetles)
#' data(beetles_col_distM)
#' output3 = iNEXTBeta.link(data = beetles, diversity = 'FD', level = seq(0.5, 0.9, 0.4), q = c(0, 1, 2), col.distM = beetles_col_distM, FDtype = "tau_values")
#' ggiNEXTBeta.link(output3, type = 'B')
#' ggiNEXTBeta.link(output3, type = 'D')
#'
#'
#' ## Functional diversity with thresholds integrating from 0 to 1
#' data(beetles)
#' data(beetles_col_distM)
#' output4 = iNEXTBeta.link(data = beetles, diversity = 'FD', level = seq(0.5, 0.9, 0.4), q = c(0, 1, 2), col.distM = beetles_col_distM, FDtype = "AUC")
#' ggiNEXTBeta.link(output4, type = 'B')
#' ggiNEXTBeta.link(output4, type = 'D')
#' @export

ggiNEXTBeta.link <- function(outcome, type = c('B', 'D'), scale = 'free'){

  if (type == 'B'){

    gamma = lapply(outcome, function(y) y[["gamma"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Gamma") %>% as_tibble()
    alpha = lapply(outcome, function(y) y[["alpha"]]) %>% do.call(rbind,.) %>% mutate(div_type = "Alpha") %>% as_tibble()
    beta =  lapply(outcome, function(y) y[["beta"]])  %>% do.call(rbind,.) %>% mutate(div_type = "Beta")  %>% as_tibble()
    beta = beta %>% filter(Method != 'Observed')
    beta[beta == 'Observed_alpha'] = 'Observed'

    df = rbind(gamma, alpha, beta)
    for (i in unique(gamma$Order.q)) df$Order.q[df$Order.q==i] = paste0('q = ', i)
    df$div_type <- factor(df$div_type, levels = c("Gamma","Alpha","Beta"))

    id_obs = which(df$Method == 'Observed')

    for (i in 1:length(id_obs)) {

      new = df[id_obs[i],]
      new$SC = new$SC - 0.0001
      new$Method = 'Rarefaction'

      newe = df[id_obs[i],]
      newe$SC = newe$SC + 0.0001
      newe$Method = 'Extrapolation'

      df = rbind(df, new, newe)

    }

    if (unique(outcome[[1]]$gamma$diversity) == 'TD') { ylab = "Taxonomic diversity" }
    if (unique(outcome[[1]]$gamma$diversity) %in% c('PD','meanPD')) { ylab = "Phylogenetic diversity" }
    if (unique(outcome[[1]]$gamma$diversity) == 'FD_tau') { ylab = "Functional diversity" }
    if (unique(outcome[[1]]$gamma$diversity) == 'FD_AUC') { ylab = "Functional diversity (AUC)" }

  }
  if (type=='D'){

    C = lapply(outcome, function(y) y[["C"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-CqN") %>% as_tibble()
    U = lapply(outcome, function(y) y[["U"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-UqN") %>% as_tibble()
    V = lapply(outcome, function(y) y[["V"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-VqN") %>% as_tibble()
    S = lapply(outcome, function(y) y[["S"]]) %>% do.call(rbind,.) %>% mutate(div_type = "1-SqN") %>% as_tibble()
    C = C %>% filter(Method != 'Observed')
    U = U %>% filter(Method != 'Observed')
    V = V %>% filter(Method != 'Observed')
    S = S %>% filter(Method != 'Observed')
    C[C == 'Observed_alpha'] = U[U == 'Observed_alpha'] = V[V == 'Observed_alpha'] = S[S == 'Observed_alpha'] = 'Observed'

    df = rbind(C, U, V, S)
    for (i in unique(C$Order.q)) df$Order.q[df$Order.q==i] = paste0('q = ', i)
    df$div_type <- factor(df$div_type, levels = c("1-CqN","1-UqN","1-VqN","1-SqN"))

    id_obs = which(df$Method == 'Observed')

    for (i in 1:length(id_obs)) {

      new = df[id_obs[i],]
      new$SC = new$SC - 0.0001
      new$Method = 'Rarefaction'

      newe = df[id_obs[i],]
      newe$SC = newe$SC + 0.0001
      newe$Method = 'Extrapolation'

      df = rbind(df, new, newe)

    }
    if (unique(outcome[[1]]$gamma$diversity) == 'TD') { ylab = "Taxonomic dissimilarity" }
    if (unique(outcome[[1]]$gamma$diversity) %in% c('PD','meanPD')) { ylab = "Phylogenetic dissimilarity" }
    if (unique(outcome[[1]]$gamma$diversity) == 'FD_tau') { ylab = "Functional dissimilarity" }
    if (unique(outcome[[1]]$gamma$diversity) == 'FD_AUC') { ylab = "Functional dissimilarity (AUC)" }

  }
  lty = c(Rarefaction = "solid", Extrapolation = "dashed")
  # lty = c(Rarefaction = "solid", Extrapolation = "twodash")

  df$Method = factor(df$Method, levels = c('Rarefaction', 'Extrapolation', 'Observed'))

  double_size = unique(df[df$Method=="Observed",]$Size)*2
  double_extrapolation = df %>% filter(Method=="Extrapolation" & round(Size) %in% double_size)
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  point_size = 2

  ggplot(data = df, aes(x = SC, y = Estimate, col = Region)) +
    geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Region, col = NULL), alpha=0.4) +
    geom_line(data = subset(df, Method!='Observed'), aes(linetype=Method), size=1.1) +
    scale_linetype_manual(values = lty) +
    scale_fill_manual(values = cbPalette) + scale_colour_manual(values = cbPalette)+
    # geom_line(lty=2) +
    geom_point(data = subset(df, Method=='Observed' & div_type=="Gamma"),shape=19, size=point_size) +
    geom_point(data = subset(df, Method=='Observed' & div_type!="Gamma"),shape=1, size=point_size,stroke=1.5)+
    geom_point(data = subset(double_extrapolation, div_type == "Gamma"),shape=17, size=point_size) +
    geom_point(data = subset(double_extrapolation, div_type!="Gamma"),shape=2, size=point_size,stroke=1.5) +
    facet_grid(div_type~Order.q, scales = scale) +
    # facet_wrap(div_type~Order.q, scales = scale, switch="both") +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    labs(x='Sample coverage', y=ylab)
}


# Spec.link -------------------------------------------------------------------
#' Estimation (Observed) of Specialization with order q
#' @param data  a \code{list} of \code{data.frames}, each \code{data.frames} represents col.species-by-row.species abundance matrix.
#' @param q a numerical vector specifying the diversity orders. Default is \code{seq(0, 2, 0.2)}.
#' @param method abinary calculation method with 'Estimated' or 'Observed'.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing
#' sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter 0 to skip the bootstrap procedures. Default is \code{30}.
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is \code{0.95}.
#' @param E.class an integer vector between 1 to 5.
#' @param C a standardized coverage for calculating specialization index. It is used when \code{method = 'Estimated'}. If \code{NULL}, \code{C = Cmax}.
#' @return A list of estimated(Observed) specialization with order q.\cr
#' Different lists represents different classes of Specialization.\cr
#' Each list is combined with order.q and sites.\cr
#'
#' @examples
#' data(beetles)
#' output = Spec.link(beetles)
#' output
#' @export

Spec.link <- function(data, q = seq(0, 2, 0.2),
                      diversity = 'TD',
                      datatype = "abundance",
                      method = "Estimated",
                      nboot = 30,
                      conf = 0.95,
                      E.class = c(1:5),
                      C = NULL){
  if (diversity == 'TD'){
    long = lapply(data, function(da){da%>%as.data.frame()%>%gather(key = "col_sp", value = "abundance")%>%.[,2]})

    Spec <- lapply(E.class, function(e){
      each_class = lapply(seq_along(long), function(i){
        res = iNEXT.4steps::Evenness(long[[i]], q = q,datatype = datatype,
                                     method = method, nboot=nboot, E.class = e, C = C)
        res['Coverage'] = NULL
        res = lapply(res, function(each_class){
          each_class%>%
            mutate(Evenness = 1-Evenness, Even.LCL = 1-Even.LCL, Even.UCL = 1-Even.UCL)%>%
            rename('Specialization'='Evenness', 'Spec.LCL' ='Even.LCL', 'Spec.UCL' ='Even.UCL')%>%
            mutate(Network = names(long)[[i]])
        })
        # if(method == "Observed") index = 1
        # if(method == "Estimated") index = 2
        # return(res[[index]]%>%mutate(Assemblage = names(long)[[i]]))
        return(res[[1]])
      })%>%do.call("rbind",.)

      each_class%>%mutate(class = paste0("1 - E",e))
    })
    names(Spec) = paste0("1 - E",E.class)
    return(Spec)

  }else if (diversity == 'PD'){


    long = lapply(data, function(da){
      da%>%as.data.frame()%>%gather(key = "col_sp", value = "abundance")%>%column_to_rownames('col_sp')
    })
    names(long) = names(data)

    Spec <- lapply(E.class, function(e){
      each_class = lapply(seq_along(data), function(i){
        res = Spec.PD(data[[i]], q = q,datatype = datatype,
                      method = method, nboot=nboot, E.class = e, C = C)

        res['Coverage'] = NULL
        res = lapply(res, function(each_class){
          each_class%>%
            mutate(Evenness = 1-Evenness, Even.LCL = 1-Even.LCL, Even.UCL = 1-Even.UCL)%>%
            rename('Specialization'='Evenness', 'Spec.LCL' ='Even.LCL', 'Spec.UCL' ='Even.UCL')%>%
            mutate(Network = names(long)[[i]])
        })
        # if(method == "Observed") index = 1
        # if(method == "Estimated") index = 2
        # return(res[[index]]%>%mutate(Assemblage = names(long)[[i]]))
        return(res[[1]])
      })%>%do.call("rbind",.)

      each_class%>%mutate(class = paste0("1-E",e))
    })



    spec_PD <- function(){
      if (datatype == "abundance") {
        qD <- Evenness.profile(data, q, "abundance", method,
                               E.class, C)
        qD <- map(qD, as.vector)
        if (nboot > 1) {
          Prob.hat <- lapply(1:length(data), function(i) iNEXT.3D:::EstiBootComm.Ind(data[[i]]))
          Abun.Mat <- lapply(1:length(data), function(i) rmultinom(nboot,
                                                                   sum(data[[i]]), Prob.hat[[i]]))
          error = apply(matrix(sapply(1:nboot, function(b) {
            dat = lapply(1:length(Abun.Mat), function(j) Abun.Mat[[j]][,
                                                                       b])
            names(dat) = paste("Site", 1:length(dat), sep = "")
            dat.qD = Evenness.profile(dat, q, "abundance",
                                      method, E.class, C)
            unlist(dat.qD)
          }), nrow = length(q) * length(E.class) * length(Abun.Mat)),
          1, sd, na.rm = TRUE)
          error = matrix(error, ncol = length(E.class))
          se = split(error, col(error))
        }
        else {
          se = lapply(1:length(E.class), function(x) NA)
        }
        out <- lapply(1:length(E.class), function(k) {
          tmp = data.frame(Order.q = rep(q, length(data)),
                           Evenness = as.vector(qD[[k]]), s.e. = as.vector(se[[k]]),
                           Even.LCL = as.vector(qD[[k]] - qnorm(1 - (1 -
                                                                       conf)/2) * se[[k]]), Even.UCL = as.vector(qD[[k]] +
                                                                                                                   qnorm(1 - (1 - conf)/2) * se[[k]]), Assemblage = rep(names(data),                                                                                                                                                                each = length(q)), Method = rep(method, length(q) *
                                                                                                                                                                                                                                                                                                                                                                                      length(data)))
          tmp$Even.LCL[tmp$Even.LCL < 0] <- 0
          tmp
        })
        if (is.null(C) == TRUE)
          C = unique(estimate3D(data, diversity = "TD", q = 0,
                                datatype = "abundance", base = "coverage", nboot = 0)$SC)
        if (method == "Estimated") {
          out <- append(C, out)
        }
      }
    }

    Spec <- lapply(E.class, function(e){
      each_class = lapply(seq_along(long), function(i){
        res = iNEXT.4steps::Evenness(long[[i]], q = q,datatype = datatype,
                                     method = method, nboot=nboot, E.class = e, C = C)
        res['Coverage'] = NULL
        res = lapply(res, function(each_class){
          each_class%>%
            mutate(Evenness = 1-Evenness, Even.LCL = 1-Even.LCL, Even.UCL = 1-Even.UCL)%>%
            rename('Specialization'='Evenness', 'Spec.LCL' ='Even.LCL', 'Spec.UCL' ='Even.UCL')%>%
            mutate(Assemblage = names(long)[[i]])
        })
        # if(method == "Observed") index = 1
        # if(method == "Estimated") index = 2
        # return(res[[index]]%>%mutate(Assemblage = names(long)[[i]]))
        return(res[[1]])
      })%>%do.call("rbind",.)

      each_class%>%mutate(class = paste0("1 - E",e))
    })
    # ### not finished yet
    # long = lapply(data, function(da){da%>%as.data.frame()%>%gather(key = "col_sp", value = "abundance")%>%.[,2]})
    #
    # Spec <- lapply(E.class, function(e){
    #   each_class = lapply(seq_along(long), function(i){
    #     res = iNEXT.4steps::Evenness(long[[i]], q = q,datatype = datatype,
    #                                 method = method, nboot=nboot, E.class = e, C = C)
    #     res['Coverage'] = NULL
    #     res = lapply(res, function(each_class){
    #       each_class%>%
    #         mutate(Evenness = 1-Evenness, Even.LCL = 1-Even.LCL, Even.UCL = 1-Even.UCL)%>%
    #         rename('Specialization'='Evenness', 'Spec.LCL' ='Even.LCL', 'Spec.UCL' ='Even.UCL')%>%
    #         mutate(Assemblage = names(long)[[i]])
    #     })
    #     # if(method == "Observed") index = 1
    #     # if(method == "Estimated") index = 2
    #     # return(res[[index]]%>%mutate(Assemblage = names(long)[[i]]))
    #     return(res[[1]])
    #   })%>%do.call("rbind",.)
    #
    #   each_class%>%mutate(class = paste0("E",e))
    # })
    # names(Spec) = paste0("E",E.class)
    # return(Spec)
  }
}

# ggSpec.link -------------------------------------------------------------------
#' ggplot for Specialization ggSpec.link The figure for estimation of Specialization with order q
#' @param output a table generated from Specialization function
#' @return a figure of estimated sample completeness with order q
#'
#' @examples
#' data(beetles)
#' output = Spec.link(beetles)
#' ggSpec.link(output)
#' @export

ggSpec.link = function (output)
{

  if (names(output[1]) == "Coverage")
    output = output[-1]
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  classdata = cbind(do.call(rbind, output))
  fig = ggplot(classdata, aes(x = Order.q, y = Specialization, colour = Network)) +
    geom_line(size = 1.2) + geom_ribbon(aes(ymin = Spec.LCL,
                                            ymax = Spec.UCL, fill = Network), alpha = 0.2, linetype = 0) +
    scale_colour_manual(values = cbPalette) + scale_fill_manual(values = cbPalette) +
    labs(x = "Order q", y = "Specialization") + theme_bw() + theme(legend.position = "bottom",
                                                                   legend.box = "vertical", legend.key.width = unit(1.2,
                                                                                                                    "cm"), legend.title = element_blank(), legend.margin = margin(0,
                                                                                                                                                                                  0, 0, 0), legend.box.margin = margin(-10, -10, -5,
                                                                                                                                                                                                                       -10), text = element_text(size = 16), plot.margin = unit(c(5.5,
                                                                                                                                                                                                                                                                                  5.5, 5.5, 5.5), "pt")) + guides(linetype = guide_legend(keywidth = 2.5))
  if (length(output) != 1)
    fig = fig + facet_wrap(~class) + theme(strip.text.x = element_text(size = 12,
                                                                       colour = "purple", face = "bold"))
  return(fig)
}


