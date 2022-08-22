ready4beta <- function(x){
  ## transform 2d matrix to vector
  ## expand each assemblage to the union of all networks
  ## replace na to zero
  data_long <- lapply(x, function(tab){
    tab = rownames_to_column(tab, "row.name")
    long = gather(data = tab,key = "col.name", value= "abundance", -row.name)
    # long = mutate(long,int_name = paste0(col.name, "x", row.name))
    long = unite(long, "int_name",row.name:col.name, sep = "*")
    long$abundance = as.numeric(long$abundance)
    return(long)
  })

  names_tab = lapply(seq_along(x), function(i){
    sets::as.set(data_long[[i]]$int_name)
  })

  res_set = sets::set()
  for(i in 1:length(names_tab)){
    res_set = sets::set_union(names_tab[[i]], res_set)
  }

  combined = data.frame(sp = sapply(res_set, as.character))
  for(i in 1:length(names_tab)){
    combined = combined%>%
      left_join( data_long[[i]], by = c("sp" = "int_name"))
  }

  colnames(combined)[2:(1+length(x))] = paste0("abundance", 1:length(x))
  combined[is.na(combined)] = 0
  combined = column_to_rownames(combined, "sp")
}

long_to_wide = function(data_long = data_gamma){

  temp = data_long%>%as.data.frame()%>%rownames_to_column("sp")%>%
    tidyr::separate("sp", into = c("row_sp", "col_sp"), sep = "\\*")%>%
    rename("abundance"=".")
  mat = temp%>%tidyr::spread(key = "col_sp", value = "abundance")%>%column_to_rownames("row_sp")
  mat[is.na(mat)] = 0
  return(mat)
}

create.aili <- function(data,row.tree = NULL,col.tree = NULL) {

  if (class(data)[1] != "matrix") { data <- as.matrix(data) }

  if ((is.null(row.tree)) == 0 & (is.null(col.tree) == 1)){
    tip <- row.tree$tip.label[-match(rownames(data),row.tree$tip.label)]
    mytree <- drop.tip(row.tree,tip)
    mytree <- iNEXT.3D:::phylo2phytree(mytree)
    tmp <- apply(data, 2, function(abun){
      chaoUtility:::phyExpandData(x=abun, labels=rownames(data), phy=mytree, datatype="abundance")
    })
    tmp <- lapply(1:length(tmp), function(x){
      tmp1 <- as.data.frame(tmp[[x]])
      tmp1$spe.c <- colnames(data)[x]
      tmp1$interaction <- paste(tmp1$label,tmp1$spe.c,sep = "*")
      tmp1
    })
    tmp <- do.call(rbind,tmp)
    out <- data.frame(branch.abun = tmp$branch.abun, branch.length = tmp$branch.length,
                      tgroup = tmp$tgroup,interaction = tmp$interaction)
  }

  if ((is.null(row.tree) == 1) & (is.null(col.tree) == 0)){
    tip <- col.tree$tip.label[-match(colnames(data),col.tree$tip.label)]
    mytree <- drop.tip(col.tree,tip)
    mytree <- iNEXT.3D:::phylo2phytree(mytree)

    tmp <- apply(data, 1, function(abun){
      # phyBranchAL_Abu(phylo = mytree, data = abun, rootExtend = T, refT = NULL)
      chaoUtility:::phyExpandData(x=abun, labels=colnames(data), phy=mytree, datatype="abundance")
    })
    tmp <- lapply(1:length(tmp), function(x){
      tmp1 <- tmp[[x]]%>%as.data.frame()
      tmp1$spe.r <- rownames(data)[x]
      tmp1$interaction <- paste(tmp1$spe.r,tmp1$label,sep = "*")
      tmp1
    })
    tmp <- do.call(rbind,tmp)
    out <- data.frame(branch.abun = tmp$branch.abun, branch.length = tmp$branch.length,tgroup = tmp$tgroup,interaction = tmp$interaction)
    out$branch.abun[is.na(out$branch.abun)] = 0
  }

  if ((is.null(row.tree) == 0) & (is.null(col.tree) == 0)){
    col.tip <- col.tree$tip.label[-match(colnames(data),col.tree$tip.label)]
    mytree.col <- drop.tip(col.tree,col.tip)
    mytree.col <- iNEXT.3D:::phylo2phytree(mytree.col)
    row.tip <- row.tree$tip.label[-match(rownames(data),row.tree$tip.label)]
    mytree.row <- drop.tip(row.tree,row.tip)
    mytree.row <- iNEXT.3D:::phylo2phytree(mytree.row)

    # create aiLi tables by col.tree (row by row)
    tmp0 <- apply(data, 1, function(abun){
      chaoUtility:::phyExpandData(x=abun, labels=colnames(data), phy=mytree.col, datatype="abundance")
    })

    # combine
    tmp = lapply(1:length(tmp0), function(i){
      tab = tmp0[[i]]
      tab$Species <- names(tmp0)[i]
      return(tab)
    })%>%do.call("rbind",.)

    # tmp = lapply(1:nrow(data), function(i){
    #   phyExpandData(x=data[i,], labels=colnames(data), phy=mytree.col, datatype="abundance")%>%
    #     mutate(Speices = rownames(data)[i])
    # })%>%do.call("rbind",.)

    # speices = row species
    t1 <- tmp%>%as.data.frame()%>%dplyr::select(label, branch.abun, Species)
    # back to 2d matrix
    t1 <- tidyr::spread(t1,Species,branch.abun)%>%column_to_rownames("label")%>%as.matrix()
    mat <- as.matrix(t1)
    label <- unique(tmp$label)
    species <- unique(tmp$Species)

    # label = colnames
    tmp_list_bylabel = lapply(label, function(lab){
      tmp[tmp$label == lab, ]%>%head(1)
    })
    names(tmp_list_bylabel) = label

    ### row by row
    tmp1 <- apply(mat, 1, function(abun){
      phyExpandData(x=abun, labels=rownames(data), phy=mytree.row, datatype="abundance")
    })
    tmp1 = lapply(1:length(tmp1), function(i){
      tab = tmp1[[i]]
      tab$Species <- names(tmp1)[i]
      return(tab)
    })%>%do.call("rbind",.)
    tmp1[is.na(tmp1)] <- 0

    t2 <- as.data.frame(tmp1[, c("branch.length", 'label','tgroup','node.age','branch.abun','Species')])
    t2$r.length <- 0
    t2$r.group <- 0
    t2$group <- 0

    t2_list_bylabel = lapply( seq_len(length(label)), function(i){
      tab = t2[t2$Species == label[i],]
      tab$r.length = tmp_list_bylabel[[label[i]]]$'branch.length'[1]
      tab$r.group =tmp_list_bylabel[[label[i]]]$'tgroup'[1]
      return(tab)
    })
    t2 = t2_list_bylabel%>%do.call('rbind',.)

    t2$group <- ifelse(t2$tgroup == t2$r.group,t2$tgroup,"Inode")
    t2$interaction <- paste(t2$label,t2$Species,sep = "*")

    out <- data.frame(branch.abun = t2$branch.abun, branch.length = t2$branch.length*t2$r.length,
                      tgroup = t2$group, interaction = t2$interaction)
  }
  # out <- out[out$branch.abun > 0,]
  # out <- na.omit(out)
  # out <- out[order(out$tgroup,decreasing = T),]
  rownames(out) <- NULL
  return(out)
}

coverage_to_size <- function (x, C, datatype = "abundance")
{
  if (datatype == "abundance") {
    n <- sum(x)
    refC <- iNEXT.3D:::Coverage(data = x, datatype = 'abundance', n)
    # f <- function(m, C) abs(iNEXT.3D:::Chat.Ind(x, m) - C)
    f <- function(m, C) abs(iNEXT.3D:::Coverage(data = x, datatype = 'abundance', m) - C)
    if (refC == C) {
      mm = n
    }
    else if (refC > C) {
      opt <- optimize(f, C = C, lower = 0, upper = sum(x))
      mm <- opt$minimum
      # mm <- round(mm)
    }
    else if (refC < C) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }
      if (f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) +
                                   2)
      }
      if (f1 == 1 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 > 0) {
        A <- 1
      }
      mm <- (log(n/f1) + log(1 - C))/log(A) - 1
      if (is.nan(mm) == TRUE)
        mm = Inf
      mm <- n + mm
      mm <- round(mm)
    }
  }
  else {
    m <- NULL
    n <- max(x)
    refC <- iNEXT.3D:::Coverage(data = x, datatype = 'incidence_raw', n)
    # f <- function(m, C) abs(iNEXT.3D:::Chat.Sam(x, m) - C)
    f <- function(m, C) abs(iNEXT.3D:::Coverage(data = x, datatype = 'incidence_raw', m)  - C)
    if (refC == C) {
      mm = n
    }
    else if (refC > C) {
      opt <- optimize(f, C = C, lower = 0, upper = max(x))
      mm <- opt$minimum
      mm <- round(mm)
    }
    else if (refC < C) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      U <- sum(x) - max(x)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }
      if (f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) +
                                   2)
      }
      if (f1 == 1 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 > 0) {
        A <- 1
      }
      mm <- (log(U/f1) + log(1 - C))/log(A) - 1
      if (is.nan(mm) == TRUE)
        mm = Inf
      mm <- n + mm
      mm <- round(mm)
    }
  }
  return(mm)
}


Evenness.profile <- function(x, q, datatype = c("abundance","incidence_freq"), method, E.class, C = NULL) {

  if (method == "Estimated") {
    estqD = estimate3D(x, class = 'TD', q, datatype, base = "coverage", level = C, nboot = 0)
    estS = estimate3D(x, class = 'TD', 0, datatype, base = "coverage", level = C, nboot = 0)

    estqD = estimate3D(x, class = 'TD', q, datatype, base = "coverage", level = C, nboot = 0)
    estS = estimate3D(x, class = 'TD', 0, datatype, base = "coverage", level = C, nboot = 0)

    out = lapply(E.class, function(i) {
      tmp = sapply(1:length(x), function(k) even.class(q, empqD[empqD$Assemblage == names(x)[k], "qD"], empS[empS$Assemblage == names(x)[k], "qD"], i, x[[k]]/sum(x[[k]])))
      if (class(tmp)[1] %in% c("numeric","integer")) {tmp = t(as.matrix(tmp, nrow = 1))}
      rownames(tmp) = q
      tmp
    })
  } else if (method == "Empirical") {

    empqD = ObsND(x, q = q, datatype = datatype, nboot = 30)
    empS = empqD%>%filter(Order.q == 0)

    out = lapply(E.class, function(i) {
      tmp = sapply(1:length(x), function(k) even.class(q,
                                                       qD =empqD[empqD$Assemblage == names(x)[k], "qD"],
                                                       S = empS[empS$Assemblage == names(x)[k], "qD"],
                                                       E.class = i))
      # x[[k]]/sum(x[[k]])))
      if (class(tmp)[1] %in% c("numeric","integer")) {tmp = t(as.matrix(tmp, nrow = 1))}
      rownames(tmp) = q
      tmp
    })
    # empqD = Obs3D(x, class = 'TD', q = q, datatype = datatype, nboot = 0)
    # empS = Obs3D(x, class = 'TD', q = 0, datatype = datatype, nboot = 0)
    #
    # out = lapply(E.class, function(i) {
    #   tmp = sapply(1:length(x), function(k) even.class(q, empqD[empqD$Assemblage == names(x)[k], "qD"], empS[empS$Assemblage == names(x)[k], "qD"], i, x[[k]]/sum(x[[k]])))
    #   if (class(tmp)[1] %in% c("numeric","integer")) {tmp = t(as.matrix(tmp, nrow = 1))}
    #   rownames(tmp) = q
    #   tmp
    # })
  }

  names(out) = paste("E", E.class, sep="")
  return(out)
}

expanddata <- function(data){
  out <- lapply(data,function(x){
    if(ncol(x) != 1){
      tmp <- as.data.frame(x)
      tmp$spe.r <- rownames(tmp)
      tmp1 <- gather(tmp,spe.c,abun,-spe.r)
      tmp2 <- data.frame(link = paste(tmp1$spe.r,tmp1$spe.c,sep = "_"),abun = tmp1$abun)
      out <- tmp2[tmp2$abun >0,]
    }
    if(ncol(x) == 1){
      tmp <- as.data.frame(x)
      # tmp$spe.r <- rownames(tmp)
      # tmp1 <- gather(tmp,spe.c,abun,-spe.r)
      tmp <- data.frame(link = rownames(x),abun = tmp$`x[[1]]`)
      # tmp2[tmp2$abun >0,]
      out <- tmp[tmp$abun >0,]
    }
    out
  })
  tmp <- do.call(rbind,out)
  link <- unique(tmp$link)
  dama <- matrix(0,length(link),length(data),dimnames = list(link,names(data)))
  for (i in 1:length(data)) {
    dama[match(out[[i]]$link,link),i] <- out[[i]]$abun
  }
  as.data.frame(dama)
}

datainffun <- function(data, row.distM = NULL,col.distM = NULL, datatype){
  if (class(data)!="dataframe") data <- as.data.frame(data)

  if(is.null(row.distM)){
    rdd = matrix(1,ncol = nrow(data),nrow = nrow(data))
    diag(rdd) = 0
    rownames(rdd) = rownames(data)
    colnames(rdd) = rownames(data)
    row.distM =  rdd}
  if(is.null(col.distM)){
    cdd = matrix(1,ncol = ncol(data),nrow = ncol(data))
    diag(cdd) = 0
    rownames(cdd) = colnames(data)
    colnames(cdd) = colnames(data)
    col.distM =  cdd}

  res <- matrix(0,10,1,dimnames=list(1:10, "value"))

  rownames(res) <- c("n", "S.obs(row)","S.obs(col)","Links.obs","Connectance", "f1", "f2", "a1'", "a2'", "threshold")

  res[1,1] <- as.integer(sum(data))
  res[2,1] <-  nrow(data[rowSums(data)>0])%>%as.integer()
  res[3,1] <-  ncol(data[colSums(data)>0])%>%as.integer()
  res[4,1] <-  sum(data>0)%>%as.integer()

  res[5,1] <-  round(sum(data>0)/ncol(data)/nrow(data),4)
  res[6,1] <-  sum(data == 1)%>%as.integer()
  res[7,1] <-  sum(data == 2)%>%as.integer()
  distM =  1-(1-row.distM[rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM))),rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM)))])*(1-col.distM[rep(1:nrow(col.distM),nrow(row.distM)),rep(1:nrow(col.distM),nrow(row.distM))])
  distM_name = paste0(rep(colnames(col.distM),3),".",rep(rownames(row.distM),rep(ncol(col.distM),nrow(row.distM))))
  colnames(distM) = distM_name
  rownames(distM) = distM_name
  tmp <- rowSums(as.matrix(as.vector(t(as.matrix(data)))))
  tmp = matrix(tmp/sum(tmp), ncol = 1)
  dmean <- sum( (tmp %*% t(tmp) ) * distM)

  aivi = iNEXT.3D:::data_transform(data,distM,tau = dmean,datatype = datatype,integer = T)
  res[8,1] <- sum(aivi$ai==1)
  res[9,1] <- sum(aivi$ai==2)
  res[10,1] <- dmean

  res = res%>%t()%>%as.data.frame()
  return(res)
}

datainfphy <- function(data, row.tree = NULL,col.tree = NULL, datatype){
  if (class(data)!="dataframe") data <- as.data.frame(data)

  rownames(data) = gsub('\\.', '_',rownames(data))
  colnames(data) = gsub('\\.', '_',colnames(data))
  # if(datatype == "abundance"){
  #   rownames(res) <- c("n", "S.obs","Links.obs","S.obs.row","S.obs.col","Connectance", "Coverage","f1","f2","f3","f4","f5","f6","f7","f8","f9","f10")
  #   rownames(tmp) <- c("n", "S.obs","Links.obs","S.obs.row","S.obs.col","Connectance", "f1*", "f2*", "g1", "g2", "observed PD", "mean_T")
  #
  # }
  # if(datatype == "incidence_freq"){
  #   rownames(tmp) <- c("n", "S.obs","Links.obs","S.obs.row","S.obs.col","Connectance", "f1*", "f2*", "g1", "g2", "observed PD", "mean_T")
  #   rownames(res) <- c("U", "S.obs","Links.obs","S.obs.row","S.obs.col","Connectance", "Coverage","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10")
  # }
  res <- matrix(0,11,1,dimnames=list(1:11, "value"))
  # rownames(res) <- c("n", "S.obs","Links.obs","S.obs.row","S.obs.col","Connectance", "f1*", "f2*", "g1", "g2", "observed PD", "mean_T")
  rownames(res) <- c("n", "S.obs(row)","S.obs(col)","Links.obs","Connectance", "f1*", "f2*", "g1", "g2", "PD.obs", "mean_T")

  res[1,1] <- as.integer(sum(data))
  res[2,1] <-  nrow(data[rowSums(data)>0])%>%as.integer()
  res[3,1] <-  ncol(data[colSums(data)>0])%>%as.integer()
  res[4,1] <-  sum(data>0)%>%as.integer()

  res[5,1] <-  round(sum(data>0)/ncol(data)/nrow(data),4)
  res[6,1] <-  sum(data == 1)%>%as.integer()
  res[7,1] <-  sum(data == 2)%>%as.integer()
  phy <- create.aili(data,row.tree = row.tree,col.tree = col.tree)
  res[8,1] <- sum(phy$branch.length[phy$branch.abun==1])
  res[9,1] <- sum(phy$branch.length[phy$branch.abun==2])
  res[10,1] <- nrow(phy)%>%as.integer()
  res[11,1] <- sum(phy$branch.length*phy$branch.abun)/sum(data)

  res = res%>%t()%>%as.data.frame()
  return(res)
}
datainf <- function(data, datatype){
  if (class(data)!="dataframe") data <- as.data.frame(data)
  res <- matrix(0,16,1,dimnames=list(1:16, "value"))

  if(datatype == "abundance"){
    rownames(res) <- c("n", "S.obs(row)","S.obs(col)","Links.obs","Connectance", "Coverage","f1","f2","f3","f4","f5","f6","f7","f8","f9","f10")
    # rownames(res) <- c("n", "S.obs","Links.obs","S.obs.row","S.obs.col","Connectance", "Coverage","f1","f2","f3","f4","f5","f6","f7","f8","f9","f10")
  }
  if(datatype == "incidence_freq"){
    rownames(res) <- c("U", "S.obs(row)","S.obs(col)","Links.obs","S.obs.col","Connectance", "Coverage","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10")
  }
  res[]

  res[1,1] <- as.integer(sum(data))
  # res[2,1] <-  sum(ncol(data),nrow(data))
  res[2,1] <-  nrow(data[rowSums(data)>0])%>%as.integer()
  res[3,1] <-  ncol(data[colSums(data)>0])%>%as.integer()
  res[4,1] <-  sum(data>0)%>%as.integer()
  res[5,1] <-  round(sum(data>0)/ncol(data)/nrow(data),4)
  res[7:16,1] <- c(sum(data==1),sum(data==2),sum(data==3),sum(data==4),sum(data==5),
                   sum(data==6),sum(data==7),sum(data==8),sum(data==9),sum(data==10))%>%as.integer()
  f1 = sum(data==1)%>%as.integer()
  f2 = sum(data==2)%>%as.integer()
  n = sum(data)%>%as.integer()
  res[6,1] <- round(1 - (f1/n)*((n-1)*f1/((n-1)*f1+2*f2)),4) #C.hat

  res = res%>%t()%>%as.data.frame()
  return(res)
}
plot.tree2 <- function(mat){
  #number of lower level must be large than or equal to the number of higher level
  t <- apply(mat,MARGIN = 1, function(x) length(unique(x)))
  if(sum((t[-1]-t[-length(t)])<0)>0) stop("number of lower level must be large than or equal to the number of higher level, please renew your structure matrix.")
  rownames(mat) <- paste0("level", 1:nrow(mat))
  colnames(mat) <- paste0("community", 1:ncol(mat))
  mat <- data.frame(t(mat), stringsAsFactors = F)
  m <- ncol(mat)
  mat$pathString <- apply(mat,1,paste,collapse="/")
  population <- as.Node(mat)
  useRtreeList <- ToListExplicit(population, unname = TRUE)

  # radialNetwork(useRtreeList, fontSize = 10, opacity = 0.9)
  diagonalNetwork(useRtreeList,fontSize = 27, opacity = 10, linkColour = "#828282", nodeStroke = "#6495ED")
}


### ke-wei
my_PhD.q.est <- function (ai, Lis, q, nt, cal)
{
  t_bars <- as.numeric(t(ai) %*% Lis/nt)
  I1 <- which(ai == 1)
  I2 <- which(ai == 2)
  f1 <- length(I1)
  f2 <- length(I2)
  if (f2 > 0) {
    A = 2 * f2/((nt - 1) * f1 + 2 * f2)
  }
  else if (f2 == 0 & f1 > 0) {
    A = 2/((nt - 1) * (f1 - 1) + 2)
  }
  else {
    A = 1
  }
  S <- length(ai)
  if (1 %in% q) {
    ai_h1_I <- ai <= (nt - 1)
    h1_pt2 <- rep(0, S)
    ai_h1 <- ai[ai_h1_I]
    h1_pt2[ai_h1_I] <- tibble(ai = ai) %>% .[ai_h1_I, ] %>%
      mutate(diga = digamma(nt) - digamma(ai)) %>% apply(.,
                                                         1, prod)/nt
  }
  if (2 %in% q) {
    q2_pt2 <- unlist(ai * (ai - 1)/nt/(nt - 1))
  }
  if (sum(abs(q - round(q)) != 0) > 0 | max(q) > 2) {
    deltas_pt2 <- sapply(0:(nt - 1), function(k) {
      ai_delt_I <- ai <= (nt - k)
      deltas_pt2 <- rep(0, S)
      deltas_pt2[ai_delt_I] <- iNEXT.3D:::delta_part2(ai = ai[ai_delt_I],
                                                      k = k, n = nt)
      deltas_pt2
    }) %>% t()
  }
  Sub <- function(q, g1, g2, PD_obs, t_bar, Li) {
    if (q == 0) {
      ans <- PD_obs + iNEXT.3D:::PDq0(nt, f1, f2, g1, g2)
    }
    else if (q == 1) {
      h2 <- iNEXT.3D:::PDq1_2(nt, g1, A)
      h1 <- sum(Li * h1_pt2)
      h <- h1 + h2
      ans <- t_bar * exp(h/t_bar)
    }
    else if (q == 2) {
      ans <- t_bar^2/sum(Li * q2_pt2)
    }
    else {
      k <- 0:(nt - 1)
      deltas <- as.numeric(deltas_pt2 %*% Li)
      a <- (choose(q - 1, k) * (-1)^k * deltas) %>% sum
      b <- ifelse(g1 == 0 | A == 1, 0, (g1 * ((1 - A)^(1 -
                                                         nt))/nt) * (A^(q - 1) - sum(choose(q - 1, k) *
                                                                                       (A - 1)^k)))
      ans <- ((a + b)/(t_bar^q))^(1/(1 - q))
    }
    return(ans)
  }
  Lis = as.data.frame(Lis)
  est <- sapply(1:ncol(Lis), function(i) {
    Li = Lis[, i]
    t_bar <- t_bars[i]
    PD_obs <- sum(Li)
    g1 <- sum(Li[I1])
    g2 <- sum(Li[I2])
    est <- sapply(q, function(q_) Sub(q = q_, g1 = g1, g2 = g2,
                                      PD_obs = PD_obs, t_bar = t_bar, Li = Li))
  })
  if (cal == "PD") {
    est <- as.numeric(est)
  }
  else if (cal == "meanPD") {
    est <- as.numeric(sapply(1:length(t_bars), function(i) {
      est[, i]/t_bars[i]
    }))
  }
  return(est)
}


my_PhD.m.est <- function (ai, Lis, m, q, nt, cal)
{
  t_bars <- as.numeric(t(ai) %*% Lis/nt)
  if (sum(m > nt) > 0) {
    EPD = function(m, obs, asy) {
      m = m - nt
      Lis = as.data.frame(Lis)
      out <- sapply(1:ncol(Lis), function(i) {
        asy_i <- asy[, i]
        obs_i <- obs[, i]
        RPD_m_i <- RPD_m[, i]
        Li <- Lis[, i]
        t_bar <- t_bars[i]
        asy_i <- sapply(1:length(q), function(j) {
          max(asy_i[j], obs_i[j])
        })
        beta <- rep(0, length(q))
        beta0plus <- which(asy_i != obs_i)
        beta[beta0plus] <- (obs_i[beta0plus] - RPD_m_i[beta0plus])/(asy_i[beta0plus] -
                                                                      RPD_m_i[beta0plus])
        outq <- sapply(1:length(q), function(i) {
          if (q[i] != 2) {
            obs_i[i] + (asy_i[i] - obs_i[i]) * (1 -
                                                  (1 - beta[i])^m)
          }
          else if (q[i] == 2) {
            1/sum((Li/(t_bar)^2) * ((1/(nt + m)) * (ai/nt) +
                                      ((nt + m - 1)/(nt + m)) * (ai * (ai -
                                                                         1)/(nt * (nt - 1)))))
          }
        })
        outq
      })
      return(out)
    }
    RPD_m <- iNEXT.3D:::RPD(ai%>%as.matrix(), Lis%>%as.matrix(), nt, nt - 1, q)
    obs <- iNEXT.3D:::RPD(ai%>%as.matrix(), Lis%>%as.matrix(), nt, nt, q)
    asy <- matrix(my_PhD.q.est(ai = ai, Lis = Lis, q = q, nt = nt,
                               cal = cal), nrow = length(q), ncol = length(t_bars))
  }
  else if (sum(m == nt) > 0) {
    obs <- iNEXT.3D:::RPD(ai%>%as.matrix(), Lis%>%as.matrix(), nt, nt, q)
  }
  if (cal == "PD") {
    out <- sapply(m, function(mm) {
      if (mm < nt) {
        ans <- iNEXT.3D:::RPD(ai = ai%>%as.matrix(), Lis = Lis%>%as.matrix(), n = nt, m = mm,
                         q = q)
      }
      else if (mm == nt) {
        ans <- obs
      }
      else {
        ans <- EPD(m = mm, obs = obs, asy = asy)
      }
      return(as.numeric(ans))
    })
  }
  else if (cal == "meanPD") {
    out <- sapply(m, function(mm) {
      if (mm < nt) {
        ans <- iNEXT.3D:::RPD(ai = ai%>%as.matrix(), Lis = Lis%>%as.matrix(), n = nt, m = mm,
                         q = q)
      }
      else if (mm == nt) {
        ans <- obs
      }
      else {
        ans <- EPD(m = mm, obs = obs, asy = asy)
      }
      ans <- sapply(1:length(t_bars), function(i) {
        ans[, i]/t_bars[i]
      })
      as.numeric(ans)
    })
  }
  out <- matrix(out, ncol = length(m))
  return(out)
}
sample.boot <- function(exdata,B) {
  gamma <- rowSums(exdata)
  n <- sum(sapply(unique(gamma), function(x){x*sum(gamma == x)}))
  f1 <- sum(gamma == 1)
  f2 <- sum(gamma == 2)
  f0.hat <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))
  adjust.pi <- function(data) {
    f1 <- sum(data == 1)
    f2 <- sum(data == 2)
    n <- sum(sapply(unique(data), function(x){x*sum(data == x)}))
    C.hat <- 1-f1/n*ifelse(f2 ==0,ifelse(f1 == 0,0,2/((n-1)*(f1-1)+2)),2*f2/((n-1)*f1+2*f2))
    lambda <- ifelse(C.hat != 1,(1-C.hat)/sum(data/n*(1-data/n)^n),0)
    f0 <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))
    ###
    f0 <- ifelse( sum(data>0) + f0 > f0.hat +nrow(exdata),nrow(exdata) + f0.hat - sum(data>0) ,f0)
    ###
    p.seen <- data/n*(1-lambda*(1-data/n)^n)
    p.seen <- p.seen[p.seen>0]
    p.unseen <- rep((1-C.hat)/f0,f0)
    list(seen = p.seen,unseen = p.unseen)
  }
  pi <- lapply(1:ncol(exdata), function(x){adjust.pi(exdata[,x])})
  xi <- function(data){
    out <- lapply(1:ncol(data),function(x) {
      tmp <- data[,x]
      ni <- sum(sapply(unique(tmp), function(y){y*sum(tmp == y)}))
      tmp[tmp>0] <- pi[[x]]$seen
      tmp[(length(tmp)+1):(length(tmp)+f0.hat)] <- 0
      tmp[sample(which(tmp == 0),length(pi[[x]]$unseen))] <- pi[[x]]$unseen
      rmultinom(1,ni,tmp)
    })
    as.data.frame(do.call(cbind,out))
  }
  lapply(1:B, function(x){xi(exdata)})
}
sample.boot.PD <- function(data,phydata,B,row.tree = NULL,col.tree = NULL) {
  gamma <- rowSums(expanddata(data))
  n <- sum(sapply(unique(gamma), function(x){x*sum(gamma == x)}))
  f1 <- sum(gamma == 1)
  f2 <- sum(gamma == 2)
  f0.hat <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))
  adjust.pi <- function(data) {
    f1 <- sum(data == 1)
    f2 <- sum(data == 2)
    n <- sum(sapply(unique(data), function(x){x*sum(data == x)}))
    C.hat <- 1-f1/n*ifelse(f2 ==0,ifelse(f1 == 0,0,2/((n-1)*(f1-1)+2)),2*f2/((n-1)*f1+2*f2))
    lambda <- ifelse(C.hat != 1,(1-C.hat)/sum(data/n*(1-data/n)^n),0)
    f0 <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))
    p.seen <- data/n*(1-lambda*(1-data/n)^n)
    p.seen <- p.seen[p.seen>0]
    p.unseen <- (1-C.hat)/f0
    list(seen = p.seen,unseen = p.unseen)
  }
  pi <- lapply(data, function(x){adjust.pi(x)})
  xi <- function(data,phydata,row.tree,col.tree){
    out <- lapply(1:length(data),function(x) {
      tmp <- data[[x]]
      tmp <- as.matrix(tmp)
      f1 <- sum(tmp == 1)
      f2 <- sum(tmp == 2)
      f0 <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))
      phy <- phydata[[x]]
      g1 <- sum(phy$branch.length[phy$branch.abun==1])
      g2 <- sum(phy$branch.length[phy$branch.abun==2])
      g0 <- ifelse(g2>g1*f2/2/f1,(n-1)/n*g1^2/2/g2,(n-1)/n*g1*(f1-1)/2/(f2+1))/f0
      ni <- sum(sapply(unique(c(tmp)), function(y){y*sum(tmp == y)}))
      tmp[tmp>0] <- pi[[x]]$seen
      #tmp[(length(tmp)+1):(length(tmp)+f0.hat)] <- 0
      la <- sample(c(which(tmp == 0),(length(tmp)+1):(length(tmp)+f0.hat)),f0)
      tmp[la[la <= length(tmp)]] <- pi[[x]]$unseen
      p <- c(tmp,rep(pi[[x]]$unseen,sum(la> length(tmp))))
      ai <- rmultinom(1,ni,p)
      tmp <- matrix(ai[1:length(tmp)],nrow(tmp),ncol(tmp),dimnames = list(rownames(tmp),colnames(tmp)))
      out <- create.aili(tmp,row.tree = row.tree,col.tree = col.tree)
      out <- rbind(out,data.frame(branch.abun = ifelse(length(ai)>length(tmp),ai[-c(1:length(tmp))],0),branch.length=g0, tgroup="Tip",interaction = paste("Unseen species",la[la >length(tmp)]-length(tmp))))
      out
    })
    names(out) <- names(data)
    out
  }
  lapply(1:B, function(x){xi(data,phydata,row.tree,col.tree)})
}

sample.boot.phy <- function(data,B,row.tree = NULL,col.tree = NULL) {
  data <- as.matrix(data)
  # n <- sum(sapply(unique(c(data)), function(x){x*sum(data == x)}))
  n <- sum(data)
  f1 <- sum(data == 1)
  f2 <- sum(data == 2)
  f0 <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))

  adjust.pi <- function(data) {
    data_straight <- c(data)
    # f1 <- sum(data_straight == 1)
    # f2 <- sum(data_straight == 2)
    ## n <- sum(sapply(unique(data_straight), function(x){x*sum(data_straight == x)}))
    # n <- sum(data_straight)
    C.hat <- 1-f1/n*ifelse(f2 ==0,ifelse(f1 == 0,0,2/((n-1)*(f1-1)+2)),2*f2/((n-1)*f1+2*f2))
    lambda <- (1-C.hat)/sum(data_straight/n*(1-data_straight/n)^n)
    # f0 <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))
    p.seen <- data_straight/n*(1-lambda*(1-data_straight/n)^n)
    p.seen <- p.seen[p.seen>0]
    p.unseen <- (1-C.hat)/f0
    list(seen = p.seen,unseen = p.unseen)
  }
  phy <- create.aili(data,row.tree,col.tree)
  g1 <- sum(phy$branch.length[phy$branch.abun==1])
  g2 <- sum(phy$branch.length[phy$branch.abun==2])
  g0 <- ifelse(g2>g1*f2/2/f1,(n-1)/n*g1^2/2/g2,(n-1)/n*g1*(f1-1)/2/(f2+1))/f0
  pi <- adjust.pi(data)
  # a_i -> p_i
  pi_matrix = data
  pi_matrix[pi_matrix>0] <- pi$seen

  # transfrom pi from S1xS1 to B1xB2
  seen_interaction_aili = create.aili(pi_matrix,row.tree,col.tree)
  unseen_interaction_aili = data.frame(branch.abun = rep(pi$unseen,f0), branch.length = g0 / f0,
                                       tgroup = "Tip",interaction = "unseen")
  p <- rbind(seen_interaction_aili,unseen_interaction_aili)%>%
    mutate(branch.abun = ifelse(tgroup == 'Root', 1, branch.abun))

  total_nodes_num = nrow(p)

  lapply(1:B, function(x){
    ai <- rbinom(total_nodes_num,n,p$branch.abun)
    # rmultinom(n = 1, size = total_nodes_num, prob =  p$branch.abun)%>%c()
    out <- cbind(ai, p[,c("branch.length", "tgroup", "interaction")])
    colnames(out)[1] = 'branch.abun'
    # out <- data.frame(branch.abun = ai,branch.length = p$branch.length,tgroup = p$tgroup,interaction = p$interaction)
    out <- out[out$branch.abun>0,]
  })
}

get.netphydiv <- function(data,q,B,row.tree = NULL,col.tree = NULL,conf, PDtype = 'PD') {
  # plan(multisession)
  if (length(data) == 1) {
    region_names = if (is.null(names(data)))
      paste0("region_", 1)else
        names(data)
  }else {
    region_names = if (is.null(names(data)))
      paste0("region_", 1:length(data))else names(data)


  }

  # phydata <- future_lapply(data,function(x){
  #   create.aili(x,row.tree = row.tree,col.tree = col.tree)%>%
  #     filter(branch.abun > 0)
  # }, future.seed=NULL)
  phydata <- lapply(data,function(x){
    create.aili(x,row.tree = row.tree,col.tree = col.tree)%>%
      filter(branch.abun > 0)
  })

  mle <- lapply(phydata, function(x){

    # PD = iNEXT.3D:::PD.qprofile(aL = x, q = q, cal = "PD", nt = sum(x[x$tgroup == "Tip","branch.abun"]))/
    #   sum(x$branch.abun*x$branch.length)*sum(x[x$tgroup == "Tip",]$branch.abun)
    PD = iNEXT.3D:::PD.Tprofile(ai = x$branch.abun, Lis = as.matrix(x[, "branch.length", drop = F]), q = q,
                                reft = sum(x$branch.abun*x$branch.length)*sum(x[x$tgroup == "Tip",]$branch.abun), cal = "PD", nt = sum(x[x$tgroup == "Tip","branch.abun"]))/
      sum(x$branch.abun*x$branch.length)*sum(x[x$tgroup == "Tip",]$branch.abun)

    if(PDtype == 'PD'){PD = PD*sum(x$branch.abun*x$branch.length)*sum(x[x$tgroup == "Tip",]$branch.abun)}
    return(PD)
  })
  est <- lapply(phydata,function(x){
    # (ai, Lis, q, nt, cal)
    # ai = x$branch.abun;
    # Lis = x$branch.length
    # q = c(0,1,2)
    # nt = nt = sum(x[x$tgroup == "Tip",]$branch.abun)
    # cal = 'PD'

    PD = my_PhD.q.est(ai = x$branch.abun,Lis = x$branch.length,q,nt = sum(x[x$tgroup == "Tip","branch.abun"]), cal = 'PD')/
      sum(x$branch.abun*x$branch.length)*sum(x[x$tgroup == "Tip",]$branch.abun)
    if(PDtype == 'PD'){PD = PD/sum(x$branch.abun*x$branch.length)*sum(x[x$tgroup == "Tip",]$branch.abun)}
    return(PD)
  })
  boot.sam <- lapply(data, function(x){
    sample.boot.phy(x,B =B,row.tree = row.tree,col.tree = col.tree)
  })
  plan(multisession)
  mle.boot <- future_lapply(boot.sam, function(x){
    lapply(x, function(y){
      # PD = iNEXT.3D:::PD.qprofile(y,q,cal = "PD",nt = sum(y[y$tgroup == "Tip",]$branch.abun))/sum(y$branch.abun*y$branch.length)*sum(y[y$tgroup == "Tip",]$branch.abun)
      PD = iNEXT.3D:::PD.Tprofile(ai = y$branch.abun, Lis = as.matrix(y[, "branch.length", drop = F]), q = q,
                                  reft = sum(y$branch.abun*y$branch.length)*sum(y[y$tgroup == "Tip",]$branch.abun), cal = "PD", nt = sum(y[y$tgroup == "Tip","branch.abun"]))/
        sum(y$branch.abun*y$branch.length)*sum(y[y$tgroup == "Tip",]$branch.abun)

      if(PDtype == 'PD'){PD = PD*sum(y$branch.abun*y$branch.length)*sum(y[y$tgroup == "Tip",]$branch.abun)}
      return(PD)

    })
  }, future.seed=NULL)
  mle.sd <- lapply(mle.boot, function(x){
    tmp <- do.call(rbind,x)
    sapply(1:length(q), function(x){
      sd(tmp[,x])
    })
  })
  est.boot <- lapply(boot.sam, function(x){
    lapply(x, function(y){
      ## NEW
      PD = my_PhD.q.est(ai = y$branch.abun,Lis = y$branch.length,q,nt = sum(y[y$tgroup == "Tip","branch.abun"]), cal = 'PD')/
        sum(y$branch.abun*y$branch.length)*sum(y[y$tgroup == "Tip",]$branch.abun)
      if(PDtype == 'PD'){PD = PD*sum(y$branch.abun*y$branch.length)*sum(y[y$tgroup == "Tip",]$branch.abun)}
      return(PD)

    })
  })
  est.sd <- lapply(est.boot, function(x){
    tmp <- do.call(rbind,x)
    sapply(1:length(q), function(x){
      sd(tmp[,x])
    })
  })
  ci <- qnorm(conf/2+0.5)
  out <- c()
  for (i in 1:length(data)) {
    tmp <- rbind(data.frame(Order.q = q,Estimate = mle[[i]][1,],Method = "Empirical",
                            s.e. = mle.sd[[i]],
                            UCL = mle[[i]][1,] + mle.sd[[i]]*ci,
                            LCL = mle[[i]][1,] - mle.sd[[i]]*ci,
                            Region = region_names[[i]],PDtype = PDtype),
                 data.frame(Order.q = q,Estimate = est[[i]],Method = "Estimated",
                            s.e. = est.sd[[i]],
                            UCL = est[[i]] + est.sd[[i]]*ci,
                            LCL = est[[i]] - est.sd[[i]]*ci,
                            Region = region_names[[i]],
                            PDtype = PDtype))

    out <- rbind(out,tmp)
  }
  return(as.data.frame(out))
}

get.netphydiv_iNE <- function(data,q,B,row.tree = NULL,col.tree = NULL,conf, knots = 40, PDtype = 'PD') {
  if(is.null(knots)){knots = 40}
  q <- unique(ceiling(q))
  ci <- qnorm(conf/2+0.5)
  inex <- function(data,q,B,row.tree = NULL,col.tree = NULL) {
    data <- as.matrix(data)
    n <- sum(data)
    m <- sort(unique(ceiling(c(seq(1,2*n,length.out = knots),n))),decreasing = F)
    phydata <- create.aili(data,row.tree = row.tree,col.tree = col.tree)
    tbar <- sum(phydata$branch.length*phydata$branch.abun)/n
    boot.sam <- sample.boot.phy(data,B,row.tree = row.tree,col.tree = col.tree)
    sc <- Coverage(data,datatype = "abundance",m,nt =n)
    # sc <- coverage(data,m)
    sc.sd <- lapply(boot.sam,function(x){
      x <- x[x$tgroup == "Tip",]$branch.abun
      Coverage(x,datatype = "abundance",m,nt =n)
    })
    sc.sd <- do.call(cbind,sc.sd)
    sc.sd <- sapply(1:length(m), function(x){
      sd(sc.sd[x,])
    })
    sc.table <- data.frame(m=m,SC = sc, SC.UCL = sc+ci * sc.sd,SC.LCL = sc - ci * sc.sd)

    out <- lapply(q, function(q_){

      PD <- lapply(m,function(y){
        my_PhD.m.est(ai = phydata$branch.abun, Lis = phydata$branch.length, m = y, q = q_, nt = n, cal = 'PD')/tbar
      })%>%unlist()

      if(PDtype == 'meanPD'){PD = PD/tbar}

      PD.sd <- lapply(boot.sam, function(z){
        tmp <- lapply(m,function(y){
          my_PhD.m.est(ai = z$branch.abun, Lis = z$branch.length, m = y, q = q_, nt = n, cal = 'PD')/tbar
        })
        unlist(tmp)
      })
      PD.sd <- do.call(cbind,PD.sd)
      PD.sd <- sapply(1:length(m), function(x){
        sd(PD.sd[x,])
      })
      PD.table <- data.frame(m=m,method = ifelse(m<n,"interpolated",ifelse(n == m,"observed","extrapolated")),
                             Order.q = q_,PD = PD, PD.UCL = PD+ci * PD.sd,PD.LCL = PD - ci * PD.sd)
      out <- left_join(PD.table,sc.table)
      out
    })%>%do.call("rbind",.)
  }

  # out <- future_lapply(data, FUN = function(x){
  #   inex(data = x,q,B,row.tree,col.tree)
  # })
  plan(multisession)
  out <- future_lapply(data, function(x){
    inex(data = x,q,B,row.tree,col.tree)
  }, future.seed = T)

  out1 <- c()
  for (i in 1:length(out)) {
    tmp <- out[[i]]
    tmp$Region <- names(out)[i]
    out1 <- rbind(out1,tmp)
  }
  return(out1)
}
iNEXTbeta.PDlink <- function(data, level, datatype='abundance', q = c(0, 1, 2),
                             nboot = 20, conf = 0.95,
                             row.tree = NULL,col.tree = NULL){
  max_alpha_coverage = F

  if(datatype=='abundance'){

    if(class(data)=="data.frame" | class(data)=="matrix" ) data = list(Region_1 = data)

    if(class(data)== "list"){
      if(is.null(names(data))) region_names = paste0("Region_", 1:length(data)) else region_names = names(data)
      Ns = sapply(data, ncol)
      data_list = data
    }

  }

  if(is.null(conf)) conf = 0.95
  tmp = qnorm(1 - (1 - conf)/2)

  for_each_region = function(data, region_name, N){

    #data
    if (datatype=='abundance') {

      n = sum(data)
      data_gamma = rowSums(data)
      data_gamma = data_gamma[data_gamma>0]
      data_alpha = as.matrix(data) %>% as.vector

      ref_gamma = iNEXT.3D:::Coverage(data_gamma, n, datatype = 'abundance')
      ref_alpha = iNEXT.3D:::Coverage(data_alpha, n, datatype = 'abundance')

      ref_alpha_max = iNEXT.3D:::Coverage(data_gamma, n*2, datatype = 'abundance')
      ref_gamma_max = iNEXT.3D:::Coverage(data_alpha, n*2, datatype = 'abundance')

      level = level[level<1]
      level = c(level, ref_gamma, ref_alpha, ref_alpha_max, ref_gamma_max) %>% sort %>% unique

      m_gamma = sapply(level, function(i) coverage_to_size(x = data_gamma, C = i, datatype='abundance'))
      m_alpha = sapply(level, function(i) coverage_to_size(data_alpha, i, datatype='abundance'))
    }

    if (datatype=='abundance') {
      ### 1. aL_table_gamma
      data_gamma_2d = long_to_wide(data_gamma)

      aL_table_gamma = create.aili(data_gamma_2d, row.tree = row.tree, col.tree = col.tree) %>%
        select(branch.abun, branch.length, tgroup)%>%
        filter(branch.abun>0)

      ### 2. aL_table_alpha
      plan(multisession)
      aL_table_alpha =  future_lapply(1:N, function(i){

        x = as.data.frame(data[data[,i]>0,i])
        rownames(x) = rownames(data)[data[,i]>0]
        names(x) = "."

        aL_table = create.aili(data = x %>%long_to_wide(),row.tree = row.tree, col.tree=col.tree )%>%
          select(branch.abun, branch.length, tgroup)%>%
          filter(branch.abun>0)
        return(aL_table)
      }, future.seed = TRUE)%>%do.call("rbind",.)


      get_phylogenetic_alpha_gamma <- function(aL_table_gamma, aL_table_alpha, n, m_gamma, m_alpha){
        Sub = function(aL_table = aL_table_gamma, m = m_gamma){
          res = sapply(m, function(mm){
            if(mm < n){
              if(mm == round(mm)){
                qPDm <- iNEXT.3D:::PhD.m.est(ai = aL_table$branch.abun,
                                             Lis = aL_table$branch.length%>%as.matrix(),m = mm,
                                             q = q,nt = n,cal = 'PD')
                return(qPDm)
              }else {
                qPDm_raw <- iNEXT.3D:::PhD.m.est(ai = aL_table$branch.abun,
                                                 Lis = aL_table$branch.length%>%as.matrix(),m = mm,
                                                 q = q,nt = n,cal = 'PD')
                qPDm_floor <- iNEXT.3D:::PhD.m.est(ai = aL_table$branch.abun,
                                                   Lis = aL_table$branch.length%>%as.matrix(),m = floor(mm),
                                                   q = q,nt = n,cal = 'PD')
                qPDm_ceil <- iNEXT.3D:::PhD.m.est(ai = aL_table$branch.abun,
                                                  Lis = aL_table$branch.length%>%as.matrix(),m = ceiling(mm),
                                                  q = q,nt = n,cal = 'PD')
                # y = [ (y2 -y1) / (x2 - x1) ] (x - x1) + y1
                qPDm_interpolated = ((qPDm_ceil - qPDm_floor) / (ceiling(mm) - floor(mm)) ) *(mm - floor(mm)) + qPDm_floor
                return(qPDm_interpolated)
              }
            }else{
              qPDm <- iNEXT.3D:::PhD.m.est(ai = aL_table$branch.abun,
                                           Lis = aL_table$branch.length%>%as.matrix(),m = mm,
                                           q = q,nt = n,cal = 'PD')
              return(qPDm)
            }

          })
        }
        qPDm_gamma = Sub(aL_table = aL_table_gamma, m = m_gamma)
        qPDm_alpha = Sub(aL_table = aL_table_alpha, m = m_alpha)
        ## gamma
        gamma = qPDm_gamma %>% t %>%as.data.frame %>%
          set_colnames(q) %>% gather(Order, Estimate) %>%
          mutate(level=rep(level, 3), Coverage_real=rep(iNEXT.3D:::Coverage(data_gamma, m_gamma, datatype = 'abundance'), 3),
                 Size=rep(m_gamma, 3))%>%
          mutate(Method = ifelse(level>=ref_gamma, ifelse(level==ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'))
        ## alpha
        qPDm_alpha = qPDm_alpha/N
        alpha = qPDm_alpha %>% t %>% as.data.frame %>%
          set_colnames(q) %>% gather(Order, Estimate) %>%
          mutate(level=rep(level, 3), Coverage_real=rep(iNEXT.3D:::Coverage(data_alpha, m_alpha, datatype = 'abundance'), 3), Size=rep(m_alpha, 3))%>%
          mutate(Method = ifelse(level>=ref_gamma, ifelse(level==ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'))

        res = list()
        res[['gamma']] = gamma
        res[['alpha']] = alpha
        return(res)
      }

      PD_results = get_phylogenetic_alpha_gamma(aL_table_gamma, aL_table_alpha,
                                                n = sum(data_gamma), m_gamma, m_alpha)
      gamma = PD_results$gamma
      alpha = PD_results$alpha
    }

    gamma = (gamma %>%
               mutate(Method = ifelse(level>=ref_gamma,
                                      ifelse(level==ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated')))[,c(2,1,6,3,5)] %>%
      set_colnames(c('Estimate', 'Order', 'Method', 'SC', 'Size'))

    if (max_alpha_coverage==T) {
      under_max_alpha = !((gamma$Order==0) & (gamma$level>ref_alpha_max))
    }else{
      under_max_alpha = level>0
    }
    gamma = gamma[under_max_alpha,]
    gamma$Order = as.numeric(gamma$Order)


    alpha = (alpha %>%
               mutate(Method = ifelse(level>=ref_alpha, ifelse(level==ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated')))[,c(2,1,6,3,5)] %>%
      set_colnames(c('Estimate', 'Order', 'Method', 'SC', 'Size'))

    alpha = alpha[under_max_alpha,]
    alpha$Order = as.numeric(alpha$Order)

    beta = alpha
    beta$Estimate = gamma$Estimate/alpha$Estimate
    beta[beta == "Observed"] = "Observed_alpha"
    beta = beta %>% rbind(., cbind(gamma %>% filter(Method == "Observed") %>% select(Estimate) / alpha %>% filter(Method == "Observed") %>% select(Estimate),
                                   Order = q, Method = "Observed", SC = NA, Size = beta[beta$Method == "Observed_alpha", 'Size']))

    C = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(1-Order) - 1)/(N^(1-Order)-1)))
    U = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(Order-1) - 1)/(N^(Order-1)-1)))
    V = beta %>% mutate(Estimate = (Estimate-1)/(N-1))
    S = beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1))

    if(nboot>1){
      ### step0. get information from data (across samples)
      get_f0_hat = function(dat){
        dat <- as.matrix(dat)
        n <- sum(dat)
        f1 <- sum(dat == 1)
        f2 <- sum(dat == 2)
        f0 <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))
        return(f0)
      }
      get_g0_hat <- function(data){
        n = sum(data)
        f1 = sum(data==1)
        f2 = sum(data==2)

        aL = create.aili(data%>%long_to_wide(), row.tree, col.tree)%>%
          select(branch.abun,branch.length)
        g1 = aL$branch.length[aL$branch.abun==1] %>% sum
        g2 = aL$branch.length[aL$branch.abun==2] %>% sum
        g0_hat = ifelse( g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) ,
                         ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))))
        return(g0_hat)
      }

      f0_pool = get_f0_hat(data_gamma)

      plan(multisession)
      datainfo_each_k = future_lapply(1:ncol(data), function(k){
        kth_net = data[,k]
        names(kth_net) = rownames(data)

        n = sum(kth_net)
        g0_k = get_g0_hat(kth_net)
        f1 <- sum(kth_net == 1)
        f2 <- sum(kth_net == 2)
        f0_k <- ceiling(ifelse(f2 > 0 ,(n-1)/n*f1^2/2/f2,(n-1)/n*f1*(f1-1)/2))
        C.hat <- 1-f1/n*ifelse(f2 ==0,ifelse(f1 == 0,0,2/((n-1)*(f1-1)+2)),2*f2/((n-1)*f1+2*f2))
        lambda <- ifelse(C.hat != 1,(1-C.hat)/sum(kth_net/n*(1-kth_net/n)^n),0)

        p.seen <- kth_net/n*(1-lambda*(1-kth_net/n)^n)
        pi_matrix = p.seen%>%long_to_wide()
        res = list()
        res[['pi_matrix']] = pi_matrix
        res[['f0_k']] = f0_k
        res[['C.hat']] = C.hat
        res[['g0_k']] = g0_k

        return(res)
      }, future.seed = T)
      # plan(multisession)

      # se = future_lapply(1:nboot, function(b){
      se = future_lapply(1:nboot, function(b){
        ## for each bootstrap samples
        ### 1. get population (random assign position from candidates)
        aiLi_k = lapply(1:ncol(data), function(k){
          pi_matrix_k = datainfo_each_k[[k]]$pi_matrix
          f0_k = datainfo_each_k[[k]]$f0_k
          C.hat = datainfo_each_k[[k]]$C.hat
          ## pi_matrix: S1*S2 (42*98)
          ## S1*S2 -> B1*B2 (8736)
          seen_interaction_aili = create.aili(pi_matrix_k,row.tree,col.tree)
          seen_interaction_aili%>%filter(tgroup == 'Tip')%>%summarise(sum(branch.abun))
          # unseen_interaction_aili = data.frame(branch.abun = rep(pi$unseen,f0_k), branch.length = g0 / f0_k,
          #                                      tgroup = "Tip",interaction = "unseen")
          unseen_interaction_aili = data.frame(branch.abun = rep(0,f0_pool), branch.length = rep(0,f0_pool),
                                               tgroup = "Tip",interaction = "unseen")%>%
            mutate(interaction = paste0(interaction, row_number()))
          seen_unseen = rbind(seen_interaction_aili, unseen_interaction_aili)
          ## 8268 candidates out of 8977
          candidates = which(seen_unseen$branch.abun == 0)
          ###
          assigned_position = sample(candidates, size = min(f0_k, length(candidates)) , replace = F)
          seen_unseen$branch.abun[assigned_position] = (1-C.hat) / f0_k
          # seen_unseen <- seen_unseen%>%filter(branch.abun != 0)
          seen_unseen_unique = seen_unseen%>%
            group_by(interaction)%>%
            summarise(tgroup = first(tgroup),
                      branch.abun = first(branch.abun),
                      branch.length = first(branch.length))
          return(seen_unseen_unique)
        })

        length_bt = lapply(aiLi_k, function(aili){
          aili[, c("branch.length", "interaction")]%>%column_to_rownames("interaction")
        })%>%do.call("cbind",.)
        colnames(length_bt) = paste0("net_", 1:ncol(length_bt))

        p_bt = lapply(aiLi_k, function(aili){
          aili[, c("branch.abun", "interaction")]%>%column_to_rownames("interaction")
        })%>%do.call("cbind",.)
        colnames(p_bt) = paste0("net_", 1:ncol(p_bt))

        p_bt2 = p_bt%>%
          cbind(tgroup = aiLi_k[[1]]$tgroup)%>%
          filter(tgroup == 'Tip')%>%select(-tgroup)%>%
          filter_all(any_vars(. != 0))

        ### 2. generate samples(x_bt)
        x_bt = sapply(1:ncol(p_bt2), function(k){
          n = sum(data[,k])
          sapply(p_bt2[,k], function(p){rbinom(1,size=n, prob = p)})
        })
        p_star_bt = p_bt%>%
          cbind(tgroup = aiLi_k[[1]]$tgroup)%>%
          filter_all(any_vars(. != 0))

        x_star_bt = sapply(1:(ncol(p_star_bt)-1), function(k){
          n = sum(data[,k])
          sapply(p_star_bt[,k], function(p){rbinom(1,size=n, prob = p)})
        })


        rownames(x_bt) = rownames(p_bt2)
        ### 3. estimate unkwown branch length
        L0_hat = sapply(1:ncol(data), function(k){
          # compare = data.frame(interaction = names(data_gamma), raw = data_gamma)%>%
          compare = data.frame(interaction = rownames(data), raw = data[,k])%>%
            inner_join(x_bt[,k]%>%as.data.frame()%>%rownames_to_column('interaction'),
                       by = 'interaction')%>%rename('sample'=".")
          # compare$sample[c(5,14)] = 3
          known_unseen_interaction = which(and((compare$raw==0) ,(compare$sample!=0)))

          # print(paste0("length of know_unseen_interaction", length(known_unseen_interaction)))
          if(length(known_unseen_interaction) == 0){ used_length = 0}else{
            used_length = compare%>%filter(raw == 0, sample != 0 )%>%
              ## note: length_bt have same value over all k
              left_join(cbind(length_bt[,c('net_1')], interaction = rownames(length_bt))%>%
                          as.data.frame()%>%rename("length"="V1"),
                        by = 'interaction')%>%pull(length)%>%as.numeric()%>%sum()
          }
          L0_hat_k = datainfo_each_k[[k]][['g0_k']] - used_length
          # L0_hat_k = g0_hat[k] - used_length
          return(L0_hat_k)
        })

        for(k in 1:ncol(data)){
          length_bt[,k] = c(length_bt[,k]%>%.[1:(length(.)-f0_pool)],
                            rep(L0_hat[k]/ f0_pool, f0_pool))

        }
        length_bt_average = length_bt%>%rowMeans()%>%as.data.frame()%>%rownames_to_column("interaction")
        ###############
        ### With obtained sample, we can start to calculate gamma, alpha, beta diversities with data_gamma and data_alpha
        ### 1. aL_table_gamma
        bootstrap_data_gamma = rowSums(x_bt)
        bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma>0]

        bt_table_data_gamma = bootstrap_data_gamma%>%
          as.data.frame()%>%rownames_to_column("interaction")

        bt_aL_table_gamma = bt_table_data_gamma%>%
          left_join(length_bt_average , by = 'interaction')%>%
          set_colnames(c("interaction", 'branch.abun', 'branch.length'))

        ### 2. aL_table_alpha
        plan(multisession)

        bt_aL_table_alpha =  future_lapply(1:N, function(k){
          bootstrap_data_alpha = cbind(interaction = rownames(x_bt)[x_bt[,k]>0],
                                       branch.abun =  x_bt[x_bt[,k]>0,k])%>%as.data.frame()%>%
            mutate(branch.abun = as.integer(branch.abun))

          aL_table = bootstrap_data_alpha%>%
            left_join(length_bt_average , by = 'interaction')%>%
            set_colnames(c("interaction", 'branch.abun', 'branch.length'))

          return(aL_table)
        }, future.seed = TRUE)%>%do.call("rbind",.)

        bootstrap_data_alpha = as.matrix(x_bt) %>% as.vector
        bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha>0]

        m_gamma_bt = sapply(level, function(i) coverage_to_size(bootstrap_data_gamma, i, datatype='abundance'))
        m_alpha_bt = sapply(level, function(i) coverage_to_size(bootstrap_data_alpha, i, datatype='abundance'))

        get_phylogenetic_alpha_gamma_bootstrap <- function(aL_table_gamma, aL_table_alpha,
                                                           n, m_gamma, m_alpha){
          ## 1. gamma
          gamma <- iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun,
                                        Lis = aL_table_gamma$branch.length%>%as.matrix(),m = m_gamma,
                                        q = q,nt = n,cal = 'PD')%>%t%>% as.vector()
          ## 2. alpha
          alpha =(iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = aL_table_alpha$branch.length%>%as.matrix(),
                                       m = m_alpha,q = q,nt = n,cal = 'PD')/N)%>%t%>% as.vector()
          ## 3. beta
          beta_obs = (iNEXT.3D:::PD.Tprofile(ai=aL_table_gamma$branch.abun,
                                             Lis=as.matrix(aL_table_gamma$branch.length),
                                             q=q, nt=n, cal="PD") /
                        (iNEXT.3D:::PD.Tprofile(ai=aL_table_alpha$branch.abun,
                                                Lis=as.matrix(aL_table_alpha$branch.length),
                                                q=q, nt=n, cal="PD") / N)) %>% unlist()%>%as.vector()

          res = list()
          res[['gamma']] = gamma
          res[['alpha']] = alpha
          res[['beta_obs']] = beta_obs
          return(res)
        }

        PD_bootstrap_results = get_phylogenetic_alpha_gamma_bootstrap(aL_table_gamma = bt_aL_table_gamma, aL_table_alpha = bt_aL_table_alpha,
                                                                      n = sum(bootstrap_data_gamma), m_gamma_bt, m_alpha_bt)
        gamma = PD_bootstrap_results$gamma
        alpha = PD_bootstrap_results$alpha
        beta_obs = PD_bootstrap_results$beta_obs

        ##
        gamma = gamma[under_max_alpha]
        alpha = alpha[under_max_alpha]
        beta = gamma/alpha

        gamma = c(gamma, rep(0, length(q)))
        alpha = c(alpha, rep(0, length(q)))
        beta = c(beta, beta_obs)

        order = rep(q, each=length(level) + 1)[under_max_alpha]

        beta = data.frame(Estimate=beta, order)

        C = (beta %>% mutate(Estimate = ifelse(order==1,log(Estimate)/log(N),(Estimate^(1-order) - 1)/(N^(1-order)-1))))$Estimate
        U = (beta %>% mutate(Estimate = ifelse(order==1,log(Estimate)/log(N),(Estimate^(order-1) - 1)/(N^(order-1)-1))))$Estimate
        V = (beta %>% mutate(Estimate = (Estimate-1)/(N-1)))$Estimate
        S = (beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1)))$Estimate

        beta = beta$Estimate
        diversity = cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
        # diversity = cbind(gamma, alpha, beta, level= c(rep(level,2), NA,NA)) %>% as.matrix
        return(diversity)

        # })%>%
        # })%>%
      }, future.seed = T)%>%
        abind(along=3) %>% apply(1:2, sd)


    } else {

      se = matrix(0, ncol = 7, nrow = nrow(gamma))
      colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
      se = as.data.frame(se)

    }
    se = as.data.frame(se)

    gamma = gamma %>% mutate(s.e. = se$gamma[1:(nrow(se)-length(q))],
                             LCL = Estimate - tmp * se$gamma[1:(nrow(se)-length(q))],
                             UCL = Estimate + tmp * se$gamma[1:(nrow(se)-length(q))],
                             Region = region_name)

    alpha = alpha %>% mutate(s.e. = se$alpha[1:(nrow(se)-length(q))],
                             LCL = Estimate - tmp * se$alpha[1:(nrow(se)-length(q))],
                             UCL = Estimate + tmp * se$alpha[1:(nrow(se)-length(q))],
                             Region = region_name)

    beta = beta %>% mutate(  s.e. = se$beta,
                             LCL = Estimate - tmp * se$beta,
                             UCL = Estimate + tmp * se$beta,
                             Region = region_name)

    C = C %>% mutate(        s.e. = se$C,
                             LCL = Estimate - tmp * se$C,
                             UCL = Estimate + tmp * se$C,
                             Region = region_name)


    U = U %>% mutate(        s.e. = se$U,
                             LCL = Estimate - tmp * se$U,
                             UCL = Estimate + tmp * se$U,
                             Region = region_name)

    V = V %>% mutate(        s.e. = se$V,
                             LCL = Estimate - tmp * se$V,
                             UCL = Estimate + tmp * se$V,
                             Region = region_name)

    S = S %>% mutate(        s.e. = se$S,
                             LCL = Estimate - tmp * se$S,
                             UCL = Estimate + tmp * se$S,
                             Region = region_name)

    list(gamma = gamma, alpha = alpha, beta = beta, C = C, U = U, V = V, S = S)
  }

  output = lapply(1:length(data), function(i) for_each_region(data = data_list[[i]],
                                                              region_name = region_names[i], N = Ns[i]))
  names(output) = region_names
  return(output)
}

### syuan-yu
iNEXTbeta = function(data, diversity = 'TD', q = c(0, 1, 2), datatype = 'abundance',
                       base = 'coverage', level = NULL, nboot = 20, conf = 0.95,
                       PDtree = NULL, PDreftime = NULL, PDtype = 'meanPD',
                       FDdistM = NULL, FDtype = 'AUC', FDtau = NULL, FDcut_number = 30) {
  max_alpha_coverage = F
  if (datatype == 'abundance') {

    if( class(data) == "data.frame" | class(data) == "matrix" ) data = list(Region_1 = data)

    if(class(data) == "list"){
      if(is.null(names(data))) region_names = paste0("Region_", 1:length(data)) else region_names = names(data)
      Ns = sapply(data, ncol)
      data_list = data
    }

  }

  if (datatype == 'incidence_raw') {

    if(is.null(names(data))) region_names = paste0("Region_", 1:length(data)) else region_names = names(data)
    Ns = sapply(data, length)
    data_list = data

  }


  if (is.null(conf)) conf = 0.95
  tmp = qnorm(1 - (1 - conf)/2)

  trunc = ifelse(is.null(level), T, F)
  if ( is.null(level) & base == 'coverage' ) level = seq(0.5, 1, 0.025) else if ( base == 'size' ) {
    if ( is.null(level) ) {

      if (datatype == "abundance") {
        endpoint <- sapply(data_list, function(x) 2*sum(x))
      } else if (datatype == "incidence_raw") {
        endpoint <- sapply(data_list, function(x) 2*ncol(x[[1]]))
      }

      level <- lapply(1:length(data_list), function(i) {

        if(datatype == "abundance") {
          ni <- sum(data_list[[i]])
        }else if(datatype == "incidence_raw"){
          ni <- ncol(data_list[[i]][[1]])
        }

        mi <- floor(c(seq(1, ni-1, length.out = 20), ni, seq(ni+1, endpoint[i], length.out = 20)))
      })

    } else {

      if (class(level) == "numeric" | class(level) == "integer" | class(level) == "double") {
        level <- list(level = level)
      }

      if (length(level) != length(data_list)) level <- lapply(1:length(data_list), function(x) level[[1]])

      level <- lapply(1:length(data_list), function(i) {

        if (datatype == "abundance") {
          ni <- sum(data_list[[i]])
        } else if (datatype == "incidence_raw"){
          ni <- ncol(data_list[[i]][[1]])
        }

        if( sum(level[[i]] == ni) == 0 ) mi <- sort(c(ni, level[[i]])) else mi <- level[[i]]
        unique(mi)
      })
    }
  }

  if (diversity == 'FD' & FDtype == 'tau_value' & is.null(FDtau) == T) {
    if (datatype == 'abundance') {
      pdata <- sapply(data_list, rowSums) %>% rowSums
      order_sp <- match(names(pdata),rownames(FDdistM))
      FDdistM <- FDdistM[order_sp,order_sp]
      pdata <- matrix(pdata/sum(pdata), ncol = 1)
    } else if (datatype == 'incidence_raw') {
      pdata <- sapply(data_list, function(x) {tmp = Reduce('+', x); tmp[tmp > 1] = 1; rowSums(tmp) }) %>% rowSums
      order_sp <- match(names(pdata),rownames(FDdistM))
      FDdistM <- FDdistM[order_sp,order_sp]
      pdata <- matrix(pdata/sum(pdata), ncol = 1)
    }
    FDtau <- sum ( (pdata %*% t(pdata) ) * FDdistM) # dmean
  }

  if (diversity == 'PD') {

    if (datatype == "abundance")
      if (length(data_list) > 1) {
        pool.data = data_list[[1]] %>% data.frame %>% rownames_to_column()
        for (i in 2:length(data_list))
          pool.data = full_join(pool.data, data_list[[i]] %>% data.frame %>% rownames_to_column(), 'rowname')
        pool.data[is.na(pool.data)] = 0
        pool.data = pool.data %>% column_to_rownames() %>% rowSums
      } else pool.data = do.call(cbind, data_list) %>% rowSums

      if (datatype == 'incidence_raw') pool.data = do.call(cbind,lapply(data_list, function(x) do.call(cbind,x)) ) %>% rowSums

      pool.name = names(pool.data[pool.data>0])
      tip = PDtree$tip.label[-match(pool.name, PDtree$tip.label)]
      mytree = drop.tip(PDtree, tip)
      H_max = get.rooted.tree.height(mytree)

      if(is.null(PDreftime)) { reft = H_max
      } else if (PDreftime <= 0) { stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.", call. = FALSE)
      } else { reft = PDreftime }

  }

  for_each_region = function(data, region_name, N) {

    #data
    if (datatype == 'abundance') {

      n = sum(data)
      data_gamma = rowSums(data)
      data_gamma = data_gamma[data_gamma>0]
      data_alpha = as.matrix(data) %>% as.vector

      ref_gamma = iNEXT.3D:::Coverage(data_gamma, 'abundance', n)
      ref_alpha = iNEXT.3D:::Coverage(data_alpha, 'abundance', n)
      ref_alpha_max = iNEXT.3D:::Coverage(data_alpha, 'abundance', n*2)
      ref_gamma_max = iNEXT.3D:::Coverage(data_gamma, 'abundance', n*2)

      level = c(level, ref_gamma, ref_alpha, ref_alpha_max, ref_gamma_max) %>% sort %>% unique
      # level = level[level<1]

      m_gamma = sapply(level, function(i) coverage_to_size(data_gamma, i, datatype='abundance'))
      m_alpha = sapply(level, function(i) coverage_to_size(data_alpha, i, datatype='abundance'))

    }

    if (datatype == 'incidence_raw') {

      sampling_units = sapply(data, ncol)
      if (length(unique(sampling_units)) > 1) stop("unsupported data structure: the sampling units of all regions must be the same.")
      if (length(unique(sampling_units)) == 1) n = unique(sampling_units)

      gamma = Reduce('+', data)
      gamma[gamma>1] = 1
      data_gamma_raw = gamma
      data_gamma_freq = c(n, rowSums(gamma))

      data_alpha_freq = sapply(data, rowSums) %>% c(n, .)

      # data_gamma_freq = data_gamma_freq[data_gamma_freq>0]
      # data_alpha_freq = data_alpha_freq[data_alpha_freq>0]

      data_2D = apply(sapply(data, rowSums), 2, function(x) c(n, x)) %>% as.data.frame

      ref_gamma = iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', n)
      ref_alpha = iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', n)
      ref_alpha_max = iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', n*2)
      ref_gamma_max = iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', n*2)

      level = c(level, ref_gamma, ref_alpha, ref_alpha_max, ref_gamma_max) %>% sort %>% unique
      # level = level[level < 1]

      m_gamma = sapply(level, function(i) coverage_to_size(data_gamma_freq, i, datatype='incidence_freq'))
      m_alpha = sapply(level, function(i) coverage_to_size(data_alpha_freq, i, datatype='incidence_freq'))

    }



    if (diversity == 'TD') {

      if (datatype == 'abundance') {

        gamma = lapply(1:length(level), function(i){
          estimate3D(as.numeric(data_gamma), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
        }) %>% do.call(rbind,.)

        alpha = lapply(1:length(level), function(i){
          estimate3D(as.numeric(data_alpha), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
        }) %>% do.call(rbind,.)

      }

      if (datatype == 'incidence_raw') {

        gamma = lapply(1:length(level), function(i){
          estimate3D(as.numeric(data_gamma_freq), diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
        }) %>% do.call(rbind,.)

        alpha = lapply(1:length(level), function(i){
          estimate3D(data_alpha_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
        }) %>% do.call(rbind,.)



      }

      gamma = (cbind(level = rep(level, each=length(q)), gamma[,-c(1,2,8,9)]) %>%
                 mutate(Method = ifelse(level>=ref_gamma, ifelse(level == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'))
      )[,c(6,5,4,1,2,3)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'level', 'Coverage_real', 'Size'))

      # for (i in 0:2) gamma$Order[gamma$Order==paste0('q = ', i)] = i
      # gamma$Order = as.numeric(gamma$Order)

      if (max_alpha_coverage == T) under_max_alpha = !((gamma$Order == 0) & (gamma$level > ref_alpha_max)) else under_max_alpha = gamma$level > 0
      gamma = gamma[under_max_alpha,]



      alpha = (cbind(level = rep(level, each = length(q)), alpha[,-c(1,2,8,9)]) %>%
                 mutate(Method = ifelse(level >= ref_alpha, ifelse(level == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'))
      )[,c(6,5,4,1,2,3)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'level', 'Coverage_real', 'Size'))

      alpha$Estimate = alpha$Estimate / N

      # for (i in 0:2) alpha$Order[alpha$Order == paste0('q = ', i)] = i
      # alpha$Order = as.numeric(alpha$Order)

      alpha = alpha[under_max_alpha,]



      beta = alpha
      beta$Estimate = gamma$Estimate/alpha$Estimate
      beta[beta == "Observed"] = "Observed_alpha"
      beta = beta %>% rbind(., cbind(gamma %>% filter(Method == "Observed") %>% select(Estimate) / alpha %>% filter(Method == "Observed") %>% select(Estimate),
                                     Order = q, Method = "Observed", level = NA, Coverage_real = NA, Size = beta[beta$Method == "Observed_alpha", 'Size']))

      C = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(1-Order) - 1)/(N^(1-Order)-1)))
      U = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(Order-1) - 1)/(N^(Order-1)-1)))
      V = beta %>% mutate(Estimate = (Estimate-1)/(N-1))
      S = beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1))

      if(nboot>1){

        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
        #                     'datatype', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))

        # plan(sequential)
        # plan(multiprocess)

        # se = parSapply(cl, 1:nboot, function(i){

        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){

          if (datatype == 'abundance') {

            bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = bootstrap_population[,k]))

            bootstrap_data_gamma = rowSums(bootstrap_sample)
            bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma > 0]
            bootstrap_data_alpha = as.matrix(bootstrap_sample) %>% as.vector
            bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha > 0]

            gamma = lapply(1:length(level), function(i){
              estimate3D(as.numeric(bootstrap_data_gamma), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
            }) %>% do.call(rbind,.)

            alpha = lapply(1:length(level), function(i){
              estimate3D(as.numeric(bootstrap_data_alpha), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
            }) %>% do.call(rbind,.)

            beta_obs = (obs3D(as.numeric(bootstrap_data_gamma), diversity = 'TD', q = q, datatype = "abundance", nboot = 0) %>% select(qD) /
                          (obs3D(as.numeric(bootstrap_data_alpha), diversity = 'TD', q = q, datatype = "abundance", nboot = 0) %>% select(qD) / N)) %>% unlist()

          }

          if (datatype == 'incidence_raw') {

            bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')

            raw = lapply(1:ncol(bootstrap_population), function(j){

              lapply(1:nrow(bootstrap_population), function(i) rbinom(n = n, size = 1, prob = bootstrap_population[i,j])) %>% do.call(rbind,.)

            })

            gamma = Reduce('+', raw)
            gamma[gamma > 1] = 1
            bootstrap_data_gamma_freq = c(n, rowSums(gamma))

            bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)

            bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq > 0]
            bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq > 0]

            gamma = lapply(1:length(level), function(i){
              estimate3D(bootstrap_data_gamma_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
            }) %>% do.call(rbind,.)

            alpha = lapply(1:length(level), function(i){
              estimate3D(bootstrap_data_alpha_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
            }) %>% do.call(rbind,.)

            beta_obs = (obs3D(as.numeric(bootstrap_data_gamma_freq), diversity = 'TD', q = q, datatype = "incidence_freq", nboot = 0) %>% select(qD) /
                          (obs3D(as.numeric(bootstrap_data_alpha_freq), diversity = 'TD', q = q, datatype = "incidence_freq", nboot = 0) %>% select(qD) / N)) %>% unlist()

          }

          gamma = gamma[,c(6,3,7)]$qD[under_max_alpha]

          alpha = alpha[,c(6,3,7)]$qD[under_max_alpha]
          alpha = alpha / N

          beta = gamma/alpha

          gamma = c(gamma, rep(0, length(q)))
          alpha = c(alpha, rep(0, length(q)))
          beta = c(beta, beta_obs)

          order = rep(q, length(level) + 1)[under_max_alpha]

          beta = data.frame(Estimate=beta, order)

          C = (beta %>% mutate(Estimate = ifelse(order == 1,log(Estimate)/log(N),(Estimate^(1 - order) - 1)/(N^(1 - order) - 1))))$Estimate
          U = (beta %>% mutate(Estimate = ifelse(order == 1,log(Estimate)/log(N),(Estimate^(order - 1) - 1)/(N^(order - 1) - 1))))$Estimate
          V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
          S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate

          beta = beta$Estimate

          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix

          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }) %>% abind(along = 3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start

        # stopCluster(cl)
        # plan(sequential)

      } else {

        se = matrix(0, ncol = 7, nrow = nrow(beta))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)

      }

      if (nboot>0 & datatype == 'abundance') {
        gamma.se = estimate3D(data_gamma, diversity = 'TD', q = q, datatype = datatype, base = "coverage", level = level, nboot = nboot) %>% arrange(., goalSC, Order.q) %>% select(s.e.) %>% unlist

        alpha.se = estimate3D(data_alpha, diversity = 'TD', q = q, datatype = datatype, base = "coverage", level = level, nboot = nboot) %>% arrange(., goalSC, Order.q) %>% select(s.e.) %>% unlist
        alpha.se = alpha.se / N

      } else if (nboot>0 & datatype == 'incidence_raw') {

        gamma.se = estimate3D(data_gamma_freq, diversity = 'TD', q = q, datatype = 'incidence_freq', base = "coverage", level = level, nboot = nboot) %>% arrange(., goalSC, Order.q) %>% select(s.e.) %>% unlist

        alpha.se = estimate3D(data_alpha_freq, diversity = 'TD', q = q, datatype = 'incidence_freq', base = "coverage", level = level, nboot = nboot) %>% arrange(., goalSC, Order.q) %>% select(s.e.) %>% unlist
        alpha.se = alpha.se / N

      } else {
        gamma.se = alpha.se = 0
      }

      se[1:( length(level) * length(q) ), 'gamma'] = gamma.se

      se[1:( length(level) * length(q) ), 'alpha'] = alpha.se

    }

    if (diversity == 'PD') {

      if (datatype == 'abundance') {

        aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = data_gamma, rootExtend = T, refT = reft)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

        gamma = iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = m_gamma, q = q, nt = n, cal = "PD") %>% t %>% as.data.frame %>%
          set_colnames(q) %>% gather(Order, Estimate) %>%
          mutate(level = rep(level, length(q)), Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma, 'abundance', m_gamma), length(q)), Size = rep(m_gamma, length(q)))


        aL_table_alpha = c()

        for (i in 1:N){

          x = data[data[,i]>0,i]
          names(x) = rownames(data)[data[,i]>0]

          aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

          aL_table_alpha = rbind(aL_table_alpha, aL_table)

        }


        qPDm = iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = m_alpha, q = q, nt = n, cal = "PD")
        qPDm = qPDm/N
        alpha = qPDm %>% t %>% as.data.frame %>%
          set_colnames(q) %>% gather(Order, Estimate) %>%
          mutate(level = rep(level, length(q)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha, 'abundance', m_alpha), length(q)), Size = rep(m_alpha, length(q)))

      }

      if (datatype == 'incidence_raw') {

        aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = as.matrix(data_gamma_raw), datatype = "incidence_raw", refT = reft, rootExtend = T)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
        gamma = iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis=as.matrix(aL_table_gamma$branch.length), m = m_gamma, q = q, nt = n, cal = "PD") %>% t %>% as.data.frame %>%
          set_colnames(q) %>% gather(Order, Estimate) %>%
          mutate(level = rep(level, length(q)), Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', m_gamma), length(q)), Size = rep(m_gamma, length(q)))

        aL_table_alpha = c()

        for (i in 1:N){

          x = data[[i]]

          aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

          aL_table_alpha = rbind(aL_table_alpha, aL_table)

        }

        alpha = (iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = m_alpha, q = q, nt = n, cal = "PD")/N) %>% t %>% as.data.frame %>%
          set_colnames(q) %>% gather(Order, Estimate) %>%
          mutate(level = rep(level, length(q)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', m_alpha), length(q)), Size = rep(m_alpha, length(q)))


      }

      gamma = (gamma %>%
                 mutate(Method = ifelse(level >= ref_gamma, ifelse(level == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated')))[,c(2,1,6,3,4,5)] %>%
        set_colnames(c('Estimate', 'Order', 'Method', 'level', 'Coverage_real', 'Size'))

      if (max_alpha_coverage == T) under_max_alpha = !((gamma$Order == 0) & (gamma$level > ref_alpha_max)) else under_max_alpha = gamma$level>0
      gamma = gamma[under_max_alpha,]
      gamma$Order = as.numeric(gamma$Order)


      alpha = (alpha %>%
                 mutate(Method = ifelse(level >= ref_alpha, ifelse(level == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated')))[,c(2,1,6,3,4,5)] %>%
        set_colnames(c('Estimate', 'Order', 'Method', 'level', 'Coverage_real', 'Size'))

      alpha = alpha[under_max_alpha,]
      alpha$Order = as.numeric(alpha$Order)

      if (PDtype == 'meanPD') {
        gamma$Estimate = gamma$Estimate/reft
        alpha$Estimate = alpha$Estimate/reft
      }

      beta = alpha
      beta$Estimate = gamma$Estimate/alpha$Estimate
      beta[beta == "Observed"] = "Observed_alpha"
      beta = beta %>% rbind(., cbind(gamma %>% filter(Method == "Observed") %>% select(Estimate) / alpha %>% filter(Method == "Observed") %>% select(Estimate),
                                     Order = q, Method = "Observed", level = NA, Coverage_real = NA, Size = beta[beta$Method == "Observed_alpha", 'Size']))

      C = beta %>% mutate(Estimate = ifelse(Order == 1, log(Estimate)/log(N), (Estimate^(1 - Order) - 1)/(N^(1 - Order) - 1)))
      U = beta %>% mutate(Estimate = ifelse(Order == 1, log(Estimate)/log(N), (Estimate^(Order - 1) - 1)/(N^(Order - 1) - 1)))
      V = beta %>% mutate(Estimate = (Estimate - 1)/(N - 1))
      S = beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1))

      if(nboot>1){

        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
        #                     'datatype', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))

        # plan(sequential)
        # plan(multiprocess)

        # se = parSapply(cl, 1:nboot, function(i){

        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){

          if (datatype == 'abundance') {

            tree_bt = PDtree

            bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data)),] %>% matrix(ncol = ncol(data))

            if ( nrow(p_bt) > nrow(data) & sum(unseen_p) > 0 ){

              unseen = unseen_p[which(rowSums(unseen_p) > 0),]
              unseen = matrix(unseen, ncol = ncol(unseen_p), byrow = T)
              p_bt = rbind(p_bt[(1:nrow(data)),], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
              rownames(p_bt) = c(rownames(data), unseen_name)

              bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
              x_bt = bootstrap_sample

              rownames(x_bt) = rownames(p_bt)

              if ( sum(x_bt[-(1:nrow(data)),])>0 ){

                g0_hat = apply(data, 2, function(x){

                  n = sum(x)
                  f1 = sum(x == 1)
                  f2 = sum(x == 2)

                  aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)

                  aL$treeNabu$branch.length = aL$BLbyT[,1]
                  aL = aL$treeNabu %>% select(branch.abun,branch.length)
                  g1 = aL$branch.length[aL$branch.abun == 1] %>% sum
                  g2 = aL$branch.length[aL$branch.abun == 2] %>% sum
                  g0_hat = ifelse( g2 > ((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
                  if(is.na(g0_hat)) {g0_hat <- 0 }
                  g0_hat

                })

                te = (x_bt[1:nrow(data),]*(data == 0))>0
                used_length = sapply(1:ncol(data), function(i) {

                  if (sum(te[,i]) == 0) return(0) else {

                    iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x_bt[1:nrow(data),i], rootExtend = T, refT = reft)$treeNabu %>%
                      subset(label %in% names(which(te[,i] == TRUE))) %>% select(branch.length) %>% sum

                  }

                })

                g0_hat = g0_hat - used_length
                g0_hat[g0_hat < 0] = 0

                unseen_sample = x_bt[-(1:nrow(data)),]
                if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol = ncol(x_bt), byrow = T)

                L0_hat = sapply(1:length(g0_hat), function(i) if(sum(unseen_sample[,i] > 0) > 0) (g0_hat[i] / nrow(unseen)) else 0 )

                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow = T) * unseen_sample)) / rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample) == 0)] = 0

                for (i in 1:length(L0_hat)){

                  tip = list(edge = matrix(c(2,1),1,2),
                             tip.label = unseen_name[i],
                             edge.length = L0_hat[i],
                             Nnode = 1)
                  class(tip) = "phylo"

                  tree_bt = tree_bt + tip

                }

              } else {

                x_bt = x_bt[1:nrow(data),]
                p_bt = p_bt[1:nrow(data),]

              }

            } else {

              p_bt = p_bt[1:nrow(data),]
              x_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
              rownames(x_bt) = rownames(data)

            }

            bootstrap_data_gamma = rowSums(x_bt)
            bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma>0]
            bootstrap_data_alpha = as.matrix(x_bt) %>% as.vector
            bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha>0]

            m_gamma = sapply(level, function(i) coverage_to_size(bootstrap_data_gamma, i, datatype='abundance'))
            m_alpha = sapply(level, function(i) coverage_to_size(bootstrap_data_alpha, i, datatype='abundance'))

            aL = iNEXT.3D:::phyBranchAL_Abu(phylo = tree_bt, data = bootstrap_data_gamma, rootExtend = T, refT = reft)
            aL$treeNabu$branch.length = aL$BLbyT[,1]
            aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

            gamma = as.vector(iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = m_gamma, q = q, nt = n, cal = "PD") %>% t)


            aL_table_alpha = c()

            for (i in 1:N){

              # x = x_bt[x_bt[,i]>0,i]
              # names(x) = rownames(p_bt)[x_bt[,i]>0]

              x = x_bt[,i]
              names(x) = rownames(p_bt)
              x = x[x_bt[,i]>0]

              aL = iNEXT.3D:::phyBranchAL_Abu(phylo = tree_bt, data = x, rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[,1]
              aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

              aL_table_alpha = rbind(aL_table_alpha, aL_table)

            }

            alpha = as.vector((iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = m_alpha, q = q, nt = n, cal = "PD")/N) %>% t)

            beta_obs = (iNEXT.3D:::PD.Tprofile(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), q = q, nt = n, cal = "PD") /
                          (iNEXT.3D:::PD.Tprofile(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), q = q, nt = n, cal = "PD") / N)) %>% unlist()
          }

          if (datatype == 'incidence_raw') {

            tree_bt = PDtree

            bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data[[1]])),] %>% matrix(ncol=N)

            if ( nrow(p_bt) > nrow(data[[1]]) & sum(unseen_p)>0 ){

              unseen = unseen_p[which(rowSums(unseen_p)>0),]
              unseen = matrix(unseen, ncol = ncol(unseen_p), byrow = T)
              p_bt = rbind(p_bt[(1:nrow(data[[1]])),], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
              rownames(p_bt) = c(rownames(data[[1]]), unseen_name)

              raw = lapply(1:ncol(p_bt), function(j){

                lapply(1:nrow(p_bt), function(i) rbinom(n=n, size=1, prob=p_bt[i,j])) %>% do.call(rbind,.)

              })

              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)

              if ( lapply(1:length(raw), function(i) raw[[i]][-(1:nrow(data[[1]])),]) %>% do.call(sum,.)>0 ){

                R0_hat = sapply(data, function(x){

                  nT = ncol(x)
                  Q1 = sum(rowSums(x)==1)
                  Q2 = sum(rowSums(x)==2)

                  aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)

                  aL$treeNabu$branch.length = aL$BLbyT[,1]
                  aL = aL$treeNabu %>% select(branch.abun,branch.length)
                  R1 = aL$branch.length[aL$branch.abun == 1] %>% sum
                  R2 = aL$branch.length[aL$branch.abun == 2] %>% sum
                  R0_hat = ifelse( R2>((R1*Q2)/(2*Q1)) , ((nT-1)/nT)*(R1^2/(2*R2)) , ((nT-1)/nT)*(R1*(Q1-1)/(2*(Q2+1))) )
                  if(is.na(R0_hat)) { R0_hat <- 0 }
                  R0_hat

                })

                te = (sapply(raw, rowSums)[1:nrow(data[[1]]),]*(sapply(data, rowSums) == 0)) > 0
                used_length = sapply(1:N, function(i) {

                  if (sum(te[,i]) == 0) return(0) else {

                    iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = raw[[i]][1:nrow(data[[1]]),], datatype = "incidence_raw", rootExtend = T, refT = reft)$treeNabu %>%
                      subset(label %in% names(which(te[,i] == TRUE))) %>% select(branch.length) %>% sum

                  }

                })

                R0_hat = R0_hat - used_length
                R0_hat[R0_hat < 0] = 0

                unseen_sample = sapply(raw, rowSums)[-(1:nrow(data[[1]])),]
                if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol = N, byrow = T)

                L0_hat = sapply(1:length(R0_hat), function(i) if(sum(unseen_sample[,i] > 0) > 0) (R0_hat[i] / nrow(unseen)) else 0 )

                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow = T) * unseen_sample)) / rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample) == 0)] = 0

                for (i in 1:length(L0_hat)){

                  tip = list(edge = matrix(c(2,1), 1, 2),
                             tip.label = unseen_name[i],
                             edge.length = L0_hat[i],
                             Nnode = 1)
                  class(tip) = "phylo"

                  tree_bt = tree_bt + tip

                }

              } else raw = lapply(raw, function(i) i[1:nrow(data[[1]]),])

            } else {

              p_bt = p_bt[1:nrow(data[[1]]),]
              raw = lapply(1:ncol(p_bt), function(j){

                lapply(1:nrow(p_bt), function(i) rbinom(n = n, size = 1, prob = p_bt[i,j])) %>% do.call(rbind,.)

              })

              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)

            }

            gamma = Reduce('+', raw)
            gamma[gamma > 1] = 1
            bootstrap_data_gamma_raw = gamma
            bootstrap_data_gamma_freq = c(n, rowSums(gamma))

            bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)

            bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq > 0]
            bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq > 0]

            m_gamma = sapply(level, function(i) coverage_to_size(bootstrap_data_gamma_freq, i, datatype = 'incidence'))
            m_alpha = sapply(level, function(i) coverage_to_size(bootstrap_data_alpha_freq, i, datatype = 'incidence'))

            aL = iNEXT.3D:::phyBranchAL_Inc(phylo = tree_bt, data = bootstrap_data_gamma_raw, datatype = "incidence_raw", rootExtend = T, refT = reft)
            aL$treeNabu$branch.length = aL$BLbyT[,1]
            aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

            gamma = as.vector(iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = m_gamma, q = q, nt = n, cal="PD") %>% t)


            aL_table_alpha = c()

            for (i in 1:N){

              x = raw[[i]]

              aL = iNEXT.3D:::phyBranchAL_Inc(phylo = tree_bt, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[,1]
              aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

              aL_table_alpha = rbind(aL_table_alpha, aL_table)

            }

            alpha = as.vector((iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = m_alpha, q=q, nt = n, cal = "PD")/N) %>% t)

            beta_obs = (iNEXT.3D:::PD.Tprofile(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), q = q, nt = n, cal = "PD") /
                          (iNEXT.3D:::PD.Tprofile(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), q = q, nt = n, cal = "PD") / N)) %>% unlist()
          }

          gamma = gamma[under_max_alpha]

          alpha = alpha[under_max_alpha]

          if (PDtype == 'meanPD') {
            gamma = gamma/reft
            alpha = alpha/reft
          }

          beta = gamma/alpha

          gamma = c(gamma, rep(0, length(q)))
          alpha = c(alpha, rep(0, length(q)))
          beta = c(beta, beta_obs)

          order = rep(q, each = length(level) + 1)[under_max_alpha]

          beta = data.frame(Estimate = beta, order)

          C = (beta %>% mutate(Estimate = ifelse(order == 1, log(Estimate)/log(N), (Estimate^(1 - order) - 1)/(N^(1 - order) - 1))))$Estimate
          U = (beta %>% mutate(Estimate = ifelse(order == 1, log(Estimate)/log(N), (Estimate^(order - 1) - 1)/(N^(order - 1) - 1))))$Estimate
          V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
          S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate

          beta = beta$Estimate

          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix

          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }) %>% abind(along = 3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start

        # stopCluster(cl)
        # plan(sequential)

      } else {

        se = matrix(0, ncol = 7, nrow = nrow(beta))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)

      }

    }

    if (diversity == 'FD') {

      FD_by_tau = function(data, distM, tau, level, datatype, by, m_gamma, m_alpha) {

        if (datatype == 'abundance') {

          zik = data
          zik = zik[rowSums(data)>0,]

          dij = distM
          dij = dij[rowSums(data)>0, rowSums(data)>0]

          dij[which(dij>tau, arr.ind = T)] = tau
          aik = (1 - dij/tau) %*% as.matrix(zik)
          positive_id = rowSums(aik)>0


          gamma_x = rowSums(zik)[positive_id]
          gamma_a = rowSums(aik)[positive_id]
          gamma_v = gamma_x/gamma_a
          # gamma_a = ifelse(gamma_a < 1, 1, round(gamma_a))

          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))

          gamma = FD.m.est_0(ai_vi_gamma, m_gamma, q, n) %>% as.vector
          # gamma = iNEXT.3D:::FD.m.est(ai_vi_gamma, m_gamma, q, n) %>% as.vector


          alpha_x = as.vector(as.matrix(zik))
          alpha_a = as.vector(aik)
          # alpha_a = ifelse(alpha_a < 1, 1, round(alpha_a))

          alpha_v = alpha_x/alpha_a
          alpha_v = rep(gamma_v,N)

          alpha_v = alpha_v[alpha_a>0]
          alpha_a = alpha_a[alpha_a>0]

          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))

          if (by == 'size') alpha = (FD.m.est_0(ai_vi_alpha, m_gamma, q, n)/N) %>% as.vector
          if (by == 'coverage') alpha = (FD.m.est_0(ai_vi_alpha, m_alpha, q, n)/N) %>% as.vector
          # if (by == 'size') alpha = (iNEXT.3D:::FD.m.est(ai_vi_alpha, m_gamma, q, n)/N) %>% as.vector
          # if (by == 'coverage') alpha = (iNEXT.3D:::FD.m.est(ai_vi_alpha, m_alpha, q, n)/N) %>% as.vector

        }

        if (datatype == 'incidence_raw') {

          data_gamma_freq = data$data_gamma_freq
          data_2D = data$data_2D

          gamma_Y = data_gamma_freq[-1]

          dij = distM
          dij = dij[gamma_Y > 0, gamma_Y > 0]
          gamma_Y = gamma_Y[gamma_Y > 0]

          dij[which(dij>tau, arr.ind = T)] = tau
          gamma_a = (1 - dij/tau) %*% as.matrix(gamma_Y)
          gamma_a[gamma_a > n] = n

          gamma_v = gamma_Y/gamma_a

          # gamma_a = ifelse(gamma_a < 1, 1, round(gamma_a))
          # gamma_a[gamma_a > n] = n

          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))

          gamma = FD.m.est_0(ai_vi_gamma, m_gamma, q, n) %>% as.vector
          # gamma = iNEXT.3D:::FD.m.est(ai_vi_gamma, m_gamma, q, n) %>% as.vector


          alpha_Y = data_2D[-1,]

          dij = distM
          dij = dij[rowSums(data_2D[-1,]) > 0, rowSums(data_2D[-1,])>0]
          alpha_Y = alpha_Y[rowSums(data_2D[-1,])>0,]

          dij[which(dij>tau, arr.ind = T)] = tau
          alpha_a = (1 - dij/tau) %*% as.matrix(alpha_Y)

          # alpha_a = ifelse(alpha_a < 1, 1, round(alpha_a))
          alpha_a[alpha_a > n] = n
          alpha_a = as.vector(alpha_a)

          alpha_v = rep(gamma_v, N)
          alpha_v = alpha_v[alpha_a > 0]
          alpha_a = alpha_a[alpha_a > 0]

          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))

          if (by == 'size') alpha = (FD.m.est_0(ai_vi_alpha, m_gamma, q, n)/N) %>% as.vector
          if (by == 'coverage') alpha = (FD.m.est_0(ai_vi_alpha, m_alpha, q, n)/N) %>% as.vector
          # if (by == 'size') alpha = (iNEXT.3D:::FD.m.est(ai_vi_alpha, m_gamma, q, n)/N) %>% as.vector
          # if (by == 'coverage') alpha = (iNEXT.3D:::FD.m.est(ai_vi_alpha, m_alpha, q, n)/N) %>% as.vector

        }

        return(data.frame(gamma,alpha))

      }

      if (FDtype == 'tau_value'){

        if (datatype == 'abundance') {

          FDdistM = FDdistM[rownames(FDdistM) %in% rownames(data), colnames(FDdistM) %in% rownames(data)]
          order_sp <- match(rownames(data),rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]

          output = FD_by_tau(data, FDdistM, FDtau, level, datatype = 'abundance', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)
          gamma = output$gamma
          alpha = output$alpha

          gamma = data.frame(level, gamma) %>%
            mutate(Method = ifelse(level>=ref_gamma, ifelse(level == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each=length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma, 'abundance', m_gamma), length(q)), Size = rep(m_gamma, length(q)))

          alpha = data.frame(level, alpha) %>%
            mutate(Method = ifelse(level>=ref_alpha, ifelse(level == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each=length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha, 'abundance', m_alpha), length(q)), Size = rep(m_alpha, length(q)))

          beta = alpha
          beta$alpha = gamma$gamma/alpha$alpha

          ## Observed Beta ##
          output_obs = FD_by_tau(data, FDdistM, FDtau, level, datatype = 'abundance', by = 'size', m_gamma = sum(data), m_alpha = sum(data))
          gamma_obs = output_obs$gamma
          alpha_obs = output_obs$alpha
          obs_beta = gamma_obs/alpha_obs

        }

        if (datatype == 'incidence_raw') {

          output = FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, FDtau, level, datatype='incidence_raw', by = 'coverage', m_gamma=m_gamma, m_alpha=m_alpha)
          gamma = output$gamma
          alpha = output$alpha

          gamma = data.frame(level, gamma) %>%
            mutate(Method = ifelse(level >= ref_gamma, ifelse(level == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each = length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', m_gamma), length(q)), Size = rep(m_gamma, length(q)))

          alpha = data.frame(level, alpha) %>%
            mutate(Method = ifelse(level>=ref_alpha, ifelse(level == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each = length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', m_alpha), length(q)), Size = rep(m_alpha, length(q)))

          beta = alpha
          beta$alpha = gamma$gamma/alpha$alpha

          ## Observed Beta ##
          output_obs = FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, FDtau, level, datatype = 'incidence_raw', by = 'size', m_gamma = data_gamma_freq[1], m_alpha = data_gamma_freq[1])
          gamma_obs = output_obs$gamma
          alpha_obs = output_obs$alpha
          obs_beta = gamma_obs/alpha_obs

        }

      }

      if (FDtype == 'AUC'){

        cut = seq(0.00000001, 1, length.out = FDcut_number)
        width = diff(cut)

        if (datatype == 'abundance') {

          FDdistM = FDdistM[rownames(FDdistM) %in% rownames(data), colnames(FDdistM) %in% rownames(data)]
          order_sp <- match(rownames(data),rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]

          gamma_alpha_over_tau = lapply(cut, function(tau) {

            FD_by_tau(data, FDdistM, tau, level, datatype = 'abundance', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)

          })

          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)

          left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)

          gamma = colSums((left_limit + right_limit)/2)

          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)

          left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)

          alpha = colSums((left_limit + right_limit)/2)

          beta_over_tau = gamma_over_tau/alpha_over_tau

          left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)

          beta = colSums((left_limit + right_limit)/2)

          gamma = data.frame(level, gamma) %>%
            mutate(Method = ifelse(level >= ref_gamma, ifelse(level == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each = length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma, 'abundance', m_gamma), length(q)), Size = rep(m_gamma, length(q)))

          alpha = data.frame(level, alpha) %>%
            mutate(Method = ifelse(level >= ref_alpha, ifelse(level == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each = length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha, 'abundance', m_alpha), length(q)), Size = rep(m_alpha, length(q)))

          beta = data.frame(level, beta) %>%
            mutate(Method = ifelse(level >= ref_alpha, ifelse(level == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each = length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha, 'abundance', m_alpha), length(q)), Size = rep(m_alpha, length(q)))

          ## Observed Beta ##
          obs_gamma_alpha_over_tau = lapply(cut, function(tau) {

            FD_by_tau(data, FDdistM, tau, level, datatype = 'abundance', by = 'size', m_gamma = sum(data), m_alpha = sum(data))

          })

          obs_beta_over_tau = sapply(obs_gamma_alpha_over_tau, function(x) x$gamma) / sapply(obs_gamma_alpha_over_tau, function(x) x$alpha)

          if (length(q) == 1) obs_beta_over_tau = matrix(obs_beta_over_tau, nrow = 1)

          obs_beta = colSums( (apply(obs_beta_over_tau, 1, function(x) x[-FDcut_number]*width) + apply(obs_beta_over_tau, 1, function(x) x[-1]*width) ) / 2)

        }

        if (datatype == 'incidence_raw') {

          FDdistM = FDdistM[rownames(FDdistM) %in% names(data_gamma_freq)[-1], colnames(FDdistM) %in% names(data_gamma_freq)[-1]]
          order_sp <- match(names(data_gamma_freq)[-1],rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]

          gamma_alpha_over_tau = lapply(cut, function(tau) {

            FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, tau, level, datatype = 'incidence_raw', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)

          })

          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)

          left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)

          gamma = colSums((left_limit + right_limit)/2)

          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)

          left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)

          alpha = colSums((left_limit + right_limit)/2)

          beta_over_tau = gamma_over_tau/alpha_over_tau

          left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)

          beta = colSums((left_limit + right_limit)/2)

          gamma = data.frame(level, gamma) %>%
            mutate(Method = ifelse(level >= ref_gamma, ifelse(level == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each = length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', m_gamma), length(q)), Size = rep(m_gamma, length(q)))

          alpha = data.frame(level, alpha) %>%
            mutate(Method = ifelse(level >= ref_alpha, ifelse(level == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each = length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', m_alpha), length(q)), Size = rep(m_alpha, length(q)))

          beta = data.frame(level, beta) %>%
            mutate(Method = ifelse(level >= ref_alpha, ifelse(level == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each = length(level)/length(q)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', m_alpha), length(q)), Size = rep(m_alpha, length(q)))

          ## Observed Beta ##
          obs_gamma_alpha_over_tau = lapply(cut, function(tau) {

            FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, tau, level, datatype = 'incidence_raw', by = 'size', m_gamma = data_gamma_freq[1], m_alpha = data_alpha_freq[1])

          })

          obs_beta_over_tau = sapply(obs_gamma_alpha_over_tau, function(x) x$gamma) / sapply(obs_gamma_alpha_over_tau, function(x) x$alpha)

          if (length(q) == 1) obs_beta_over_tau = matrix(obs_beta_over_tau, nrow = 1)

          obs_beta = colSums( (apply(obs_beta_over_tau, 1, function(x) x[-FDcut_number]*width) + apply(obs_beta_over_tau, 1, function(x) x[-1]*width) ) / 2)

        }

      }

      gamma = gamma[,c(2,4,3,1,5,6)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'level', 'Coverage_real', 'Size'))

      if (max_alpha_coverage == T) under_max_alpha = !((gamma$Order == 0) & (gamma$level > ref_alpha_max)) else under_max_alpha = gamma$level > 0
      gamma = gamma[under_max_alpha,]


      alpha = alpha[,c(2,4,3,1,5,6)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'level', 'Coverage_real', 'Size'))

      alpha = alpha[under_max_alpha,]

      beta = beta[,c(2,4,3,1,5,6)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'level', 'Coverage_real', 'Size'))

      beta = beta[under_max_alpha,]

      beta[beta == "Observed"] = "Observed_alpha"
      beta = beta %>%
        rbind(., data.frame(Estimate = obs_beta, Order = q, Method = "Observed", level = NA, Coverage_real = NA, Size = beta[beta$Method == "Observed_alpha", 'Size']))

      C = beta %>% mutate(Estimate = ifelse(Order == 1, log(Estimate)/log(N), (Estimate^(1 - Order) - 1)/(N^(1 - Order) - 1)))
      U = beta %>% mutate(Estimate = ifelse(Order == 1, log(Estimate)/log(N), (Estimate^(Order - 1) - 1)/(N^(Order - 1) - 1)))
      V = beta %>% mutate(Estimate = (Estimate - 1)/(N - 1))
      S = beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1))

      if(nboot > 1){

        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
        #                     'datatype', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))

        # plan(sequential)
        # plan(multiprocess, workers=7)

        # se = parSapply(cl, 1:nboot, function(i){

        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){

          if (datatype == 'abundance') {

            p_bt = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            f0_hat = nrow(p_bt) - nrow(data)

            distance_matrix_bt = Bootstrap_distance_matrix(rowSums(data), FDdistM, f0_hat, 'abundance')

            data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))

            data_gamma = rowSums(data_bt)
            data_gamma = data_gamma[data_gamma>0]
            data_alpha = as.matrix(data_bt) %>% as.vector

            if (FDtype == 'tau_value'){

              m_gamma = sapply(level, function(i) coverage_to_size(data_gamma, i, datatype='abundance'))
              m_alpha = sapply(level, function(i) coverage_to_size(data_alpha, i, datatype='abundance'))

              output = FD_by_tau(data_bt, distance_matrix_bt, FDtau, level, datatype='abundance', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)
              gamma = output$gamma
              alpha = output$alpha

              beta=gamma/alpha

              ## Observed Beta ##
              output_obs = FD_by_tau(data_bt, distance_matrix_bt, FDtau, level, datatype = 'abundance', by = 'size', m_gamma = sum(data_bt), m_alpha = sum(data_bt))
              gamma_obs = output_obs$gamma
              alpha_obs = output_obs$alpha
              beta_obs = gamma_obs/alpha_obs

            }

            if (FDtype == 'AUC'){

              m_gamma = sapply(level, function(i) coverage_to_size(data_gamma, i, datatype='abundance'))
              m_alpha = sapply(level, function(i) coverage_to_size(data_alpha, i, datatype='abundance'))

              gamma_alpha_over_tau = lapply(cut, function(tau) {

                FD_by_tau(data_bt, distance_matrix_bt, tau, level, datatype = 'abundance', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)

              })

              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)

              left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)

              gamma = colSums((left_limit + right_limit)/2)

              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)

              left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)

              alpha = colSums((left_limit + right_limit)/2)

              beta_over_tau = gamma_over_tau/alpha_over_tau

              left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)

              beta = colSums((left_limit + right_limit)/2)

              ## Observed Beta ##
              gamma_alpha_over_tau = lapply(cut, function(tau) {

                FD_by_tau(data_bt, distance_matrix_bt, tau, level, datatype = 'abundance', by = 'size', m_gamma = sum(data_bt), m_alpha = sum(data_bt))

              })

              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)

              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)

              beta_over_tau = gamma_over_tau/alpha_over_tau

              if (length(q) == 1) beta_over_tau = matrix(beta_over_tau, nrow = 1)

              left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)

              beta_obs = colSums((left_limit + right_limit)/2)

            }

          }

          if (datatype == 'incidence_raw') {

            p_bt = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            f0_hat = nrow(p_bt) - nrow(data_2D[-1,])

            distance_matrix_bt = Bootstrap_distance_matrix(c(n,rowSums(data_gamma_raw)), FDdistM, f0_hat, 'incidence_freq')

            # p_bt = p_bt[rowSums(p_bt)>0,]

            raw = lapply(1:ncol(p_bt), function(j){

              lapply(1:nrow(p_bt), function(i) rbinom(n = n, size = 1, prob = p_bt[i,j])) %>% do.call(rbind,.)

            })

            gamma = Reduce('+', raw)
            gamma[gamma>1] = 1
            data_gamma_raw_bt = gamma
            data_gamma_freq_bt = c(n, rowSums(gamma))

            data_alpha_freq_bt = sapply(raw, rowSums) %>% c(n, .)

            # data_gamma_freq_bt = data_gamma_freq_bt[data_gamma_freq_bt > 0]
            # data_alpha_freq_bt = data_alpha_freq_bt[data_alpha_freq_bt > 0]

            data_2D_bt = apply(sapply(raw, rowSums), 2, function(x) c(n, x)) %>% as.data.frame

            if (FDtype == 'tau_value'){

              m_gamma = sapply(level, function(i) coverage_to_size(data_gamma_freq_bt, i, datatype='incidence_freq'))
              m_alpha = sapply(level, function(i) coverage_to_size(data_alpha_freq_bt, i, datatype='incidence_raw'))

              output = FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, FDtau, level, datatype = 'incidence_raw', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)
              gamma = output$gamma
              alpha = output$alpha

              beta = gamma/alpha

              ## Observed Beta ##
              output_obs = FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, FDtau, level, datatype = 'incidence_raw', by = 'size', m_gamma = data_gamma_freq_bt[1], m_alpha = data_gamma_freq_bt[1])
              gamma_obs = output_obs$gamma
              alpha_obs = output_obs$alpha
              beta_obs = gamma_obs/alpha_obs

            }

            if (FDtype == 'AUC'){

              m_gamma = sapply(level, function(i) coverage_to_size(data_gamma_freq_bt, i, datatype='incidence_freq'))
              m_alpha = sapply(level, function(i) coverage_to_size(data_alpha_freq_bt, i, datatype='incidence_raw'))

              gamma_alpha_over_tau = lapply(cut, function(tau) {

                FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, level, datatype = 'incidence_raw', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)

              })

              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)

              left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)

              gamma = colSums((left_limit + right_limit)/2)

              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)

              left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)

              alpha = colSums((left_limit + right_limit)/2)

              beta_over_tau = gamma_over_tau/alpha_over_tau

              left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)

              beta = colSums((left_limit + right_limit)/2)

              ## Observed Beta ##
              gamma_alpha_over_tau = lapply(cut, function(tau) {

                FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, level, datatype = 'incidence_raw', by = 'size', m_gamma = data_gamma_freq_bt[1], m_alpha = data_gamma_freq_bt[1])

              })

              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)

              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)

              beta_over_tau = gamma_over_tau/alpha_over_tau

              if (length(q) == 1) beta_over_tau = matrix(beta_over_tau, nrow = 1)

              left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)

              beta_obs = colSums((left_limit + right_limit)/2)

            }

          }

          gamma = gamma[under_max_alpha]
          alpha = alpha[under_max_alpha]
          beta = beta[under_max_alpha]

          beta = gamma/alpha

          gamma = c(gamma, rep(0, length(q)))
          alpha = c(alpha, rep(0, length(q)))
          beta = c(beta, beta_obs)

          order = rep(q, each=length(level) + 1)[under_max_alpha]

          beta = data.frame(Estimate = beta, order)

          C = (beta %>% mutate(Estimate = ifelse(order == 1, log(Estimate)/log(N), (Estimate^(1 - order) - 1)/(N^(1 - order) - 1))))$Estimate
          U = (beta %>% mutate(Estimate = ifelse(order == 1, log(Estimate)/log(N), (Estimate^(order - 1) - 1)/(N^(order - 1) - 1))))$Estimate
          V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
          S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate

          beta = beta$Estimate

          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix

          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }) %>% abind(along = 3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start

        # stopCluster(cl)
        # plan(sequential)

      } else {

        se = matrix(0, ncol = 7, nrow = nrow(beta))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)

      }

    }


    se = as.data.frame(se)

    gamma = gamma %>% mutate(s.e. = se$gamma[1:(nrow(se) - length(q))],
                             LCL = Estimate - tmp * se$gamma[1:(nrow(se) - length(q))],
                             UCL = Estimate + tmp * se$gamma[1:(nrow(se) - length(q))],
                             Region = region_name)

    alpha = alpha %>% mutate(s.e. = se$alpha[1:(nrow(se) - length(q))],
                             LCL = Estimate - tmp * se$alpha[1:(nrow(se) - length(q))],
                             UCL = Estimate + tmp * se$alpha[1:(nrow(se) - length(q))],
                             Region = region_name)

    beta = beta %>% mutate(  s.e. = se$beta,
                             LCL = Estimate - tmp * se$beta,
                             UCL = Estimate + tmp * se$beta,
                             Region = region_name)

    C = C %>% mutate(        s.e. = se$C,
                             LCL = Estimate - tmp * se$C,
                             UCL = Estimate + tmp * se$C,
                             Region = region_name)


    U = U %>% mutate(        s.e. = se$U,
                             LCL = Estimate - tmp * se$U,
                             UCL = Estimate + tmp * se$U,
                             Region = region_name)

    V = V %>% mutate(        s.e. = se$V,
                             LCL = Estimate - tmp * se$V,
                             UCL = Estimate + tmp * se$V,
                             Region = region_name)

    S = S %>% mutate(        s.e. = se$S,
                             LCL = Estimate - tmp * se$S,
                             UCL = Estimate + tmp * se$S,
                             Region = region_name)

    if (trunc) {

      gamma = gamma %>% filter(!(Order==0 & round(Size)>2*n))

      alpha = alpha %>% filter(!(Order==0 & round(Size)>2*n))

      beta  = beta  %>% filter(!(Order==0 & round(Size)>2*n))

      C    =  C    %>% filter(!(Order==0 & round(Size)>2*n))

      U    =  U    %>% filter(!(Order==0 & round(Size)>2*n))

      V    =  V    %>% filter(!(Order==0 & round(Size)>2*n))

      S    =  S    %>% filter(!(Order==0 & round(Size)>2*n))

    }

    list(gamma = gamma, alpha = alpha, beta = beta, C = C, U = U, V = V, S = S)

  }

  for_each_region.size = function(data, region_name, N, level) {

    #data
    if (datatype == 'abundance') {

      n = sum(data)
      data_gamma = rowSums(data)
      data_gamma = data_gamma[data_gamma>0]
      data_alpha = as.matrix(data) %>% as.vector

      ref_gamma = n
      ref_alpha = n

    }

    if (datatype == 'incidence_raw') {

      sampling_units = sapply(data, ncol)
      if (length(unique(sampling_units)) > 1) stop("unsupported data structure: the sampling units of all regions must be the same.")
      if (length(unique(sampling_units)) == 1) n = unique(sampling_units)

      gamma = Reduce('+', data)
      gamma[gamma>1] = 1
      data_gamma_raw = gamma
      data_gamma_freq = c(n, rowSums(gamma))

      data_alpha_freq = sapply(data, rowSums) %>% c(n, .)


      data_2D = apply(sapply(data, rowSums), 2, function(x) c(n, x)) %>% as.data.frame

      ref_gamma = n
      ref_alpha = n

    }



    if (diversity == 'TD') {

      if (datatype == 'abundance') {

        gamma = estimate3D(as.numeric(data_gamma), diversity = 'TD', q = q, datatype = "abundance", base = "size", level = level, nboot = nboot) %>% arrange(., SC, Order.q)

        alpha = estimate3D(as.numeric(data_alpha), diversity = 'TD', q = q, datatype = "abundance", base = "size", level = level, nboot = nboot) %>% arrange(., SC, Order.q)

      }

      if (datatype == 'incidence_raw') {

        gamma = estimate3D(as.numeric(data_gamma_freq), diversity = 'TD', q = q, datatype = "incidence_freq", base = "size", level = level, nboot = nboot) %>% arrange(., SC, Order.q)

        alpha = estimate3D(data_alpha_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "size", level = level, nboot = nboot) %>% arrange(., SC, Order.q)

      }

      se = cbind(gamma$s.e., alpha$s.e. / N)
      colnames(se) = c("gamma", "alpha")
      se = as.data.frame(se)
      se[is.na(se)] = 0

      gamma = (cbind(Size = rep(level, each=length(q)), gamma[,-c(1,2,8,9)]) %>%
                 mutate(Method = ifelse(Size>=ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'))
      )[,c(5,3,2,4,1)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_real', 'Size'))


      alpha = (cbind(Size = rep(level, each = length(q)), alpha[,-c(1,2,8,9)]) %>%
                 mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'))
      )[,c(5,3,2,4,1)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_real', 'Size'))

      alpha$Estimate = alpha$Estimate / N

      # beta = alpha
      # beta$Estimate = gamma$Estimate/alpha$Estimate
      #
      # C = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(1-Order) - 1)/(N^(1-Order)-1)))
      # U = beta %>% mutate(Estimate = ifelse(Order==1, log(Estimate)/log(N), (Estimate^(Order-1) - 1)/(N^(Order-1)-1)))
      # V = beta %>% mutate(Estimate = (Estimate-1)/(N-1))
      # S = beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1))

      # if(nboot>1){
      #
      #   # cl = makeCluster(cluster_numbers)
      #   # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
      #   #                     'datatype', 'data_2D'))
      #   # clusterEvalQ(cl, library(tidyverse, magrittr))
      #
      #   # plan(sequential)
      #   # plan(multiprocess)
      #
      #   # se = parSapply(cl, 1:nboot, function(i){
      #
      #   # start = Sys.time()
      #   se = future_lapply(1:nboot, function(i){
      #
      #     if (datatype == 'abundance') {
      #
      #       bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
      #       bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = bootstrap_population[,k]))
      #
      #       bootstrap_data_gamma = rowSums(bootstrap_sample)
      #       bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma > 0]
      #       bootstrap_data_alpha = as.matrix(bootstrap_sample) %>% as.vector
      #       bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha > 0]
      #
      #       gamma = lapply(1:length(level), function(i){
      #         estimate3D(as.numeric(bootstrap_data_gamma), diversity = 'TD', q = q, datatype = "abundance", base = "size", level = level[i], nboot = 0)
      #       }) %>% do.call(rbind,.)
      #
      #       alpha = lapply(1:length(level), function(i){
      #         estimate3D(as.numeric(bootstrap_data_alpha), diversity = 'TD', q = q, datatype = "abundance", base = "size", level = level[i], nboot = 0)
      #       }) %>% do.call(rbind,.)
      #
      #     }
      #
      #     if (datatype == 'incidence_raw') {
      #
      #       bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
      #
      #       raw = lapply(1:ncol(bootstrap_population), function(j){
      #
      #         lapply(1:nrow(bootstrap_population), function(i) rbinom(n = n, size = 1, prob = bootstrap_population[i,j])) %>% do.call(rbind,.)
      #
      #       })
      #
      #       gamma = Reduce('+', raw)
      #       gamma[gamma > 1] = 1
      #       bootstrap_data_gamma_freq = c(n, rowSums(gamma))
      #
      #       bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)
      #
      #       bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq > 0]
      #       bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq > 0]
      #
      #       gamma = lapply(1:length(level), function(i){
      #         estimate3D(bootstrap_data_gamma_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "size", level = level[i], nboot = 0)
      #       }) %>% do.call(rbind,.)
      #
      #       alpha = lapply(1:length(level), function(i){
      #         estimate3D(bootstrap_data_alpha_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "size", level = level[i], nboot = 0)
      #       }) %>% do.call(rbind,.)
      #
      #     }
      #
      #     gamma = gamma$qD
      #
      #     alpha = alpha$qD
      #     alpha = alpha / N
      #
      #     # beta = gamma/alpha
      #     #
      #     # order = rep(q, length(level))
      #     #
      #     # beta = data.frame(Estimate=beta, order)
      #     #
      #     # C = (beta %>% mutate(Estimate = ifelse(order == 1,log(Estimate)/log(N),(Estimate^(1 - order) - 1)/(N^(1 - order) - 1))))$Estimate
      #     # U = (beta %>% mutate(Estimate = ifelse(order == 1,log(Estimate)/log(N),(Estimate^(order - 1) - 1)/(N^(order - 1) - 1))))$Estimate
      #     # V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
      #     # S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate
      #     #
      #     # beta = beta$Estimate
      #     #
      #     # cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
      #     cbind(gamma, alpha) %>% as.matrix
      #
      #     # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
      #   }) %>% abind(along = 3) %>% apply(1:2, sd)
      #   # end = Sys.time()
      #   # end - start
      #
      #   # stopCluster(cl)
      #   # plan(sequential)
      #
      # } else {
      #
      #   # se = matrix(0, ncol = 7, nrow = nrow(beta))
      #   # colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
      #   # se = as.data.frame(se)
      #
      #   se = matrix(0, ncol = 2, nrow = nrow(gamma))
      #   colnames(se) = c("gamma", "alpha")
      #   se = as.data.frame(se)
      #
      # }

    }

    if (diversity == 'PD') {

      if (datatype == 'abundance') {

        aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = data_gamma, rootExtend = T, refT = reft)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

        gamma = iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = level, q = q, nt = n, cal = "PD") %>% t %>% as.data.frame %>%
          set_colnames(q) %>% gather(Order, Estimate) %>%
          mutate(Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma, 'abundance', level), length(q)), Size = rep(level, length(q)))


        aL_table_alpha = c()

        for (i in 1:N){

          x = data[data[,i]>0,i]
          names(x) = rownames(data)[data[,i]>0]

          aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

          aL_table_alpha = rbind(aL_table_alpha, aL_table)

        }


        qPDm = iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = level, q = q, nt = n, cal = "PD")
        qPDm = qPDm/N
        alpha = qPDm %>% t %>% as.data.frame %>%
          set_colnames(q) %>% gather(Order, Estimate) %>%
          mutate(Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha, 'abundance', level), length(q)), Size = rep(level, length(q)))

      }

      if (datatype == 'incidence_raw') {

        aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = as.matrix(data_gamma_raw), datatype = "incidence_raw", refT = reft, rootExtend = T)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
        gamma = iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis=as.matrix(aL_table_gamma$branch.length), m = level, q = q, nt = n, cal = "PD") %>% t %>% as.data.frame %>%
          set_colnames(q) %>% gather(Order, Estimate) %>%
          mutate( Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', level), length(q)), Size = rep(level, length(q)))

        aL_table_alpha = c()

        for (i in 1:N){

          x = data[[i]]

          aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

          aL_table_alpha = rbind(aL_table_alpha, aL_table)

        }

        alpha = (iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = level, q = q, nt = n, cal = "PD")/N) %>% t %>% as.data.frame %>%
          set_colnames(q) %>% gather(Order, Estimate) %>%
          mutate(Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', level), length(q)), Size = rep(level, length(q)))


      }

      gamma = (gamma %>%
                 mutate(Method = ifelse(Size >= ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated')))[,c(2,1,5,3,4)] %>%
        set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_real', 'Size'))

      gamma$Order = as.numeric(gamma$Order)


      alpha = (alpha %>%
                 mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated')))[,c(2,1,5,3,4)] %>%
        set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_real', 'Size'))

      alpha$Order = as.numeric(alpha$Order)

      if (PDtype == 'meanPD') {
        gamma$Estimate = gamma$Estimate/reft
        alpha$Estimate = alpha$Estimate/reft
      }

      # beta = alpha
      # beta$Estimate = gamma$Estimate/alpha$Estimate
      #
      # C = beta %>% mutate(Estimate = ifelse(Order == 1, log(Estimate)/log(N), (Estimate^(1 - Order) - 1)/(N^(1 - Order) - 1)))
      # U = beta %>% mutate(Estimate = ifelse(Order == 1, log(Estimate)/log(N), (Estimate^(Order - 1) - 1)/(N^(Order - 1) - 1)))
      # V = beta %>% mutate(Estimate = (Estimate - 1)/(N - 1))
      # S = beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1))

      if(nboot>1){

        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
        #                     'datatype', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))

        # plan(sequential)
        # plan(multiprocess)

        # se = parSapply(cl, 1:nboot, function(i){

        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){

          if (datatype == 'abundance') {

            tree_bt = PDtree

            bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data)),] %>% matrix(ncol = ncol(data))

            if ( nrow(p_bt) > nrow(data) & sum(unseen_p) > 0 ){

              unseen = unseen_p[which(rowSums(unseen_p) > 0),]
              unseen = matrix(unseen, ncol = ncol(unseen_p), byrow = T)
              p_bt = rbind(p_bt[(1:nrow(data)),], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
              rownames(p_bt) = c(rownames(data), unseen_name)

              bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
              x_bt = bootstrap_sample

              rownames(x_bt) = rownames(p_bt)

              if ( sum(x_bt[-(1:nrow(data)),])>0 ){

                g0_hat = apply(data, 2, function(x){

                  n = sum(x)
                  f1 = sum(x == 1)
                  f2 = sum(x == 2)

                  aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)

                  aL$treeNabu$branch.length = aL$BLbyT[,1]
                  aL = aL$treeNabu %>% select(branch.abun,branch.length)
                  g1 = aL$branch.length[aL$branch.abun == 1] %>% sum
                  g2 = aL$branch.length[aL$branch.abun == 2] %>% sum
                  g0_hat = ifelse( g2 > ((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
                  g0_hat

                })

                te = (x_bt[1:nrow(data),]*(data == 0))>0
                used_length = sapply(1:ncol(data), function(i) {

                  if (sum(te[,i]) == 0) return(0) else {

                    iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x_bt[1:nrow(data),i], rootExtend = T, refT = reft)$treeNabu %>%
                      subset(label %in% names(which(te[,i] == TRUE))) %>% select(branch.length) %>% sum

                  }

                })

                g0_hat = g0_hat - used_length
                g0_hat[g0_hat < 0] = 0

                unseen_sample = x_bt[-(1:nrow(data)),]
                if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol = ncol(x_bt), byrow = T)

                L0_hat = sapply(1:length(g0_hat), function(i) if(sum(unseen_sample[,i] > 0) > 0) (g0_hat[i] / nrow(unseen)) else 0 )

                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow = T) * unseen_sample)) / rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample) == 0)] = 0

                for (i in 1:length(L0_hat)){

                  tip = list(edge = matrix(c(2,1),1,2),
                             tip.label = unseen_name[i],
                             edge.length = L0_hat[i],
                             Nnode = 1)
                  class(tip) = "phylo"

                  tree_bt = tree_bt + tip

                }

              } else {

                x_bt = x_bt[1:nrow(data),]
                p_bt = p_bt[1:nrow(data),]

              }

            } else {

              p_bt = p_bt[1:nrow(data),]
              x_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
              rownames(x_bt) = rownames(data)

            }

            bootstrap_data_gamma = rowSums(x_bt)
            bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma>0]
            bootstrap_data_alpha = as.matrix(x_bt) %>% as.vector
            bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha>0]

            aL = iNEXT.3D:::phyBranchAL_Abu(phylo = tree_bt, data = bootstrap_data_gamma, rootExtend = T, refT = reft)
            aL$treeNabu$branch.length = aL$BLbyT[,1]
            aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

            gamma = as.vector(iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = level, q = q, nt = n, cal = "PD") %>% t)


            aL_table_alpha = c()

            for (i in 1:N){

              # x = x_bt[x_bt[,i]>0,i]
              # names(x) = rownames(p_bt)[x_bt[,i]>0]

              x = x_bt[,i]
              names(x) = rownames(p_bt)
              x = x[x_bt[,i]>0]

              aL = iNEXT.3D:::phyBranchAL_Abu(phylo = tree_bt, data = x, rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[,1]
              aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

              aL_table_alpha = rbind(aL_table_alpha, aL_table)

            }

            alpha = as.vector((iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = level, q = q, nt = n, cal = "PD")/N) %>% t)

          }

          if (datatype == 'incidence_raw') {

            tree_bt = PDtree

            bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data[[1]])),] %>% matrix(ncol=N)

            if ( nrow(p_bt) > nrow(data[[1]]) & sum(unseen_p)>0 ){

              unseen = unseen_p[which(rowSums(unseen_p)>0),]
              unseen = matrix(unseen, ncol = ncol(unseen_p), byrow = T)
              p_bt = rbind(p_bt[(1:nrow(data[[1]])),], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
              rownames(p_bt) = c(rownames(data[[1]]), unseen_name)

              raw = lapply(1:ncol(p_bt), function(j){

                lapply(1:nrow(p_bt), function(i) rbinom(n=n, size=1, prob=p_bt[i,j])) %>% do.call(rbind,.)

              })

              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)

              if ( lapply(1:length(raw), function(i) raw[[i]][-(1:nrow(data[[1]])),]) %>% do.call(sum,.)>0 ){

                R0_hat = sapply(data, function(x){

                  nT = ncol(x)
                  Q1 = sum(rowSums(x)==1)
                  Q2 = sum(rowSums(x)==2)

                  aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)

                  aL$treeNabu$branch.length = aL$BLbyT[,1]
                  aL = aL$treeNabu %>% select(branch.abun,branch.length)
                  R1 = aL$branch.length[aL$branch.abun == 1] %>% sum
                  R2 = aL$branch.length[aL$branch.abun == 2] %>% sum
                  R0_hat = ifelse( R2>((R1*Q2)/(2*Q1)) , ((nT-1)/nT)*(R1^2/(2*R2)) , ((nT-1)/nT)*(R1*(Q1-1)/(2*(Q2+1))) )
                  R0_hat

                })

                te = (sapply(raw, rowSums)[1:nrow(data[[1]]),]*(sapply(data, rowSums) == 0)) > 0
                used_length = sapply(1:N, function(i) {

                  if (sum(te[,i]) == 0) return(0) else {

                    iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = raw[[i]][1:nrow(data[[1]]),], datatype = "incidence_raw", rootExtend = T, refT = reft)$treeNabu %>%
                      subset(label %in% names(which(te[,i] == TRUE))) %>% select(branch.length) %>% sum

                  }

                })

                R0_hat = R0_hat - used_length
                R0_hat[R0_hat < 0] = 0

                unseen_sample = sapply(raw, rowSums)[-(1:nrow(data[[1]])),]
                if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol = N, byrow = T)

                L0_hat = sapply(1:length(R0_hat), function(i) if(sum(unseen_sample[,i] > 0) > 0) (R0_hat[i] / nrow(unseen)) else 0 )

                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow = T) * unseen_sample)) / rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample) == 0)] = 0

                for (i in 1:length(L0_hat)){

                  tip = list(edge = matrix(c(2,1), 1, 2),
                             tip.label = unseen_name[i],
                             edge.length = L0_hat[i],
                             Nnode = 1)
                  class(tip) = "phylo"

                  tree_bt = tree_bt + tip

                }

              } else raw = lapply(raw, function(i) i[1:nrow(data[[1]]),])

            } else {

              p_bt = p_bt[1:nrow(data[[1]]),]
              raw = lapply(1:ncol(p_bt), function(j){

                lapply(1:nrow(p_bt), function(i) rbinom(n = n, size = 1, prob = p_bt[i,j])) %>% do.call(rbind,.)

              })

              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)

            }

            gamma = Reduce('+', raw)
            gamma[gamma > 1] = 1
            bootstrap_data_gamma_raw = gamma
            bootstrap_data_gamma_freq = c(n, rowSums(gamma))

            bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)

            bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq > 0]
            bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq > 0]

            aL = iNEXT.3D:::phyBranchAL_Inc(phylo = tree_bt, data = bootstrap_data_gamma_raw, datatype = "incidence_raw", rootExtend = T, refT = reft)
            aL$treeNabu$branch.length = aL$BLbyT[,1]
            aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

            gamma = as.vector(iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = level, q = q, nt = n, cal="PD") %>% t)


            aL_table_alpha = c()

            for (i in 1:N){

              x = raw[[i]]

              aL = iNEXT.3D:::phyBranchAL_Inc(phylo = tree_bt, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[,1]
              aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)

              aL_table_alpha = rbind(aL_table_alpha, aL_table)

            }

            alpha = as.vector((iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = level, q=q, nt = n, cal = "PD")/N) %>% t)

          }

          if (PDtype == 'meanPD') {
            gamma = gamma/reft
            alpha = alpha/reft
          }

          # beta = gamma/alpha
          #
          # order = rep(q, each = length(level))
          #
          # beta = data.frame(Estimate = beta, order)
          #
          # C = (beta %>% mutate(Estimate = ifelse(order == 1, log(Estimate)/log(N), (Estimate^(1 - order) - 1)/(N^(1 - order) - 1))))$Estimate
          # U = (beta %>% mutate(Estimate = ifelse(order == 1, log(Estimate)/log(N), (Estimate^(order - 1) - 1)/(N^(order - 1) - 1))))$Estimate
          # V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
          # S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate
          #
          # beta = beta$Estimate
          #
          # cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
          cbind(gamma, alpha) %>% as.matrix

          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }) %>% abind(along = 3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start

        # stopCluster(cl)
        # plan(sequential)

      } else {

        # se = matrix(0, ncol = 7, nrow = nrow(beta))
        # colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        # se = as.data.frame(se)

        se = matrix(0, ncol = 2, nrow = nrow(gamma))
        colnames(se) = c("gamma", "alpha")
        se = as.data.frame(se)

      }

    }

    if (diversity == 'FD') {

      FD_by_tau = function(data, distM, tau, datatype, m_gamma, m_alpha) {

        if (datatype == 'abundance') {

          zik = data
          zik = zik[rowSums(data)>0,]

          dij = distM
          dij = dij[rowSums(data)>0, rowSums(data)>0]

          dij[which(dij>tau, arr.ind = T)] = tau
          aik = (1 - dij/tau) %*% as.matrix(zik)
          positive_id = rowSums(aik)>0


          gamma_x = rowSums(zik)[positive_id]
          gamma_a = rowSums(aik)[positive_id]
          gamma_v = gamma_x/gamma_a
          # gamma_a = ifelse(gamma_a < 1, 1, round(gamma_a))

          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))

          gamma = FD.m.est_0(ai_vi_gamma, m_gamma, q, n) %>% as.vector
          # gamma = iNEXT.3D:::FD.m.est(ai_vi_gamma, m_gamma, q, n) %>% as.vector


          alpha_x = as.vector(as.matrix(zik))
          alpha_a = as.vector(aik)
          # alpha_a = ifelse(alpha_a < 1, 1, round(alpha_a))

          alpha_v = alpha_x/alpha_a
          alpha_v = rep(gamma_v,N)

          alpha_v = alpha_v[alpha_a>0]
          alpha_a = alpha_a[alpha_a>0]

          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))

          alpha = (FD.m.est_0(ai_vi_alpha, m_gamma, q, n)/N) %>% as.vector

        }

        if (datatype == 'incidence_raw') {

          data_gamma_freq = data$data_gamma_freq
          data_2D = data$data_2D

          gamma_Y = data_gamma_freq[-1]

          dij = distM
          dij = dij[gamma_Y > 0, gamma_Y > 0]
          gamma_Y = gamma_Y[gamma_Y > 0]

          dij[which(dij>tau, arr.ind = T)] = tau
          gamma_a = (1 - dij/tau) %*% as.matrix(gamma_Y)
          gamma_a[gamma_a > n] = n

          gamma_v = gamma_Y/gamma_a

          # gamma_a = ifelse(gamma_a < 1, 1, round(gamma_a))
          # gamma_a[gamma_a > n] = n

          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))

          gamma = FD.m.est_0(ai_vi_gamma, m_gamma, q, n) %>% as.vector
          # gamma = iNEXT.3D:::FD.m.est(ai_vi_gamma, m_gamma, q, n) %>% as.vector


          alpha_Y = data_2D[-1,]

          dij = distM
          dij = dij[rowSums(data_2D[-1,]) > 0, rowSums(data_2D[-1,])>0]
          alpha_Y = alpha_Y[rowSums(data_2D[-1,])>0,]

          dij[which(dij>tau, arr.ind = T)] = tau
          alpha_a = (1 - dij/tau) %*% as.matrix(alpha_Y)

          # alpha_a = ifelse(alpha_a < 1, 1, round(alpha_a))
          alpha_a[alpha_a > n] = n
          alpha_a = as.vector(alpha_a)

          alpha_v = rep(gamma_v, N)
          alpha_v = alpha_v[alpha_a > 0]
          alpha_a = alpha_a[alpha_a > 0]

          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))
          alpha = (FD.m.est_0(ai_vi_alpha, m_gamma, q, n)/N) %>% as.vector

        }

        return(data.frame(gamma,alpha))

      }

      if (FDtype == 'tau_value'){

        if (datatype == 'abundance') {

          output = FD_by_tau(data, FDdistM, FDtau, datatype = 'abundance', m_gamma = level, m_alpha = level)
          gamma = output$gamma
          alpha = output$alpha

          gamma = data.frame(Size = level, gamma) %>%
            mutate(Method = ifelse(Size >= ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma, 'incidence_freq', level), length(q)))

          alpha = data.frame(Size = level, alpha) %>%
            mutate(Method = ifelse(Size>=ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha, 'incidence_freq', level), length(q)))

          beta = alpha
          beta$alpha = gamma$gamma/alpha$alpha

        }

        if (datatype == 'incidence_raw') {

          output = FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, FDtau, datatype='incidence_raw', m_gamma=level, m_alpha=level)
          gamma = output$gamma
          alpha = output$alpha

          gamma = data.frame(Size = level, gamma) %>%
            mutate(Method = ifelse(Size >= ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', level), length(q)))

          alpha = data.frame(Size = level, alpha) %>%
            mutate(Method = ifelse(Size>=ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', level), length(q)))

          beta = alpha
          beta$alpha = gamma$gamma/alpha$alpha

        }

      }

      if (FDtype == 'AUC'){

        cut = seq(0.00000001, 1, length.out = FDcut_number)
        width = diff(cut)

        if (datatype == 'abundance') {

          gamma_alpha_over_tau = lapply(cut, function(tau) {

            FD_by_tau(data, FDdistM, tau, datatype = 'abundance', m_gamma = level, m_alpha = level)

          })

          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)

          left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)

          gamma = colSums((left_limit + right_limit)/2)

          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)

          left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)

          alpha = colSums((left_limit + right_limit)/2)

          beta_over_tau = gamma_over_tau/alpha_over_tau

          left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)

          beta = colSums((left_limit + right_limit)/2)

          gamma = data.frame(Size = level, gamma) %>%
            mutate(Method = ifelse(Size >= ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma, 'incidence_freq', level), length(q)))

          alpha = data.frame(Size = level, alpha) %>%
            mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha, 'incidence_freq', level), length(q)))

          beta = data.frame(Size = level, beta) %>%
            mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha, 'incidence_freq', level), length(q)))

        }

        if (datatype == 'incidence_raw') {

          gamma_alpha_over_tau = lapply(cut, function(tau) {

            FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, tau, datatype = 'incidence_raw', m_gamma = level, m_alpha = level)

          })

          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)

          left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)

          gamma = colSums((left_limit + right_limit)/2)

          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)

          left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)

          alpha = colSums((left_limit + right_limit)/2)

          beta_over_tau = gamma_over_tau/alpha_over_tau

          left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)

          beta = colSums((left_limit + right_limit)/2)

          gamma = data.frame(Size = level, gamma) %>%
            mutate(Method = ifelse(Size >= ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', level), length(q)))

          alpha = data.frame(Size = level, alpha) %>%
            mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', level), length(q)))

          beta = data.frame(Size = level, beta) %>%
            mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolated'), 'Interpolated'),
                   Order = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', level), length(q)))

        }

      }

      gamma = gamma[,c(2,4,3,5,1)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_real', 'Size'))


      alpha = alpha[,c(2,4,3,5,1)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_real', 'Size'))


      # beta = beta[,c(2,4,3,5,1)] %>% set_colnames(c('Estimate', 'Order', 'Method', 'Coverage_real', 'Size'))
      #
      # C = beta %>% mutate(Estimate = ifelse(Order == 1, log(Estimate)/log(N), (Estimate^(1 - Order) - 1)/(N^(1 - Order) - 1)))
      # U = beta %>% mutate(Estimate = ifelse(Order == 1, log(Estimate)/log(N), (Estimate^(Order - 1) - 1)/(N^(Order - 1) - 1)))
      # V = beta %>% mutate(Estimate = (Estimate - 1)/(N - 1))
      # S = beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1))

      if(nboot > 1){

        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
        #                     'datatype', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))

        # plan(sequential)
        # plan(multiprocess, workers=7)

        # se = parSapply(cl, 1:nboot, function(i){

        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){

          if (datatype == 'abundance') {

            p_bt = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            f0_hat = nrow(p_bt) - nrow(data)

            distance_matrix_bt = Bootstrap_distance_matrix(rowSums(data), FDdistM, f0_hat, 'abundance')

            data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))

            data_gamma = rowSums(data_bt)
            data_gamma = data_gamma[data_gamma>0]
            data_alpha = as.matrix(data_bt) %>% as.vector

            if (FDtype == 'tau_value'){

              output = FD_by_tau(data_bt, distance_matrix_bt, FDtau, datatype='abundance', m_gamma = level, m_alpha = level)
              gamma = output$gamma
              alpha = output$alpha

              beta=gamma/alpha

            }

            if (FDtype == 'AUC'){

              gamma_alpha_over_tau = lapply(cut, function(tau) {

                FD_by_tau(data_bt, distance_matrix_bt, tau, datatype = 'abundance', m_gamma = level, m_alpha = level)

              })

              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)

              left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)

              gamma = colSums((left_limit + right_limit)/2)

              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)

              left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)

              alpha = colSums((left_limit + right_limit)/2)

              beta_over_tau = gamma_over_tau/alpha_over_tau

              left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)

              beta = colSums((left_limit + right_limit)/2)

            }

          }

          if (datatype == 'incidence_raw') {

            p_bt = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            f0_hat = nrow(p_bt) - nrow(data_2D[-1,])

            distance_matrix_bt = Bootstrap_distance_matrix(c(n,rowSums(data_gamma_raw)), FDdistM, f0_hat, 'incidence_freq')

            # p_bt = p_bt[rowSums(p_bt)>0,]

            raw = lapply(1:ncol(p_bt), function(j){

              lapply(1:nrow(p_bt), function(i) rbinom(n = n, size = 1, prob = p_bt[i,j])) %>% do.call(rbind,.)

            })

            gamma = Reduce('+', raw)
            gamma[gamma>1] = 1
            data_gamma_raw_bt = gamma
            data_gamma_freq_bt = c(n, rowSums(gamma))

            data_alpha_freq_bt = sapply(raw, rowSums) %>% c(n, .)

            # data_gamma_freq_bt = data_gamma_freq_bt[data_gamma_freq_bt > 0]
            # data_alpha_freq_bt = data_alpha_freq_bt[data_alpha_freq_bt > 0]

            data_2D_bt = apply(sapply(raw, rowSums), 2, function(x) c(n, x)) %>% as.data.frame

            if (FDtype == 'tau_value'){

              output = FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, FDtau, datatype = 'incidence_raw', m_gamma = level, m_alpha = level)
              gamma = output$gamma
              alpha = output$alpha

              beta = gamma/alpha

            }

            if (FDtype == 'AUC'){

              gamma_alpha_over_tau = lapply(cut, function(tau) {

                FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, datatype = 'incidence_raw', m_gamma = level, m_alpha = level)

              })

              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)

              left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)

              gamma = colSums((left_limit + right_limit)/2)

              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)

              left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)

              alpha = colSums((left_limit + right_limit)/2)

              beta_over_tau = gamma_over_tau/alpha_over_tau

              left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)

              beta = colSums((left_limit + right_limit)/2)

            }

          }

          # beta = gamma/alpha
          #
          # order = rep(q, each=length(level))
          #
          # beta = data.frame(Estimate = beta, order)
          #
          # C = (beta %>% mutate(Estimate = ifelse(order == 1, log(Estimate)/log(N), (Estimate^(1 - order) - 1)/(N^(1 - order) - 1))))$Estimate
          # U = (beta %>% mutate(Estimate = ifelse(order == 1, log(Estimate)/log(N), (Estimate^(order - 1) - 1)/(N^(order - 1) - 1))))$Estimate
          # V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
          # S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate
          #
          # beta = beta$Estimate
          #
          # cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
          cbind(gamma, alpha) %>% as.matrix

          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }) %>% abind(along = 3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start

        # stopCluster(cl)
        # plan(sequential)

      } else {

        # se = matrix(0, ncol = 7, nrow = nrow(beta))
        # colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        # se = as.data.frame(se)

        se = matrix(0, ncol = 2, nrow = nrow(gamma))
        colnames(se) = c("gamma", "alpha")
        se = as.data.frame(se)

      }

    }


    se = as.data.frame(se)

    gamma = gamma %>% mutate(s.e. = se$gamma,
                             LCL = Estimate - tmp * se$gamma,
                             UCL = Estimate + tmp * se$gamma,
                             Region = region_name)

    alpha = alpha %>% mutate(s.e. = se$alpha,
                             LCL = Estimate - tmp * se$alpha,
                             UCL = Estimate + tmp * se$alpha,
                             Region = region_name)

    # beta = beta %>% mutate(  s.e. = se$beta,
    #                          LCL = Estimate - tmp * se$beta,
    #                          UCL = Estimate + tmp * se$beta,
    #                          Region = region_name)
    #
    # C = C %>% mutate(        s.e. = se$C,
    #                          LCL = Estimate - tmp * se$C,
    #                          UCL = Estimate + tmp * se$C,
    #                          Region = region_name)
    #
    #
    # U = U %>% mutate(        s.e. = se$U,
    #                          LCL = Estimate - tmp * se$U,
    #                          UCL = Estimate + tmp * se$U,
    #                          Region = region_name)
    #
    # V = V %>% mutate(        s.e. = se$V,
    #                          LCL = Estimate - tmp * se$V,
    #                          UCL = Estimate + tmp * se$V,
    #                          Region = region_name)
    #
    # S = S %>% mutate(        s.e. = se$S,
    #                          LCL = Estimate - tmp * se$S,
    #                          UCL = Estimate + tmp * se$S,
    #                          Region = region_name)
    #
    # list(gamma = gamma, alpha = alpha, beta = beta, C = C, U = U, V = V, S = S)

    if (datatype == 'incidence_raw') {
      colnames(gamma)[colnames(gamma) == 'Size'] = 'nT'
      colnames(alpha)[colnames(alpha) == 'Size'] = 'nT'
    }

    list(gamma = gamma, alpha = alpha)

  }

  if (base == 'coverage') output = lapply(1:length(data_list), function(i) for_each_region(data = data_list[[i]], region_name = region_names[i], N = Ns[i]))
  if (base == 'size') output = lapply(1:length(data_list), function(i) for_each_region.size(data = data_list[[i]], region_name = region_names[i], N = Ns[i], level = level[[i]]))
  names(output) = region_names

  return(output)

}
iNEXTPDlink = function (data, datatype = "abundance", col.tree = NULL, row.tree = NULL, q = c(0, 1, 2), reftime = NULL, type = "meanPD", endpoint = NULL,
                        knots = 40, size = NULL, nboot = 50, conf = 0.95)
{

  if (!is.null(col.tree)){
    if (sum(c(duplicated(col.tree$tip.label), duplicated(col.tree$node.label[col.tree$node.label !=
                                                                             ""]))) > 0)
      stop("The column phylo tree should not contains duplicated tip or node labels, please remove them.",
           call. = FALSE)
  }
  if (!is.null(row.tree)){
    if (sum(c(duplicated(row.tree$tip.label), duplicated(row.tree$node.label[row.tree$node.label !=
                                                                             ""]))) > 0)
      stop("The row phylo tree should not contains duplicated tip or node labels, please remove them.",
           call. = FALSE)
  }


  DATATYPE <- c("abundance")

  if (is.na(pmatch(datatype, DATATYPE)) == T)
    stop("Invalid datatype", call. = FALSE)
  divtype <- c("PD", "meanPD")
  if (is.na(pmatch(type, divtype)) == T)
    stop("Incorrect type of desired diversity type, please use either PD or meanPD.",
         call. = FALSE)
  if (class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value",
         call. = FALSE)
  if ((conf < 0) | (conf > 1) | (is.numeric(conf) == F))
    stop("conf (confidence level) must be a numerical value between 0 and 1, We use \"conf\" = 0.95 to calculate!",
         call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot) == F))
    stop("nboot must be a nonnegative integer, We use \"nboot\" = 50 to calculate!",
         call. = FALSE)


  if (length(data) == 1) {
    dat = as.matrix(as.vector(t(data[[1]])))
    region_names = if (is.null(names(data)))
      paste0("region_", 1)else
        names(data)
    colnames(dat) = region_names
    rownames(dat) = paste0(rep(colnames(data[[1]]),nrow(data[[1]])),"*",rep(rownames(data[[1]]),rep(ncol(data[[1]]),nrow(data[[1]]))))
    data1 = dat

  }else {
    region_names = if (is.null(names(data)))
      paste0("region_", 1:length(data))else names(data)
    data2 = lapply(data, function(i) {

      int_name = paste0(rep(colnames(i),nrow(i)),"*",rep(rownames(i),rep(ncol(i),nrow(i))))
      i = as.vector(t(i))

      i = data.frame(species = int_name, i)
      return(i)
    })
    data1 = data2[[1]]
    for (i in 2:length(data2)) {
      data1 = data.frame(full_join(data1, data2[[i]],
                                   by = "species"))
    }
    data1[is.na(data1)] = 0
    rownames(data1) = data1$species
    data1 = data1[!colnames(data1) == "species"]
    colnames(data1) = region_names


  }


  data1 <- data1[rowSums(data1) > 0, , drop = FALSE]
  pool.name <- rownames(data1)
  mydata = list()
  if (is.null(colnames(data1))) {
    colnames(data1) <- paste0("assemblage", 1:ncol(data1))
  }
  mydata <- lapply(1:ncol(data1), function(i) {
    x <- data1[, i]
    names(x) <- pool.name
    x
  })
  names(mydata) = colnames(data1)



  # tip <- tree$tip.label[-match(pool.name, tree$tip.label)]
  # mytree <- drop.tip(tree, tip)
  # H_max = ifelse(is.ultrametric(mytree), get.rooted.tree.height(mytree),
  #                max(ape::node.depth.edgelength(mytree)))
  # if (is.null(reftime))
  #   reftime <- H_max
  # else reftime <- reftime
  reftime <- sort(unique(reftime))
  # if (sum(reftime <= 0) > 0) {
  #   stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",
  #        call. = FALSE)
  # }

  if (length(knots) != length(mydata))
    knots <- rep(knots, length(mydata))

  if (is.null(size)) {
    if (is.null(endpoint)) {
      if (datatype == "abundance") {
        endpoint <- sapply(mydata, function(x) 2 * sum(x))
      }
      else if (datatype == "incidence_raw") {
        endpoint <- sapply(mydata, function(x) 2 * ncol(x))
      }
    }
    else {
      if (length(endpoint) != length(mydata)) {
        endpoint <- rep(endpoint, length(mydata))
      }
    }
    size <- lapply(1:length(mydata), function(i) {
      if (datatype == "abundance") {
        ni <- sum(mydata[[i]])
      }
      else if (datatype == "incidence_raw") {
        ni <- ncol(mydata[[i]])
      }
      if (endpoint[i] <= ni) {
        mi <- floor(seq(1, endpoint[i], length.out = knots[i]))
      }
      else {
        mi <- floor(c(seq(1, ni, length.out = floor(knots[i]/2)),
                      seq(ni + 1, endpoint[i], length.out = knots[i] -
                            floor(knots[i]/2))))
      }
      unique(mi)
    })
  }
  else {
    if (class(size) == "numeric" | class(size) == "integer" |
        class(size) == "double") {
      size <- list(size = size)
    }
    if (length(size) != length(mydata))
      size <- lapply(1:length(mydata), function(x) size[[1]])
    size <- lapply(1:length(mydata), function(i) {
      if (datatype == "abundance") {
        ni <- sum(mydata[[i]])
      }
      else if (datatype == "incidence_raw") {
        ni <- ncol(mydata[[i]])
      }
      if (sum(size[[i]] == ni) == 0)
        mi <- sort(c(ni, size[[i]]))
      else mi <- size[[i]]
      unique(mi)
    })
  }
  FUN <- function(e) {
    if (class(data) == "list") {
      inextPDlink(datalist = data, datatype = datatype,
                  col.tree = col.tree,row.tree = row.tree, q = q, reft = reftime, m = size,
                  cal = type, nboot = nboot, conf = conf, unconditional_var = TRUE)
    }
    else {
      return(NULL)
    }
  }
  out <- tryCatch(FUN(e), error = function(e) {
    return()
  })
  index <- AO.link(data = data, diversity = "PD",row.tree = row.tree, col.tree  = col.tree,
                   q = c(0, 1, 2), datatype = datatype, PDtype = type, nboot = nboot, conf = 0.95)
  index = index[order(index$Network), ]
  LCL <- index$qD.LCL[index$Method == "Estimated"]
  UCL <- index$qD.UCL[index$Method == "Estimated"]
  index <- dcast(index, formula = Network + Order.q ~ Method,
                 value.var = "qD")
  index <- cbind(index, se = (UCL - index$Estimated)/qnorm(1 -
                                                             (1 - conf)/2), LCL, UCL)
  if (nboot > 0)
    index$LCL[index$LCL < index$Empirical & index$Order.q ==
                0] <- index$Empirical[index$LCL < index$Empirical &
                                        index$Order.q == 0]
  index$Order.q <- c("Species richness", "Shannon diversity",
                     "Simpson diversity")
  colnames(index) <- c("Assemblage", "Phylogenetic Diversity",
                       "Phylogenetic Observed", "Phylogenetic Estimator", "s.e.",
                       "LCL", "UCL")
  info <- DataInfo.link(data = data, diversity = "PD", datatype = datatype,
                        row.tree = row.tree, col.tree = col.tree)

  return(list(PDInfo = info, PDiNextEst = out, PDAsyEst = index))
}




inextPDlink = function (datalist, datatype,col.tree,row.tree, q, reft, m, cal, nboot,
                        conf = 0.95, unconditional_var = TRUE)
{

  nms <- names(datalist)
  qtile <- qnorm(1 - (1 - conf)/2)
  rt = c()
  for(i in 1:length(datalist)){
    aL <- create.aili(data = datalist[[i]],row.tree = row.tree,col.tree = col.tree)

    x <- as.vector(datalist[[i]]) %>% .[. > 0]
    n <- sum(x)
    rt[i] = sum(aL$branch.abun*aL$branch.length/n)
  }
  reft = max(rt)
  if (datatype == "abundance") {
    Estoutput <- lapply(1:length(datalist), function(i) {

      aL <- create.aili(data = datalist[[i]],row.tree = row.tree,col.tree = col.tree)

      x <- as.vector(datalist[[i]]) %>% .[. > 0]
      n <- sum(x)
      aL = filter(aL,branch.abun >0)
      qPDm <- iNEXT.3D:::PhD.m.est(ai = aL$branch.abun%>%as.matrix(),
                                   Lis = aL$branch.length%>%as.matrix(), m = m[[i]], q = q, nt = n, reft = reft,
                                   cal = cal) %>% as.numeric()
      covm = iNEXT.3D:::Coverage(x, datatype, m[[i]])
      if (unconditional_var) {
        goalSC <- unique(covm)
        lis = aL$branch.length%>%as.matrix()
        colnames(lis) = paste0("T",reft)
        qPD_unc <- unique(iNEXT.3D:::invChatPD_abu(x = x, ai = aL$branch.abun%>%as.matrix(),
                                                   Lis = lis, q = q, Cs = goalSC, n = n,
                                                   reft = reft, cal = cal))
        qPD_unc$Method[round(qPD_unc$m) == n] = "Observed"
      }
      if (nboot > 1) {
        boot.sam <- sample.boot.phy(datalist[[i]],nboot,row.tree = row.tree,col.tree = col.tree)


        if (unconditional_var) {
          ses <- sapply(1:nboot, function(B) {

            ai_B <- boot.sam[[B]]$branch.abun %>% as.matrix()
            Li_b <- boot.sam[[B]]$branch.length %>% as.matrix()
            colnames(Li_b) = paste0("T",reft)
            isn0 <- ai_B > 0
            qPDm_b <- iNEXT.3D:::PhD.m.est(ai = ai_B[isn0], Lis = Li_b[isn0,
                                                                       , drop = F], m = m[[i]], q = q, nt = n,
                                           reft = reft, cal = cal) %>% as.numeric()
            covm_b <- iNEXT.3D:::Coverage(ai_B[isn0], datatype, m[[i]])
            qPD_unc_b <- unique(iNEXT.3D:::invChatPD_abu(x = ai_B[isn0], ai = ai_B[isn0], Lis = Li_b[isn0,
                                                                                                     , drop = F], q = q, Cs = goalSC, n = n,
                                                         reft = reft, cal = cal))$qPD
            return(c(qPDm_b, covm_b, qPD_unc_b))
          }) %>% apply(., 1, sd)
        }
        else {
          ses <- sapply(1:nboot, function(B) {
            ai_B <- boot.sam[[B]]$branch.abun %>% as.matrix()
            Li_b <- boot.sam[[B]]$branch.length %>% as.matrix()
            colnames(Li_b) = paste0("T",reft)
            isn0 <- ai_B > 0
            qPDm_b <- iNEXT.3D:::PhD.m.est(ai = ai_B[isn0], Lis = Li_b[isn0,
                                                                       , drop = F], m = m[[i]], q = q, nt = n,
                                           reft = reft, cal = cal) %>% as.numeric()
            covm_b <- iNEXT.3D:::Coverage(ai_B[isn0], datatype, m[[i]])
            return(c(qPDm_b, covm_b))
          }) %>% apply(., 1, sd)
        }
      }
      else {
        if (unconditional_var) {
          ses <- rep(NA, length(c(qPDm, covm, qPD_unc$qPD)))
        }
        else {
          ses <- rep(NA, length(c(qPDm, covm)))
        }
      }
      ses_pd <- ses[1:length(qPDm)]
      ses_cov <- ses[(length(qPDm) + 1):(length(qPDm) +
                                           length(covm))]
      m_ <- rep(m[[i]], each = length(q) * length(reft))
      method <- ifelse(m[[i]] > n, "Extrapolation", ifelse(m[[i]] <
                                                             n, "Rarefaction", "Observed"))
      method <- rep(method, each = length(q) * length(reft))
      orderq <- rep(q, length(reft) * length(m[[i]]))
      SC_ <- rep(covm, each = length(q) * length(reft))
      SC.se <- rep(ses_cov, each = length(q) * length(reft))
      SC.LCL_ <- rep(covm - qtile * ses_cov, each = length(q) *
                       length(reft))
      SC.UCL_ <- rep(covm + qtile * ses_cov, each = length(q) *
                       length(reft))
      reft_ <- rep(rep(reft, each = length(q)), length(m[[i]]))
      out_m <- tibble(Assemblage = nms[i], m = m_, Method = method,
                      Order.q = orderq, qPD = qPDm, s.e. = ses_pd,
                      qPD.LCL = qPDm - qtile * ses_pd, qPD.UCL = qPDm +
                        qtile * ses_pd, SC = SC_, SC.s.e. = SC.se,
                      SC.LCL = SC.LCL_, SC.UCL = SC.UCL_, Reftime = reft_,
                      Type = cal) %>% arrange(Reftime, Order.q, m)
      out_m$qPD.LCL[out_m$qPD.LCL < 0] <- 0
      out_m$SC.LCL[out_m$SC.LCL < 0] <- 0
      out_m$SC.UCL[out_m$SC.UCL > 1] <- 1
      if (unconditional_var) {
        ses_pd_unc <- ses[-(1:(length(qPDm) + length(covm)))]
        out_C <- qPD_unc %>% mutate(qPD.LCL = qPD -
                                      qtile * ses_pd_unc, qPD.UCL = qPD + qtile *
                                      ses_pd_unc, s.e. = ses_pd_unc, Type = cal,
                                    Assemblage = nms[i])
        id_C <- match(c("Assemblage", "goalSC", "SC",
                        "m", "Method", "Order.q", "qPD", "s.e.", "qPD.LCL",
                        "qPD.UCL", "Reftime", "Type"), names(out_C),
                      nomatch = 0)
        out_C <- out_C[, id_C] %>% arrange(Reftime,
                                           Order.q, m)
        out_C$qPD.LCL[out_C$qPD.LCL < 0] <- 0
      }
      else {
        out_C <- NULL
      }
      return(list(size_based = out_m, coverage_based = out_C))
    })
  }
  else if (datatype == "incidence_raw") {
    Estoutput <- lapply(1:length(datalist), function(i) {
      aL <- phyBranchAL_Inc(phylo = phylotr, data = datalist[[i]],
                            datatype, refT = reft)
      x <- datalist[[i]] %>% .[rowSums(.) > 0, colSums(.) >
                                 0]
      n <- ncol(x)
      qPDm <- iNEXT.3D:::PhD.m.est(ai = aL$treeNabu$branch.abun,
                                   Lis = aL$BLbyT, m = m[[i]], q = q, nt = n, reft = reft,
                                   cal = cal) %>% as.numeric()
      covm = Coverage(x, datatype, m[[i]])
      if (unconditional_var) {
        goalSC <- unique(covm)
        qPD_unc <- unique(invChatPD_inc(x = rowSums(x),
                                        ai = aL$treeNabu$branch.abun, Lis = aL$BLbyT,
                                        q = q, Cs = goalSC, n = n, reft = reft, cal = cal))
        qPD_unc$Method[round(qPD_unc$nt) == n] = "Observed"
      }
      if (nboot > 1) {
        Boots <- Boots.one(phylo = phylotr, aL$treeNabu,
                           datatype, nboot, reft, aL$BLbyT, n)
        Li_b <- Boots$Li
        refinfo <- colnames(Li_b)
        Li_b <- sapply(1:length(reft), function(l) {
          tmp <- Li_b[, l]
          tmp[tmp > reft[l]] <- reft[l]
          tmp
        })
        colnames(Li_b) <- refinfo
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip", nrow(x) + f0), rep("Inode",
                                                    nrow(Li_b) - nrow(x) - f0))
        if (unconditional_var) {
          ses <- sapply(1:nboot, function(B) {
            ai_B <- Boots$boot_data[, B]
            isn0 <- ai_B > 0
            qPDm_b <- iNEXT.3D:::PhD.m.est(ai = ai_B[isn0], Lis = Li_b[isn0,
                                                                       , drop = F], m = m[[i]], q = q, nt = n,
                                           reft = reft, cal = cal) %>% as.numeric()
            covm_b = Coverage(c(n, ai_B[isn0 & tgroup_B ==
                                          "Tip"]), "incidence_freq", m[[i]])
            qPD_unc_b <- unique(invChatPD_inc(x = ai_B[isn0 &
                                                         tgroup_B == "Tip"], ai = ai_B[isn0], Lis = Li_b[isn0,
                                                                                                         , drop = F], q = q, Cs = goalSC, n = n,
                                              reft = reft, cal = cal))$qPD
            return(c(qPDm_b, covm_b, qPD_unc_b))
          }) %>% apply(., 1, sd)
        }
        else {
          ses <- sapply(1:nboot, function(B) {
            ai_B <- Boots$boot_data[, B]
            isn0 <- ai_B > 0
            qPDm_b <- iNEXT.3D:::PhD.m.est(ai = ai_B[isn0], Lis = Li_b[isn0,
                                                                       , drop = F], m = m[[i]], q = q, nt = n,
                                           reft = reft, cal = cal) %>% as.numeric()
            covm_b = Coverage(c(n, ai_B[isn0 & tgroup_B ==
                                          "Tip"]), "incidence_freq", m[[i]])
            return(c(qPDm_b, covm_b))
          }) %>% apply(., 1, sd)
        }
      }
      else {
        if (unconditional_var) {
          ses <- rep(NA, length(c(qPDm, covm, qPD_unc$qPD)))
        }
        else {
          ses <- rep(NA, length(c(qPDm, covm)))
        }
      }
      ses_pd <- ses[1:length(qPDm)]
      ses_cov <- ses[(length(qPDm) + 1):(length(qPDm) +
                                           length(covm))]
      m_ <- rep(m[[i]], each = length(q) * length(reft))
      method <- ifelse(m[[i]] > n, "Extrapolation", ifelse(m[[i]] <
                                                             n, "Rarefaction", "Observed"))
      method <- rep(method, each = length(q) * length(reft))
      orderq <- rep(q, length(reft) * length(m[[i]]))
      SC_ <- rep(covm, each = length(q) * length(reft))
      SC.se <- rep(ses_cov, each = length(q) * length(reft))
      SC.LCL_ <- rep(covm - qtile * ses_cov, each = length(q) *
                       length(reft))
      SC.UCL_ <- rep(covm + qtile * ses_cov, each = length(q) *
                       length(reft))
      reft_ = rep(rep(reft, each = length(q)), length(m[[i]]))
      out_m <- tibble(Assemblage = nms[i], nt = m_, Method = method,
                      Order.q = orderq, qPD = qPDm, s.e. = ses_pd,
                      qPD.LCL = qPDm - qtile * ses_pd, qPD.UCL = qPDm +
                        qtile * ses_pd, SC = SC_, SC.s.e. = SC.se,
                      SC.LCL = SC.LCL_, SC.UCL = SC.UCL_, Reftime = reft_,
                      Type = cal) %>% arrange(Reftime, Order.q, nt)
      out_m$qPD.LCL[out_m$qPD.LCL < 0] <- 0
      out_m$SC.LCL[out_m$SC.LCL < 0] <- 0
      out_m$SC.UCL[out_m$SC.UCL > 1] <- 1
      if (unconditional_var) {
        ses_pd_unc <- ses[-(1:(length(qPDm) + length(covm)))]
        out_C <- qPD_unc %>% mutate(qPD.LCL = qPD -
                                      qtile * ses_pd_unc, qPD.UCL = qPD + qtile *
                                      ses_pd_unc, s.e. = ses_pd_unc, Type = cal,
                                    Assemblage = nms[i])
        id_C <- match(c("Assemblage", "goalSC", "SC",
                        "nt", "Method", "Order.q", "qPD", "s.e.",
                        "qPD.LCL", "qPD.UCL", "Reftime", "Type"),
                      names(out_C), nomatch = 0)
        out_C <- out_C[, id_C] %>% arrange(Reftime,
                                           Order.q, nt)
        out_C$qPD.LCL[out_C$qPD.LCL < 0] <- 0
      }
      else {
        out_C <- NULL
      }
      return(list(size_based = out_m, coverage_based = out_C))
    })
  }
  if (unconditional_var) {
    ans <- list(size_based = do.call(rbind, lapply(Estoutput,
                                                   function(x) x$size_based)), coverage_based = do.call(rbind,
                                                                                                        lapply(Estoutput, function(x) x$coverage_based)))
  }
  else {
    ans <- list(size_based = do.call(rbind, lapply(Estoutput,
                                                   function(x) x$size_based)))
  }
  return(ans)
}

iNEXTlinkFD = function (data, row.distM = NULL, col.distM = NULL , datatype = "abundance", q = c(0, 1, 2),
                        endpoint = NULL, knots = 40, size = NULL, conf = 0.95, nboot = 50,
                        threshold = NULL, nT = NULL)
{

  if (length(data) == 1) {
    dat = as.matrix(as.vector(t(data[[1]])))
    region_names = if (is.null(names(data)))
      paste0("region_", 1)else
        names(data)
    colnames(dat) = region_names
    rownames(dat) = paste0(rep(colnames(data[[1]]),3),".",rep(rownames(data[[1]]),rep(ncol(data[[1]]),nrow(data[[1]]))))
    data1 = dat
    if(is.null(row.distM)){
      rdd = matrix(1,ncol = nrow(data[[1]]),nrow = nrow(data[[1]]))
      diag(rdd) = 0
      rownames(rdd) = rownames(data[[1]])
      colnames(rdd) = rownames(data[[1]])
      row.distM =  rdd}
    if(is.null(col.distM)){
      cdd = matrix(1,ncol = ncol(data[[1]]),nrow = ncol(data[[1]]))
      diag(cdd) = 0
      rownames(cdd) = colnames(data[[1]])
      colnames(cdd) = colnames(data[[1]])
      col.distM =  cdd}
  }else {
    region_names = if (is.null(names(data)))
      paste0("region_", 1:length(data))else names(data)
    data2 = lapply(data, function(i) {

      int_name = paste0(rep(colnames(i),3),".",rep(rownames(i),rep(ncol(i),nrow(i))))
      i = as.vector(t(i))

      i = data.frame(species = int_name, i)
      return(i)
    })
    data1 = data2[[1]]
    for (i in 2:length(data2)) {
      data1 = data.frame(full_join(data1, data2[[i]],
                                   by = "species"))
    }
    data1[is.na(data1)] = 0
    rownames(data1) = data1$species
    data1 = data1[!colnames(data1) == "species"]
    colnames(data1) = region_names
    row_sp = c()
    col_sp = c()
    for(i in 1:length(data)){
      row_sp = c(row_sp,rownames(data[[i]]))
      col_sp = c(col_sp,colnames(data[[i]]))

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

  }

  row.distM = as.matrix(row.distM)
  col.distM = as.matrix(col.distM)
  # if (datatype == "incidence_raw") {
  #   data = as.incfreq(data, nT = nT)
  #   datatype = "incidence_freq"
  # }
  distM =  1-(1-row.distM[rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM))),rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM)))])*(1-col.distM[rep(1:nrow(col.distM),nrow(row.distM)),rep(1:nrow(col.distM),nrow(row.distM))])
  distM_name = paste0(rep(colnames(col.distM),3),".",rep(rownames(row.distM),rep(ncol(col.distM),nrow(row.distM))))
  colnames(distM) = distM_name
  rownames(distM) = distM_name
  # DATATYPE <- c("abundance", "incidence_freq")
  # if (is.na(pmatch(datatype, DATATYPE)) == T)
  #   stop("invalid datatype", call. = FALSE)
  if (class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value",
         call. = FALSE)
  if ((conf < 0) | (conf > 1) | (is.numeric(conf) == F))
    stop("conf (confidence level) must be a numerical value between 0 and 1, We use \"conf\" = 0.95 to calculate!",
         call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot) == F))
    stop("nboot must be a nonnegative integer, We use \"nboot\" = 50 to calculate!",
         call. = FALSE)
  # if (datatype == "incidence_freq") {
  #   nT <- data[1, ]
  #   data <- data[-1, , drop = FALSE]
  # }

  distM = distM[rownames(distM) %in% rownames(data1), colnames(distM) %in%
                  rownames(data1)]

  if (nrow(data1) != nrow(distM))
    stop("The number of species in data should equal to that in distance matrix",
         call. = FALSE)

  if (sum(rownames(data1) %in% rownames(distM)) != nrow(distM))
    stop("Data and distance matrix contain unmatched species",
         call. = FALSE)

  order_sp <- match(rownames(data1), rownames(distM))
  distM <- distM[order_sp, order_sp]



  # if (datatype == "incidence_freq") {
  #   data <- rbind(nT, data)
  # }
  name_sp <- rownames(data1)

  dat <- lapply(1:ncol(data1), function(k) {
    x <- data1[, k]
    names(x) <- name_sp
    x
  })
  if (is.null(colnames(data1))) {
    names(dat) <- paste0("site", 1:length(dat))
  }else {
    names(dat) = colnames(data1)
  }
  if (is.null(threshold)) {
    if (datatype == "abundance") {
      tmp = sapply(dat, function(x) x) %>% apply(., 1,
                                                 sum)
      tmp <- matrix(tmp/sum(tmp), ncol = 1)
    }
    else if (datatype == "incidence_freq") {
      tmp = sapply(dat, function(x) x) %>% apply(., 1,
                                                 sum)
      tmp <- matrix(tmp[-1]/sum(tmp[-1]), ncol = 1)
    }
    dmean <- sum((tmp %*% t(tmp)) * distM)
    dmin <- min(distM[distM > 0])
    threshold <- dmean
  }else if (sum(threshold < 0) > 0 | sum(threshold > 1) > 0) {
    stop("Threshold must be a number between 0 and 1. Use NULL to set it to dmean/2.",
         call. = FALSE)
  }
  if (length(knots) != length(dat))
    knots <- rep(knots, length(dat))
  if (is.null(size)) {
    if (is.null(endpoint)) {
      if (datatype == "abundance") {
        endpoint <- sapply(dat, function(x) 2 * sum(x))
      }
      else if (datatype == "incidence_freq") {
        endpoint <- sapply(dat, function(x) 2 * x[1])
      }
    }
    else {
      if (length(endpoint) != length(dat)) {
        endpoint <- rep(endpoint, length(dat))
      }
    }
    size <- lapply(1:length(dat), function(i) {
      if (datatype == "abundance") {
        ni <- sum(dat[[i]])
      }
      else if (datatype == "incidence_freq") {
        ni <- dat[[i]][1]
      }
      if (endpoint[i] <= ni) {
        mi <- floor(seq(1, endpoint[i], length.out = knots[i]))
      }
      else {
        mi <- floor(c(seq(1, ni, length.out = floor(knots[i]/2)),
                      seq(ni + 1, endpoint[i], length.out = knots[i] -
                            floor(knots[i]/2))))
      }
      unique(mi)
    })
  }else {
    if (class(size) == "numeric" | class(size) == "integer" |
        class(size) == "double") {
      size <- list(size = size)
    }
    if (length(size) != length(dat))
      size = lapply(1:length(dat), function(x) size[[1]])
    size <- lapply(1:length(dat), function(i) {
      if (datatype == "abundance")
        ni <- sum(dat[[i]])
      else ni <- (dat[[i]])[1]
      if (sum(size[[i]] == ni) == 0)
        mi <- sort(c(ni, size[[i]]))
      else mi <- size[[i]]
      unique(mi)
    })
  }
  FUN <- function(e) {
    if (class(dat) == "list") {
      temp1 = iNEXT.3D:::iNextFD(datalist = dat, dij = distM, q = q,
                                 datatype = datatype, tau = threshold, nboot = nboot,
                                 conf = conf, m = size)
      temp1$qFD.LCL[temp1$qFD.LCL < 0] <- 0
      temp1$SC.LCL[temp1$SC.LCL < 0] <- 0
      temp1$SC.UCL[temp1$SC.UCL > 1] <- 1
      if (datatype == "incidence_freq")
        colnames(temp1)[colnames(temp1) == "m"] = "nt"
      temp2 <- lapply(1:length(dat), function(i) iNEXT.3D:::invChatFD(datalist = dat[i],
                                                                      dij = distM, q = q, datatype = datatype, level = iNEXT.3D:::Coverage(data = dat[[i]],
                                                                                                                                           datatype = datatype, m = size[[i]]), nboot = nboot,
                                                                      conf = conf, tau = threshold)) %>% do.call(rbind,
                                                                                                                 .)
      temp2$qFD.LCL[temp2$qFD.LCL < 0] <- 0
      if (datatype == "incidence_freq")
        colnames(temp2)[colnames(temp2) == "m"] = "nt"
      temp1$Type = "FD"
      temp2$Type = "FD"
      ans <- list(size_based = temp1, coverage_based = temp2)
      return(ans)
    }
    else {
      return(NULL)
    }
  }
  out <- tryCatch(FUN(e), error = function(e) {
    return()
  })
  index <- iNEXT.3D:::AO3D(data = data1, diversity = "FD", FDdistM = distM, q = c(0, 1, 2), datatype = datatype, nboot = nboot, conf = 0.95,
                           FDtype = 'tau_values', FDtau = NULL)
  index <- index %>% arrange(., Assemblage)
  LCL <- index$qFD.LCL[index$Method == "Asymptotic"]
  UCL <- index$qFD.UCL[index$Method == "Asymptotic"]
  index <- dcast(index, formula = Assemblage + Order.q ~ Method,
                 value.var = "qFD")
  index <- cbind(index, se = (UCL - index$Asymptotic)/qnorm(1 -
                                                              (1 - conf)/2), LCL, UCL)
  if (nboot > 0)
    index$LCL[index$LCL < index$Empirical & index$Order.q ==
                0] <- index$Empirical[index$LCL < index$Empirical &
                                        index$Order.q == 0]
  index$Order.q <- c("Species richness", "Shannon diversity",
                     "Simpson diversity")
  index[, 3:4] = index[, 4:3]
  colnames(index) <- c("Assemblage", "Functional Diversity",
                       "Functional Observed", "Functional Estimator", "s.e.",
                       "LCL", "UCL")
  info <- iNEXT.3D:::DataInfo3D(data1, diversity = "FD", datatype = datatype,
                                FDdistM = distM, FDtype = "tau_values", FDtau = threshold,
                                nT = nT)
  info$n = lapply(dat, function(x) sum(x))
  info$SC = lapply(dat, function(x) {
    n = sum(x)
    f1 = sum(x == 1)
    f2 = sum(x == 2)
    f0.hat <- ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2,
                     (n - 1)/n * f1^2/2/f2)
    A <- ifelse(f1 > 0, n * f0.hat/(n * f0.hat + f1), 1)
    Chat <- round(1 - f1/n * A, 4)
  })

  return(list(FDInfo = info, FDiNextEst = out, FDAsyEst = index))
}


iNEXTlinkAUC = function (data, row.distM = NULL, col.distM = NULL , datatype = "abundance", q = c(0, 1, 2),
                         endpoint = NULL, knots = 40, size = NULL, conf = 0.95, nboot = 50,
                         threshold = NULL, nT = NULL)
{
  if (length(data) == 1) {
    dat = as.matrix(as.vector(t(data[[1]])))
    region_names = if (is.null(names(data)))
      paste0("region_", 1)else
        names(data)
    colnames(dat) = region_names
    rownames(dat) = paste0(rep(colnames(data[[1]]),3),".",rep(rownames(data[[1]]),rep(ncol(data[[1]]),nrow(data[[1]]))))
    data1 = dat
    if(is.null(row.distM)){
      rdd = matrix(1,ncol = nrow(data[[1]]),nrow = nrow(data[[1]]))
      diag(rdd) = 0
      rownames(rdd) = rownames(data[[1]])
      colnames(rdd) = rownames(data[[1]])
      row.distM =  rdd}
    if(is.null(col.distM)){
      cdd = matrix(1,ncol = ncol(data[[1]]),nrow = ncol(data[[1]]))
      diag(cdd) = 0
      rownames(cdd) = colnames(data[[1]])
      colnames(cdd) = colnames(data[[1]])
      col.distM =  cdd}
  }else {
    region_names = if (is.null(names(data)))
      paste0("region_", 1:length(data))else names(data)
    data2 = lapply(data, function(i) {

      int_name = paste0(rep(colnames(i),3),".",rep(rownames(i),rep(ncol(i),nrow(i))))
      i = as.vector(t(i))

      i = data.frame(species = int_name, i)
      return(i)
    })
    data1 = data2[[1]]
    for (i in 2:length(data2)) {
      data1 = data.frame(full_join(data1, data2[[i]],
                                   by = "species"))
    }
    data1[is.na(data1)] = 0
    rownames(data1) = data1$species
    data1 = data1[!colnames(data1) == "species"]
    colnames(data1) = region_names
    row_sp = c()
    col_sp = c()
    for(i in 1:length(data)){
      row_sp = c(row_sp,rownames(data[[i]]))
      col_sp = c(col_sp,colnames(data[[i]]))

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

  }
  row.distM = as.matrix(row.distM)
  col.distM = as.matrix(col.distM)
  # if (datatype == "incidence_raw") {
  #   data = as.incfreq(data, nT = nT)
  #   datatype = "incidence_freq"
  # }
  distM =  1-(1-row.distM[rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM))),rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM)))])*(1-col.distM[rep(1:nrow(col.distM),nrow(row.distM)),rep(1:nrow(col.distM),nrow(row.distM))])
  distM_name = paste0(rep(colnames(col.distM),3),".",rep(rownames(row.distM),rep(ncol(col.distM),nrow(row.distM))))
  colnames(distM) = distM_name
  rownames(distM) = distM_name

  # DATATYPE <- c("abundance", "incidence_freq")
  # if (is.na(pmatch(datatype, DATATYPE)) == T)
  #   stop("invalid datatype", call. = FALSE)
  if (class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value",
         call. = FALSE)
  if ((conf < 0) | (conf > 1) | (is.numeric(conf) == F))
    stop("conf (confidence level) must be a numerical value between 0 and 1, We use \"conf\" = 0.95 to calculate!",
         call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot) == F))
    stop("nboot must be a nonnegative integer, We use \"nboot\" = 50 to calculate!",
         call. = FALSE)
  # if (datatype == "incidence_freq") {
  #   nT <- data[1, ]
  #   data <- data[-1, , drop = FALSE]
  # }

  distM = distM[rownames(distM) %in% rownames(data1), colnames(distM) %in%
                  rownames(data1)]
  if (nrow(data1) != nrow(distM))
    stop("The number of species in data should equal to that in distance matrix",
         call. = FALSE)

  if (sum(rownames(data1) %in% rownames(distM)) != nrow(distM))
    stop("Data and distance matrix contain unmatched species",
         call. = FALSE)

  order_sp <- match(rownames(data1), rownames(distM))
  distM <- distM[order_sp, order_sp]



  # if (datatype == "incidence_freq") {
  #   data <- rbind(nT, data)
  # }
  name_sp <- rownames(data1)

  dat <- lapply(1:ncol(data1), function(k) {
    x <- data1[, k]
    names(x) <- name_sp
    x
  })
  if (is.null(colnames(data1))) {
    names(dat) <- paste0("site", 1:length(dat))
  }else {
    names(dat) = colnames(data1)
  }
  if (is.null(threshold)) {
    if (datatype == "abundance") {
      tmp = sapply(dat, function(x) x) %>% apply(., 1,
                                                 sum)
      tmp <- matrix(tmp/sum(tmp), ncol = 1)
    }
    else if (datatype == "incidence_freq") {
      tmp = sapply(dat, function(x) x) %>% apply(., 1,
                                                 sum)
      tmp <- matrix(tmp[-1]/sum(tmp[-1]), ncol = 1)
    }
    dmean <- sum((tmp %*% t(tmp)) * distM)
    dmin <- min(distM[distM > 0])
    threshold <- dmean
  }else if (sum(threshold < 0) > 0 | sum(threshold > 1) > 0) {
    stop("Threshold must be a number between 0 and 1. Use NULL to set it to dmean/2.",
         call. = FALSE)
  }
  if (length(knots) != length(dat))
    knots <- rep(knots, length(dat))
  if (is.null(size)) {
    if (is.null(endpoint)) {
      if (datatype == "abundance") {
        endpoint <- sapply(dat, function(x) 2 * sum(x))
      }
      else if (datatype == "incidence_freq") {
        endpoint <- sapply(dat, function(x) 2 * x[1])
      }
    }
    else {
      if (length(endpoint) != length(dat)) {
        endpoint <- rep(endpoint, length(dat))
      }
    }
    size <- lapply(1:length(dat), function(i) {
      if (datatype == "abundance") {
        ni <- sum(dat[[i]])
      }
      else if (datatype == "incidence_freq") {
        ni <- dat[[i]][1]
      }
      if (endpoint[i] <= ni) {
        mi <- floor(seq(1, endpoint[i], length.out = knots[i]))
      }
      else {
        mi <- floor(c(seq(1, ni, length.out = floor(knots[i]/2)),
                      seq(ni + 1, endpoint[i], length.out = knots[i] -
                            floor(knots[i]/2))))
      }
      unique(mi)
    })
  }else {
    if (class(size) == "numeric" | class(size) == "integer" |
        class(size) == "double") {
      size <- list(size = size)
    }
    if (length(size) != length(dat))
      size = lapply(1:length(dat), function(x) size[[1]])
    size <- lapply(1:length(dat), function(i) {
      if (datatype == "abundance")
        ni <- sum(dat[[i]])
      else ni <- (dat[[i]])[1]
      if (sum(size[[i]] == ni) == 0)
        mi <- sort(c(ni, size[[i]]))
      else mi <- size[[i]]
      unique(mi)
    })
  }

  FUN <- function(e) {
    if (class(dat) == "list") {
      tau <- seq(0, 1, length.out = 30)
      temp1 = iNEXT.3D:::AUCtable_iNextFD(datalist = dat, dij = distM, q = q,
                                          datatype = datatype, tau = tau, nboot = nboot,
                                          conf = conf, m = size)


      temp1$qAUC.LCL[temp1$qAUC.LCL < 0] <- 0
      temp1$SC.LCL[temp1$SC.LCL < 0] <- 0
      temp1$SC.UCL[temp1$SC.UCL > 1] <- 1
      if (datatype == "incidence_freq")
        colnames(temp1)[colnames(temp1) == "m"] = "nt"
      temp2 <- lapply(1:length(dat), function(i) iNEXT.3D:::AUCtable_invFD(datalist = dat[i],
                                                                           dij = distM, q = q, datatype = datatype, level = iNEXT.3D:::Coverage(data = dat[[i]],
                                                                                                                                                datatype = datatype, m = size[[i]]), nboot = nboot,
                                                                           conf = conf, tau =tau)) %>% do.call(rbind,
                                                                                                               .)
      temp2$qAUC.LCL[temp2$qAUC.LCL < 0] <- 0
      if (datatype == "incidence_freq")
        colnames(temp2)[colnames(temp2) == "m"] = "nt"
      temp1$Type = "FD"
      temp2$Type = "FD"
      ans <- list(size_based = temp1, coverage_based = temp2)
      return(ans)
    }
    else {
      return(NULL)
    }
  }
  out <- tryCatch(FUN(e), error = function(e) {
    return()
  })
  index <- iNEXT.3D:::AO3D(data = data1, diversity = "FD", FDdistM = distM, q = c(0, 1, 2), datatype = datatype, nboot = nboot, conf = 0.95)
  index = index[order(index$Assemblage), ]
  LCL <- index$qAUC.LCL[index$Method == "Asymptotic"]
  UCL <- index$qAUC.UCL[index$Method == "Asymptotic"]
  index <- dcast(index, formula = Assemblage + Order.q ~ Method,
                 value.var = "qAUC")
  index <- cbind(index, se = (UCL - index$Asymptotic)/qnorm(1 -
                                                              (1 - conf)/2), LCL, UCL)
  if (nboot > 0)
    index$LCL[index$LCL < index$Empirical & index$Order.q ==
                0] <- index$Empirical[index$LCL < index$Empirical &
                                        index$Order.q == 0]
  index$Order.q <- c("Species richness", "Shannon diversity",
                     "Simpson diversity")
  index[, 3:4] = index[, 4:3]
  colnames(index) <- c("Assemblage", "Functional Diversity",
                       "Functional Observed", "Functional Estimator", "s.e.",
                       "LCL", "UCL")
  info <- iNEXT.3D:::DataInfo3D(data1, diversity = "FD", datatype = datatype,
                                FDdistM = distM, FDtype = "AUC", nT = nT)



  return(list(AUCInfo = info, AUCiNextEst = out, AUCAsyEst = index))
}

AsylinkTD = function (data, diversity = 'TD', q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95)
{
  out = lapply(1:length(data), function(i) {
    x = data[[i]]
    assemblage = names(data)[[i]]
    # please note that nboot has to larger than 0
    res = MakeTable_Proposeprofile(data = x, B = nboot, q, conf = conf)%>%
      mutate(Network = assemblage, method = "Estimated")%>%
      filter(Target == "Diversity")%>%select(-Target, -`s.e.`)%>%
      rename("qD"="Estimate", "qD.LCL"="LCL", "qD.UCL"="UCL", 'Method' = 'method')

    return(res)
  })%>%do.call("rbind",.)

  out
}

ObslinkTD = function (data, diversity = 'TD', q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95)
{
  out = lapply(1:length(data), function(i){
    x = data[[i]]
    assemblage = names(data)[[i]]
    tmp <- c(as.matrix(x))
    ## nboot has to larger than 0
    res = MakeTable_Empericalprofile(data = x, B = nboot, q, conf = conf)%>%
      mutate(Network = assemblage, method = "Empirical")%>%
      filter(Target == "Diversity")%>%select(-Target, -`s.e.`)%>%
      rename("qD"="Emperical", "qD.LCL"="LCL", "qD.UCL"="UCL", 'Method' = 'method')
    return(res)
  })%>%do.call("rbind",.)

  out
}

AsylinkFD = function (data, row.distM = NULL, col.distM = NULL, datatype = "abundance", q = seq(0, 2,
                                                                                                by = 0.25), nboot = 50, conf = 0.95, threshold = NULL, nT = NULL)
{
  if (length(data) == 1) {
    dat = as.matrix(as.vector(t(data[[1]])))
    region_names = if (is.null(names(data)))
      paste0("region_", 1)else
        names(data)
    colnames(dat) = region_names
    rownames(dat) = paste0(rep(colnames(data[[1]]),3),".",rep(rownames(data[[1]]),rep(ncol(data[[1]]),nrow(data[[1]]))))
    data1 = dat
    if(is.null(row.distM)){
      rdd = matrix(1,ncol = nrow(data[[1]]),nrow = nrow(data[[1]]))
      diag(rdd) = 0
      rownames(rdd) = rownames(data[[1]])
      colnames(rdd) = rownames(data[[1]])
      row.distM =  rdd}
    if(is.null(col.distM)){
      cdd = matrix(1,ncol = ncol(data[[1]]),nrow = ncol(data[[1]]))
      diag(cdd) = 0
      rownames(cdd) = colnames(data[[1]])
      colnames(cdd) = colnames(data[[1]])
      col.distM =  cdd}
  }else {
    region_names = if (is.null(names(data)))
      paste0("region_", 1:length(data))else names(data)
    data2 = lapply(data, function(i) {

      int_name = paste0(rep(colnames(i),3),".",rep(rownames(i),rep(ncol(i),nrow(i))))
      i = as.vector(t(i))

      i = data.frame(species = int_name, i)
      return(i)
    })
    data1 = data2[[1]]
    for (i in 2:length(data2)) {
      data1 = data.frame(full_join(data1, data2[[i]],
                                   by = "species"))
    }
    data1[is.na(data1)] = 0
    rownames(data1) = data1$species
    data1 = data1[!colnames(data1) == "species"]
    colnames(data1) = region_names
    row_sp = c()
    col_sp = c()
    for(i in 1:length(data)){
      row_sp = c(row_sp,rownames(data[[i]]))
      col_sp = c(col_sp,colnames(data[[i]]))

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

  }
  row.distM = as.matrix(row.distM)
  col.distM = as.matrix(col.distM)
  # if (datatype == "incidence_raw") {
  #   data = as.incfreq(data, nT = nT)
  #   datatype = "incidence_freq"
  # }
  distM =  1-(1-row.distM[rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM))),rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM)))])*(1-col.distM[rep(1:nrow(col.distM),nrow(row.distM)),rep(1:nrow(col.distM),nrow(row.distM))])
  distM_name = paste0(rep(colnames(col.distM),3),".",rep(rownames(row.distM),rep(ncol(col.distM),nrow(row.distM))))
  colnames(distM) = distM_name
  rownames(distM) = distM_name


  DATATYPE <- c("abundance", "incidence_freq")
  if (is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  if (class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value",
         call. = FALSE)
  if ((conf < 0) | (conf > 1) | (is.numeric(conf) == F))
    stop("conf (confidence level) must be a numerical value between 0 and 1, We use \"conf\" = 0.95 to calculate!",
         call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot) == F))
    stop("nboot must be a nonnegative integer, We use \"nboot\" = 50 to calculate!",
         call. = FALSE)
  if (datatype == "incidence_freq") {
    nT <- data[1, ]
    data <- data[-1, , drop = FALSE]
  }
  distM = distM[rownames(distM) %in% rownames(data1), colnames(distM) %in%
                  rownames(data1)]
  if (nrow(data1) != nrow(distM))
    stop("The number of species in data should equal to that in distance matrix",
         call. = FALSE)


  order_sp <- match(rownames(data1), rownames(distM))
  distM <- distM[order_sp, order_sp]

  if (datatype == "incidence_freq") {
    data <- rbind(nT, data)
  }
  name_sp <- rownames(data1)
  dat <- lapply(1:ncol(data1), function(k) {
    x <- data1[, k]
    names(x) <- name_sp
    x
  })
  if (is.null(colnames(data1))) {
    names(dat) <- paste0("site", 1:length(dat))
  }else {
    names(dat) = colnames(data1)
  }
  if (is.null(threshold)) {
    if (datatype == "abundance") {
      tmp = sapply(dat, function(x) x) %>% apply(., 1,
                                                 sum)
      tmp <- matrix(tmp/sum(tmp), ncol = 1)
    }
    else if (datatype == "incidence_freq") {
      tmp = sapply(dat, function(x) x) %>% apply(., 1,
                                                 sum)
      tmp <- matrix(tmp[-1]/sum(tmp[-1]), ncol = 1)
    }
    dmean <- sum((tmp %*% t(tmp)) * distM)
    dmin <- min(distM[distM > 0])
    threshold <- dmean
  }else if (sum(threshold < 0) > 0 | sum(threshold > 1) > 0) {
    stop("Threshold must be a number between 0 and 1. Use NULL to set it to dmean/2.",
         call. = FALSE)
  }
  out <- iNEXT.3D:::FDtable_est(datalist = dat, dij = distM, q = q, datatype = datatype,
                                nboot = nboot, conf = conf, tau = threshold)
  out
}

AsylinkAUC = function (data, row.distM = NULL, col.distM = NULL, datatype = "abundance", q = seq(0, 2,by = 0.25), nboot = 20, conf = 0.95, tau = NULL, nT = NULL)
{
  if (length(data) == 1) {
    dat = as.matrix(as.vector(t(data[[1]])))
    region_names = if (is.null(names(data)))
      paste0("region_", 1)else
        names(data)
    colnames(dat) = region_names
    rownames(dat) = paste0(rep(colnames(data[[1]]),3),".",rep(rownames(data[[1]]),rep(ncol(data[[1]]),nrow(data[[1]]))))
    data1 = dat
    if(is.null(row.distM)){
      rdd = matrix(1,ncol = nrow(data[[1]]),nrow = nrow(data[[1]]))
      diag(rdd) = 0
      rownames(rdd) = rownames(data[[1]])
      colnames(rdd) = rownames(data[[1]])
      row.distM =  rdd}
    if(is.null(col.distM)){
      cdd = matrix(1,ncol = ncol(data[[1]]),nrow = ncol(data[[1]]))
      diag(cdd) = 0
      rownames(cdd) = colnames(data[[1]])
      colnames(cdd) = colnames(data[[1]])
      col.distM =  cdd}
  }else {
    region_names = if (is.null(names(data)))
      paste0("region_", 1:length(data))else names(data)
    data2 = lapply(data, function(i) {

      int_name = paste0(rep(colnames(i),3),".",rep(rownames(i),rep(ncol(i),nrow(i))))
      i = as.vector(t(i))

      i = data.frame(species = int_name, i)
      return(i)
    })
    data1 = data2[[1]]
    for (i in 2:length(data2)) {
      data1 = data.frame(full_join(data1, data2[[i]],
                                   by = "species"))
    }
    data1[is.na(data1)] = 0
    rownames(data1) = data1$species
    data1 = data1[!colnames(data1) == "species"]
    colnames(data1) = region_names
    row_sp = c()
    col_sp = c()
    for(i in 1:length(data)){
      row_sp = c(row_sp,rownames(data[[i]]))
      col_sp = c(col_sp,colnames(data[[i]]))

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

  }
  row.distM = as.matrix(row.distM)
  col.distM = as.matrix(col.distM)
  # if (datatype == "incidence_raw") {
  #   data = as.incfreq(data, nT = nT)
  #   datatype = "incidence_freq"
  # }
  distM =  1-(1-row.distM[rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM))),rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM)))])*(1-col.distM[rep(1:nrow(col.distM),nrow(row.distM)),rep(1:nrow(col.distM),nrow(row.distM))])
  distM_name = paste0(rep(colnames(col.distM),3),".",rep(rownames(row.distM),rep(ncol(col.distM),nrow(row.distM))))
  colnames(distM) = distM_name
  rownames(distM) = distM_name

  DATATYPE <- c("abundance", "incidence_freq")
  if (is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  if (class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value",
         call. = FALSE)
  if ((conf < 0) | (conf > 1) | (is.numeric(conf) == F))
    stop("conf (confidence level) must be a numerical value between 0 and 1, We use \"conf\" = 0.95 to calculate!",
         call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot) == F))
    stop("nboot must be a nonnegative integer, We use \"nboot\" = 50 to calculate!",
         call. = FALSE)
  if (datatype == "incidence_freq") {
    nT <- data[1, ]
    data <- data[-1, , drop = FALSE]
  }
  distM = distM[rownames(distM) %in% rownames(data1), colnames(distM) %in%
                  rownames(data1)]
  if (nrow(data1) != nrow(distM))
    stop("The number of species in data should equal to that in distance matrix",
         call. = FALSE)

  order_sp <- match(rownames(data1), rownames(distM))
  distM <- distM[order_sp, order_sp]

  if (datatype == "incidence_freq") {
    data <- rbind(nT, data)
  }
  name_sp <- rownames(data1)
  dat <- lapply(1:ncol(data1), function(k) {
    x <- data1[, k]
    names(x) <- name_sp
    x
  })
  if (is.null(colnames(data1))) {
    names(dat) <- paste0("site", 1:length(dat))
  }else {
    names(dat) = colnames(data1)
  }
  out <- iNEXT.3D:::AUCtable_est(datalist = dat, dij = distM, q = q,
                                 datatype = datatype, nboot = nboot, conf = conf, tau = NULL)
  out
}


ObslinkFD = function (data, row.distM = NULL, col.distM = NULL, datatype = "abundance", q = seq(0, 2,
                                                                                                by = 0.25), nboot = 50, conf = 0.95, threshold = NULL, nT = NULL)
{
  if (length(data) == 1) {
    dat = as.matrix(as.vector(t(data[[1]])))
    region_names = if (is.null(names(data)))
      paste0("region_", 1)else
        names(data)
    colnames(dat) = region_names
    rownames(dat) = paste0(rep(colnames(data[[1]]),3),".",rep(rownames(data[[1]]),rep(ncol(data[[1]]),nrow(data[[1]]))))
    data1 = dat
    if(is.null(row.distM)){
      rdd = matrix(1,ncol = nrow(data[[1]]),nrow = nrow(data[[1]]))
      diag(rdd) = 0
      rownames(rdd) = rownames(data[[1]])
      colnames(rdd) = rownames(data[[1]])
      row.distM =  rdd}
    if(is.null(col.distM)){
      cdd = matrix(1,ncol = ncol(data[[1]]),nrow = ncol(data[[1]]))
      diag(cdd) = 0
      rownames(cdd) = colnames(data[[1]])
      colnames(cdd) = colnames(data[[1]])
      col.distM =  cdd}
  }else {
    region_names = if (is.null(names(data)))
      paste0("region_", 1:length(data))else names(data)
    data2 = lapply(data, function(i) {

      int_name = paste0(rep(colnames(i),3),".",rep(rownames(i),rep(ncol(i),nrow(i))))
      i = as.vector(t(i))

      i = data.frame(species = int_name, i)
      return(i)
    })
    data1 = data2[[1]]
    for (i in 2:length(data2)) {
      data1 = data.frame(full_join(data1, data2[[i]],
                                   by = "species"))
    }
    data1[is.na(data1)] = 0
    rownames(data1) = data1$species
    data1 = data1[!colnames(data1) == "species"]
    colnames(data1) = region_names
    row_sp = c()
    col_sp = c()
    for(i in 1:length(data)){
      row_sp = c(row_sp,rownames(data[[i]]))
      col_sp = c(col_sp,colnames(data[[i]]))

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

  }
  row.distM = as.matrix(row.distM)
  col.distM = as.matrix(col.distM)
  # if (datatype == "incidence_raw") {
  #   data = as.incfreq(data, nT = nT)
  #   datatype = "incidence_freq"
  # }
  distM =  1-(1-row.distM[rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM))),rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM)))])*(1-col.distM[rep(1:nrow(col.distM),nrow(row.distM)),rep(1:nrow(col.distM),nrow(row.distM))])
  distM_name = paste0(rep(colnames(col.distM),3),".",rep(rownames(row.distM),rep(ncol(col.distM),nrow(row.distM))))
  colnames(distM) = distM_name
  rownames(distM) = distM_name

  DATATYPE <- c("abundance", "incidence_freq")
  if (is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  if (class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value",
         call. = FALSE)
  if ((conf < 0) | (conf > 1) | (is.numeric(conf) == F))
    stop("conf (confidence level) must be a numerical value between 0 and 1, We use \"conf\" = 0.95 to calculate!",
         call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot) == F))
    stop("nboot must be a nonnegative integer, We use \"nboot\" = 50 to calculate!",
         call. = FALSE)
  if (datatype == "incidence_freq") {
    nT <- data[1, ]
    data <- data[-1, , drop = FALSE]
  }
  distM = distM[rownames(distM) %in% rownames(data1), colnames(distM) %in%
                  rownames(data1)]
  if (nrow(data1) != nrow(distM))
    stop("The number of species in data should equal to that in distance matrix",
         call. = FALSE)


  order_sp <- match(rownames(data1), rownames(distM))
  distM <- distM[order_sp, order_sp]

  if (datatype == "incidence_freq") {
    data <- rbind(nT, data)
  }
  name_sp <- rownames(data1)
  dat <- lapply(1:ncol(data1), function(k) {
    x <- data1[, k]
    names(x) <- name_sp
    x
  })
  if (is.null(colnames(data1))) {
    names(dat) <- paste0("site", 1:length(dat))
  }else {
    names(dat) = colnames(data1)
  }
  if (is.null(threshold)) {
    if (datatype == "abundance") {
      tmp = sapply(dat, function(x) x) %>% apply(., 1,
                                                 sum)
      tmp <- matrix(tmp/sum(tmp), ncol = 1)
    }
    else if (datatype == "incidence_freq") {
      tmp = sapply(dat, function(x) x) %>% apply(., 1,
                                                 sum)
      tmp <- matrix(tmp[-1]/sum(tmp[-1]), ncol = 1)
    }
    dmean <- sum((tmp %*% t(tmp)) * distM)
    dmin <- min(distM[distM > 0])
    threshold <- dmean
  }else if (sum(threshold < 0) > 0 | sum(threshold > 1) > 0) {
    stop("Threshold must be a number between 0 and 1. Use NULL to set it to dmean/2.",
         call. = FALSE)
  }
  out <- iNEXT.3D:::FDtable_mle(datalist = dat, dij = distM, q = q, datatype = datatype,
                                nboot = nboot, conf = conf, tau = threshold)
  out
}

ObslinkAUC = function (data, row.distM = NULL, col.distM = NULL, datatype = "abundance", q = seq(0, 2,
                                                                                                 by = 0.25), nboot = 50, conf = 0.95, tau = NULL, nT = NULL)
{
  if (length(data) == 1) {
    dat = as.matrix(as.vector(t(data[[1]])))
    region_names = if (is.null(names(data)))
      paste0("region_", 1)else
        names(data)
    colnames(dat) = region_names
    rownames(dat) = paste0(rep(colnames(data[[1]]),3),".",rep(rownames(data[[1]]),rep(ncol(data[[1]]),nrow(data[[1]]))))
    data1 = dat
    if(is.null(row.distM)){
      rdd = matrix(1,ncol = nrow(data[[1]]),nrow = nrow(data[[1]]))
      diag(rdd) = 0
      rownames(rdd) = rownames(data[[1]])
      colnames(rdd) = rownames(data[[1]])
      row.distM =  rdd}
    if(is.null(col.distM)){
      cdd = matrix(1,ncol = ncol(data[[1]]),nrow = ncol(data[[1]]))
      diag(cdd) = 0
      rownames(cdd) = colnames(data[[1]])
      colnames(cdd) = colnames(data[[1]])
      col.distM =  cdd}
  }else {
    region_names = if (is.null(names(data)))
      paste0("region_", 1:length(data))else names(data)
    data2 = lapply(data, function(i) {

      int_name = paste0(rep(colnames(i),3),".",rep(rownames(i),rep(ncol(i),nrow(i))))
      i = as.vector(t(i))

      i = data.frame(species = int_name, i)
      return(i)
    })
    data1 = data2[[1]]
    for (i in 2:length(data2)) {
      data1 = data.frame(full_join(data1, data2[[i]],
                                   by = "species"))
    }
    data1[is.na(data1)] = 0
    rownames(data1) = data1$species
    data1 = data1[!colnames(data1) == "species"]
    colnames(data1) = region_names
    row_sp = c()
    col_sp = c()
    for(i in 1:length(data)){
      row_sp = c(row_sp,rownames(data[[i]]))
      col_sp = c(col_sp,colnames(data[[i]]))

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

  }
  row.distM = as.matrix(row.distM)
  col.distM = as.matrix(col.distM)
  # if (datatype == "incidence_raw") {
  #   data = as.incfreq(data, nT = nT)
  #   datatype = "incidence_freq"
  # }
  distM =  1-(1-row.distM[rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM))),rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM)))])*(1-col.distM[rep(1:nrow(col.distM),nrow(row.distM)),rep(1:nrow(col.distM),nrow(row.distM))])
  distM_name = paste0(rep(colnames(col.distM),3),".",rep(rownames(row.distM),rep(ncol(col.distM),nrow(row.distM))))
  colnames(distM) = distM_name
  rownames(distM) = distM_name

  DATATYPE <- c("abundance", "incidence_freq")
  if (is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  if (class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value",
         call. = FALSE)
  if ((conf < 0) | (conf > 1) | (is.numeric(conf) == F))
    stop("conf (confidence level) must be a numerical value between 0 and 1, We use \"conf\" = 0.95 to calculate!",
         call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot) == F))
    stop("nboot must be a nonnegative integer, We use \"nboot\" = 50 to calculate!",
         call. = FALSE)
  if (datatype == "incidence_freq") {
    nT <- data[1, ]
    data <- data[-1, , drop = FALSE]
  }
  distM = distM[rownames(distM) %in% rownames(data1), colnames(distM) %in%
                  rownames(data1)]
  if (nrow(data1) != nrow(distM))
    stop("The number of species in data should equal to that in distance matrix",
         call. = FALSE)

  order_sp <- match(rownames(data1), rownames(distM))
  distM <- distM[order_sp, order_sp]

  if (datatype == "incidence_freq") {
    data <- rbind(nT, data)
  }
  name_sp <- rownames(data1)
  dat <- lapply(1:ncol(data1), function(k) {
    x <- data1[, k]
    names(x) <- name_sp
    x
  })
  if (is.null(colnames(data1))) {
    names(dat) <- paste0("site", 1:length(dat))
  }else {
    names(dat) = colnames(data1)
  }
  out <- iNEXT.3D:::AUCtable_mle(datalist = dat, dij = distM, q = q,
                                 datatype = datatype, nboot = nboot, conf = conf, tau = tau)
  out
}

estimatelinkFD = function (data, row.distM = NULL, col.distM = NULL, datatype = "abundance", q = c(0, 1, 2),
                           base = "coverage", threshold = NULL, level = NULL, nboot = 50,
                           conf = 0.95, nT = NULL)
{
  if (length(data) == 1) {
    dat = as.matrix(as.vector(t(data[[1]])))
    region_names = if (is.null(names(data)))
      paste0("region_", 1)else
        names(data)
    colnames(dat) = region_names
    rownames(dat) = paste0(rep(colnames(data[[1]]),3),".",rep(rownames(data[[1]]),rep(ncol(data[[1]]),nrow(data[[1]]))))
    data1 = dat
    if(is.null(row.distM)){
      rdd = matrix(1,ncol = nrow(data[[1]]),nrow = nrow(data[[1]]))
      diag(rdd) = 0
      rownames(rdd) = rownames(data[[1]])
      colnames(rdd) = rownames(data[[1]])
      row.distM =  rdd}
    if(is.null(col.distM)){
      cdd = matrix(1,ncol = ncol(data[[1]]),nrow = ncol(data[[1]]))
      diag(cdd) = 0
      rownames(cdd) = colnames(data[[1]])
      colnames(cdd) = colnames(data[[1]])
      col.distM =  cdd}
  }else {
    region_names = if (is.null(names(data)))
      paste0("region_", 1:length(data))else names(data)
    data2 = lapply(data, function(i) {

      int_name = paste0(rep(colnames(i),3),".",rep(rownames(i),rep(ncol(i),nrow(i))))
      i = as.vector(t(i))

      i = data.frame(species = int_name, i)
      return(i)
    })
    data1 = data2[[1]]
    for (i in 2:length(data2)) {
      data1 = data.frame(full_join(data1, data2[[i]],
                                   by = "species"))
    }
    data1[is.na(data1)] = 0
    rownames(data1) = data1$species
    data1 = data1[!colnames(data1) == "species"]
    colnames(data1) = region_names
    row_sp = c()
    col_sp = c()
    for(i in 1:length(data)){
      row_sp = c(row_sp,rownames(data[[i]]))
      col_sp = c(col_sp,colnames(data[[i]]))

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

  }
  row.distM = as.matrix(row.distM)
  col.distM = as.matrix(col.distM)
  # if (datatype == "incidence_raw") {
  #   data = as.incfreq(data, nT = nT)
  #   datatype = "incidence_freq"
  # }
  distM =  1-(1-row.distM[rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM))),rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM)))])*(1-col.distM[rep(1:nrow(col.distM),nrow(row.distM)),rep(1:nrow(col.distM),nrow(row.distM))])
  distM_name = paste0(rep(colnames(col.distM),3),".",rep(rownames(row.distM),rep(ncol(col.distM),nrow(row.distM))))
  colnames(distM) = distM_name
  rownames(distM) = distM_name


  DATATYPE <- c("abundance", "incidence_freq")
  if (is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  if (class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value",
         call. = FALSE)
  if ((conf < 0) | (conf > 1) | (is.numeric(conf) == F))
    stop("conf (confidence level) must be a numerical value between 0 and 1, We use \"conf\" = 0.95 to calculate!",
         call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot) == F))
    stop("nboot must be a nonnegative integer, We use \"nboot\" = 50 to calculate!",
         call. = FALSE)
  # if (datatype == "incidence_freq") {
  #   nT <- data[1, ]
  #   data <- data[-1, , drop = FALSE]
  # }
  distM = distM[rownames(distM) %in% rownames(data1), colnames(distM) %in%
                  rownames(data1)]
  if (nrow(data1) != nrow(distM))
    stop("The number of species in data should equal to that in distance matrix",
         call. = FALSE)


  order_sp <- match(rownames(data1), rownames(distM))
  distM <- distM[order_sp, order_sp]

  if (datatype == "incidence_freq") {
    data <- rbind(nT, data)
  }
  name_sp <- rownames(data1)
  dat <- lapply(1:ncol(data1), function(k) {
    x <- data1[, k]
    names(x) <- name_sp
    x
  })
  if (is.null(colnames(data1))) {
    names(dat) <- paste0("site", 1:length(dat))
  }else {
    names(dat) = colnames(data1)
  }
  if (is.null(threshold)) {
    if (datatype == "abundance") {
      tmp = sapply(dat, function(x) x) %>% apply(., 1,
                                                 sum)
      tmp <- matrix(tmp/sum(tmp), ncol = 1)
    }
    else if (datatype == "incidence_freq") {
      tmp = sapply(dat, function(x) x) %>% apply(., 1,
                                                 sum)
      tmp <- matrix(tmp[-1]/sum(tmp[-1]), ncol = 1)
    }
    dmean <- sum((tmp %*% t(tmp)) * distM)
    dmin <- min(distM[distM > 0])
    threshold <- dmean
  }else if (sum(threshold < 0) > 0 | sum(threshold > 1) > 0) {
    stop("Threshold must be a number between 0 and 1. Use NULL to set it to dmean/2.",
         call. = FALSE)
  }
  if (is.null(level) & base == "size") {
    if (datatype == "abundance") {
      level <- sapply(dat, function(x) 2 * sum(x))
    }
    else if (datatype == "incidence_freq") {
      level <- sapply(dat, function(x) 2 * x[1])
    }
    level <- lapply(dat, function(i) min(level))
  }
  else if (is.null(level) & base == "coverage") {
    if (datatype == "abundance") {
      level <- sapply(dat, function(x) {
        ni <- sum(x)
        iNEXT.3D:::Coverage(data = x, datatype = datatype, m = 2 *
                              ni)
      })
    }
    else if (datatype == "incidence_freq") {
      level <- sapply(dat, function(x) {
        ni <- x[1]
        iNEXT.3D:::Coverage(data = x, datatype = datatype, m = 2 *
                              ni)
      })
    }
    level <- min(level)
  }
  if (base == "size") {
    out = iNEXT.3D:::iNextFD(datalist = dat, dij = distM, q = q, datatype = datatype,
                             tau = threshold, nboot = nboot, conf = conf, m = level) %>%
      select(-c("SC.s.e.", "SC.LCL", "SC.UCL"))
    out$qFD.LCL[out$qFD.LCL < 0] <- 0
    if (datatype == "incidence_freq")
      colnames(out)[colnames(out) == "m"] = "nt"
  }
  else if (base == "coverage") {
    out <- iNEXT.3D:::invChatFD(datalist = dat, dij = distM, q = q,
                                datatype = datatype, level = level, nboot = nboot,
                                conf = conf, tau = threshold)
    out$qFD.LCL[out$qFD.LCL < 0] <- 0
    if (datatype == "incidence_freq")
      colnames(out)[colnames(out) == "m"] = "nt"
  }
  return(out)
}

estimatelinkAUC = function (data, row.distM = NULL, col.distM = NULL, datatype = "abundance", q = c(0, 1, 2),
                            base = "coverage", level = NULL, nboot = 50, conf = 0.95,
                            tau = NULL, nT = NULL)
{

  if (length(data) == 1) {
    dat = as.matrix(as.vector(t(data[[1]])))
    region_names = if (is.null(names(data)))
      paste0("region_", 1)else
        names(data)
    colnames(dat) = region_names
    rownames(dat) = paste0(rep(colnames(data[[1]]),3),".",rep(rownames(data[[1]]),rep(ncol(data[[1]]),nrow(data[[1]]))))
    data1 = dat
    if(is.null(row.distM)){
      rdd = matrix(1,ncol = nrow(data[[1]]),nrow = nrow(data[[1]]))
      diag(rdd) = 0
      rownames(rdd) = rownames(data[[1]])
      colnames(rdd) = rownames(data[[1]])
      row.distM =  rdd}
    if(is.null(col.distM)){
      cdd = matrix(1,ncol = ncol(data[[1]]),nrow = ncol(data[[1]]))
      diag(cdd) = 0
      rownames(cdd) = colnames(data[[1]])
      colnames(cdd) = colnames(data[[1]])
      col.distM =  cdd}
  }else {
    region_names = if (is.null(names(data)))
      paste0("region_", 1:length(data))else names(data)
    data2 = lapply(data, function(i) {

      int_name = paste0(rep(colnames(i),3),".",rep(rownames(i),rep(ncol(i),nrow(i))))
      i = as.vector(t(i))

      i = data.frame(species = int_name, i)
      return(i)
    })
    data1 = data2[[1]]
    for (i in 2:length(data2)) {
      data1 = data.frame(full_join(data1, data2[[i]],
                                   by = "species"))
    }
    data1[is.na(data1)] = 0
    rownames(data1) = data1$species
    data1 = data1[!colnames(data1) == "species"]
    colnames(data1) = region_names
    row_sp = c()
    col_sp = c()
    for(i in 1:length(data)){
      row_sp = c(row_sp,rownames(data[[i]]))
      col_sp = c(col_sp,colnames(data[[i]]))

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

  }
  row.distM = as.matrix(row.distM)
  col.distM = as.matrix(col.distM)
  # if (datatype == "incidence_raw") {
  #   data = as.incfreq(data, nT = nT)
  #   datatype = "incidence_freq"
  # }
  distM =  1-(1-row.distM[rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM))),rep(1:nrow(row.distM),rep(nrow(col.distM),nrow(row.distM)))])*(1-col.distM[rep(1:nrow(col.distM),nrow(row.distM)),rep(1:nrow(col.distM),nrow(row.distM))])
  distM_name = paste0(rep(colnames(col.distM),3),".",rep(rownames(row.distM),rep(ncol(col.distM),nrow(row.distM))))
  colnames(distM) = distM_name
  rownames(distM) = distM_name

  DATATYPE <- c("abundance", "incidence_freq")
  if (is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  if (class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value",
         call. = FALSE)
  if ((conf < 0) | (conf > 1) | (is.numeric(conf) == F))
    stop("conf (confidence level) must be a numerical value between 0 and 1, We use \"conf\" = 0.95 to calculate!",
         call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot) == F))
    stop("nboot must be a nonnegative integer, We use \"nboot\" = 50 to calculate!",
         call. = FALSE)
  # if (datatype == "incidence_freq") {
  #   nT <- data[1, ]
  #   data <- data[-1, , drop = FALSE]
  # }
  distM = distM[rownames(distM) %in% rownames(data1), colnames(distM) %in%
                  rownames(data1)]
  if (nrow(data1) != nrow(distM))
    stop("The number of species in data should equal to that in distance matrix",
         call. = FALSE)

  order_sp <- match(rownames(data1), rownames(distM))
  distM <- distM[order_sp, order_sp]

  if (datatype == "incidence_freq") {
    data <- rbind(nT, data)
  }
  name_sp <- rownames(data1)
  dat <- lapply(1:ncol(data1), function(k) {
    x <- data1[, k]
    names(x) <- name_sp
    x
  })
  if (is.null(colnames(data1))) {
    names(dat) <- paste0("site", 1:length(dat))
  }else {
    names(dat) = colnames(data1)
  }
  if (is.null(level) & base == "size") {
    if (datatype == "abundance") {
      level <- sapply(dat, function(x) 2 * sum(x))
    }
    else if (datatype == "incidence_freq") {
      level <- sapply(dat, function(x) 2 * x[1])
    }
    level <- lapply(dat, function(i) min(level))
  }
  else if (is.null(level) & base == "coverage") {
    if (datatype == "abundance") {
      level <- sapply(dat, function(x) {
        ni <- sum(x)
        iNEXT.3D:::Coverage(data = x, datatype = datatype, m = 2 *
                              ni)
      })
    }
    else if (datatype == "incidence_freq") {
      level <- sapply(dat, function(x) {
        ni <- x[1]
        iNEXT.3D:::Coverage(data = x, datatype = datatype, m = 2 *
                              ni)
      })
    }
    level <- min(level)
  }
  if (base == "size") {
    out = iNEXT.3D:::AUCtable_iNextFD(datalist = dat, dij = distM,
                                      q = q, datatype = datatype, tau = tau, nboot = nboot,
                                      conf = conf, m = level) %>% select(-c("SC.s.e.",
                                                                            "SC.LCL", "SC.UCL"))
    out$qAUC.LCL[out$qAUC.LCL < 0] <- 0
  }
  else if (base == "coverage") {

    out <- iNEXT.3D:::AUCtable_invFD(datalist = dat, dij = distM, q = q,
                                     datatype = datatype, level = level, nboot = nboot,
                                     conf = conf, tau = tau)
    if (datatype == "incidence_freq")
      colnames(out)[colnames(out) == "m"] = "nt"
  }
  return(out)
}
##
Diversity_profile <- function(x,q){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  r <- 1:(n-1)
  Sub <- function(q){
    if(q==0){
      sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }
    else if(q==1){
      A <- sum(x/n*(digamma(n)-digamma(x)))
      B <- ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(-log(p1)-sum((1-p1)^r/r)))
      exp(A+B)
    }else if(abs(q-round(q))==0){
      A <- sum(exp(lchoose(x,q)-lchoose(n,q)))
      #ifelse(A==0,NA,A^(1/(1-q)))
      A^(1/(1-q))
    }else {
      sort.data = sort(unique(x))
      tab = table(x)
      term = sapply(sort.data,function(z){
        k=0:(n-z)
        sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
      })
      r <- 0:(n-1)
      A = sum(tab*term)
      B = ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(p1^(q-1)-sum(choose(q-1,r)*(p1-1)^r)))
      (A+B)^(1/(1-q))
    }
  }
  sapply(q, Sub)
}
Diversity_profile_MLE <- function(x,q){
  p <- x[x>0]/sum(x)
  Sub <- function(q){
    if(q==0) sum(p>0)
    else if(q==1) exp(-sum(p*log(p)))
    else exp(1/(1-q)*log(sum(p^q)))
  }
  sapply(q, Sub)
}

Diversity_Tsallis <- function(x,q){
  qD = Diversity_profile(x, q)
  ans = rep(0,length(qD))
  ans[which(q==1)] <- log(qD[which(q==1)])
  q1 = q[q!=1]
  ans[which(q!=1)] <- (qD[which(q!=1)]^(1-q1)-1)/(1-q1)
  return(ans)
}
Diversity_Tsallis_MLE <- function(x,q){
  qD = Diversity_profile_MLE(x,q)
  ans = rep(0,length(qD))
  ans[which(q==1)] <- log(qD[which(q==1)])
  q1 = q[q!=1]
  ans[which(q!=1)] <- (qD[which(q!=1)]^(1-q1)-1)/(1-q1)
  return(ans)
}
bootstrap_forq = function(data,B,q,conf,FUNNAME){
  data <- data[data!=0]
  n <- sum(data)
  f1 = sum(data==1); f2 = sum(data==2)
  f0 = ceiling(ifelse( f2>0, (n-1)*f1^2/n/2/f2, (n-1)*f1*(f1-1)/2/n ))
  C_hat = ifelse(f2>0, 1-f1*(n-1)*f1/n/((n-1)*f1+2*f2), 1-f1*(n-1)*(f1-1)/n/((n-1)*(f1-1)+2))
  lamda_hat = (1-C_hat)/sum((data/n)*(1-data/n)^n)
  pi_hat = (data/n)*(1-lamda_hat*((1-data/n)^n))
  p_hat = c( pi_hat, rep( (1-C_hat)/f0, f0 ))
  random = rmultinom( B, n, p_hat )
  #Bt_estimate <- sapply(c(1:B),function(i) FUNNAME(random[,i],q))
  Bt_estimate <- apply(random,MARGIN = 2,function(i) FUNNAME(i,q))
  estimate <- FUNNAME(data,q)
  #Interval_mean = apply( Bt_estimate, 1, mean)
  Interval_mean = rowMeans(Bt_estimate)
  Interval_sd = apply(Bt_estimate, 1, sd)
  Interval_quantileL = apply( Bt_estimate, 1, quantile, p=(1-conf)/2)
  Interval_quantileU = apply( Bt_estimate, 1, quantile, p=1-(1-conf)/2)
  Upper_bound = estimate+Interval_quantileU-Interval_mean
  Lower_bound = estimate+Interval_quantileL-Interval_mean
  result <- cbind("estimate"=estimate,"sd"=Interval_sd,"LCL"=Lower_bound,"UCL"=Upper_bound)
  result
}

MakeTable_Proposeprofile = function(data, B, q, conf){
  Diversity = bootstrap_forq(data, B, q, conf, Diversity_profile)
  Entropy = bootstrap_forq(data, B, q, conf, Diversity_Tsallis)
  # tmp <- Diversity_Tsallis(Diversity[,1],q)
  # Entropy = cbind("estimate"=tmp,"sd"=Diversity[,2],"LCL"=tmp-Diversity[,2],"UCL"=tmp+Diversity[,2])
  output = rbind(data.frame("Order.q" = q,"Target"="Diversity","Estimate"=Diversity[,1],"s.e."=Diversity[,2],"LCL"=Diversity[,3],"UCL"=Diversity[,4]),
                 data.frame("Order.q" = q,"Target"="Entropy","Estimate"=Entropy[,1],"s.e."=Entropy[,2],"LCL"=Entropy[,3],"UCL"=Entropy[,4]))
  output[,c(1,3,4,5,6)] = round(output[,c(1,3,4,5,6)],9)

  return(output)
}
MakeTable_Empericalprofile = function(data, B, q, conf){
  Diversity = bootstrap_forq( data, B,q,conf,Diversity_profile_MLE)
  Entropy = bootstrap_forq( data, B,q,conf,Diversity_Tsallis_MLE)
  # tmp <- Diversity_Tsallis(Diversity[,1],q)
  # Entropy = cbind("estimate"=tmp,"sd"=Diversity[,2],"LCL"=tmp-Diversity[,2],"UCL"=tmp+Diversity[,2])
  output = rbind(data.frame("Order.q" = q,"Target"="Diversity","Emperical"=Diversity[,1],"s.e."=Diversity[,2],"LCL"=Diversity[,3],"UCL"=Diversity[,4]),
                 data.frame("Order.q" = q,"Target"="Entropy","Emperical"=Entropy[,1],"s.e."=Entropy[,2],"LCL"=Entropy[,3],"UCL"=Entropy[,4]))
  output[,c(1,3,4,5,6)] = round(output[,c(1,3,4,5,6)],9)

  return(output)
}
coverage_to_size = function (x, C, datatype = "abundance")
{
  if (datatype == "abundance") {
    n <- sum(x)
    refC <- iNEXT.3D:::Coverage(x, "abundance", n)
    f <- function(m, C) abs(iNEXT.3D:::Coverage(x, "abundance",
                                                m) - C)
    if (refC == C) {
      mm = n
    }
    else if (refC > C) {
      opt <- optimize(f, C = C, lower = 0, upper = sum(x))
      mm <- opt$minimum
    }
    else if (refC < C) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }
      if (f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) +
                                   2)
      }
      if (f1 == 1 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 > 0) {
        A <- 1
      }
      mm <- (log(n/f1) + log(1 - C))/log(A) - 1
      if (is.nan(mm) == TRUE)
        mm = Inf
      mm <- n + mm
    }
  }
  else {
    m <- NULL
    n <- max(x)
    refC <- iNEXT.3D:::Coverage(x, "incidence_freq", n)
    f <- function(m, C) abs(iNEXT.3D:::Coverage(x, "incidence_freq",
                                                m) - C)
    if (refC == C) {
      mm = n
    }
    else if (refC > C) {
      opt <- optimize(f, C = C, lower = 0, upper = max(x))
      mm <- opt$minimum
    }
    else if (refC < C) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      U <- sum(x) - max(x)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }
      if (f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) +
                                   2)
      }
      if (f1 == 1 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 > 0) {
        A <- 1
      }
      mm <- (log(U/f1) + log(1 - C))/log(A) - 1
      if (is.nan(mm) == TRUE)
        mm = Inf
      mm <- n + mm
    }
  }
  return(mm)
}

bootstrap_population_multiple_assemblage = function(data, data_gamma, datatype){

  if (datatype == 'abundance'){

    S_obs = sum(data_gamma > 0)
    n = sum(data_gamma)
    f1 = sum(data_gamma == 1)
    f2 = sum(data_gamma == 2)
    f0_hat = ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2, (n - 1)/n * f1^2/2/f2) %>% ceiling()

    output = apply(data, 2, function(x){

      p_i_hat = iNEXT.3D:::EstiBootComm.Ind(Spec = x)

      if(length(p_i_hat) != length(x)){

        p_i_hat_unobs = p_i_hat[(length(x)+1):length(p_i_hat)]
        p_i_hat_obs = p_i_hat[1:length(x)]
        p_i_hat = c(p_i_hat_obs, rep(0, f0_hat))
        candidate = which(p_i_hat==0)

        chosen = sample(x = candidate, size = min(length(p_i_hat_unobs), length(candidate)), replace = F)
        p_i_hat[chosen] = (1-sum(p_i_hat))/length(chosen)

        p_i_hat

      } else {

        p_i_hat = c(p_i_hat, rep(0, f0_hat))
        p_i_hat

      }
    })

  }

  if (datatype == 'incidence'){

    S_obs = sum(data_gamma > 0)
    t = data_gamma[1]
    Q1 = sum(data_gamma == 1)
    Q2 = sum(data_gamma == 2)
    Q0_hat = if ( Q2 == 0 ){( (t-1)/t ) * ( Q1*(Q1-1)/2 )} else {( (t-1)/t ) * ( (Q1^2) / (2*Q2) )} %>% ceiling

    output = apply(data, 2, function(x){

      pi_i_hat = iNEXT.3D:::EstiBootComm.Sam(Spec = x)

      if(length(pi_i_hat) != (length(x) - 1)){

        pi_i_hat_unobs = pi_i_hat[length(x):length(pi_i_hat)]
        pi_i_hat_obs = pi_i_hat[1:(length(x)-1)]
        pi_i_hat = c(pi_i_hat_obs, rep(0, Q0_hat))
        candidate = which(pi_i_hat == 0)
        chosen = sample(x = candidate, size = min(length(pi_i_hat_unobs), length(candidate)), replace = F)
        pi_i_hat[chosen] = pi_i_hat_unobs

        pi_i_hat

      } else {

        pi_i_hat = c(pi_i_hat, rep(0, Q0_hat))
        pi_i_hat

      }
    })

  }

  return(output)

}

Bootstrap_distance_matrix = function(data, distance_matrix, f0.hat, datatype){

  if (datatype == "incidence_freq") {
    n = data[1]
    X = data[-1]
    u = sum(data)
  } else if (datatype == "abundance") {
    n = sum(data)
    X = data
  }

  # n = sum(data)
  distance = as.matrix(distance_matrix)
  dij = distance
  # X = data

  F.1 <- sum(dij[, X==1]) ; F.2 <- sum(dij[, X==2])
  F11 <- sum(dij[X==1, X==1]) ; F22 <- sum(dij[X==2, X==2])

  if (datatype == "abundance") {
    F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
    F00hat <- ifelse(F22 > 0, ((n-2)* (n-3)* (F11^2)/(4* n* (n-1)* F22)), ((n-2)* (n-3)* (F11*(F11-0.01))/(4 *n * (n-1))) )
  } else if (datatype == "incidence_freq") {
    F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
    F00hat <- ifelse(F22 > 0, ((n-1)^2 * (F11^2)/(4* n* n* F22)), ((n-1)* (n-1)* (F11*(F11-0.01))/(4 *n * n)) )
  }

  if (f0.hat == 0) {
    d = dij
  } else if (f0.hat == 1) {
    d.0bar <- matrix(rep(F.0hat/length(X)/f0.hat, length(X)*f0.hat), length(X), f0.hat)

    d00 = matrix(0, f0.hat, f0.hat)
    d <- cbind(dij, d.0bar )
    aa <- cbind(t(d.0bar), d00 )
    d <- rbind(d, aa)
    diag(d) = 0
  } else {
    d.0bar <- matrix(rep(F.0hat/length(X)/f0.hat, length(X)*f0.hat), length(X), f0.hat)

    fo.num = (f0.hat * (f0.hat-1) )/2
    d00 = matrix(0, f0.hat, f0.hat)
    d00[upper.tri(d00)] = (F00hat/2)/fo.num
    d00 <- pmax(d00, t(d00))###signmatrix
    d <- cbind(dij, d.0bar )
    aa <- cbind(t(d.0bar), d00 )
    d <- rbind(d, aa)
    diag(d) = 0
  }

  return(d)

}

FD.m.est_0 = function (ai_vi, m, q, nT) {
  EFD = function(m, qs, obs, asy, beta, av) {
    m = m - nT
    out <- sapply(1:length(qs), function(i) {
      if (qs[i] != 2) {
        obs[i] + (asy[i] - obs[i]) * (1 - (1 - beta[i])^m)
      }
      else if (qs[i] == 2) {
        V_bar^2/sum((av[, 2]) * ((1/(nT + m)) * (av[, 1]/nT) + ((nT + m - 1)/(nT + m)) * (av[, 1] * (av[, 1] - 1)/(nT * (nT - 1)))))
      }
    })
    return(out)
  }
  V_bar <- sum(ai_vi$ai[, 1] * ai_vi$vi[, 1])/nT
  asy <- iNEXT.3D:::FD_est(ai_vi, q, nT)$est
  obs <- iNEXT.3D:::FD_mle(ai_vi, q)
  out <- sapply(1:ncol(ai_vi$ai), function(i) {
    ai <- ai_vi$ai[, i]
    ai[ai < 1] <- 1
    av = cbind(ai = round(ai), vi = ai_vi$vi[, i])
    RFD_m = iNEXT.3D:::RFD(av, nT, nT - 1, q, V_bar)
    beta <- rep(0, length(q))
    asy_i <- asy[, i]
    obs_i <- obs[, i]
    asy_i <- sapply(1:length(q), function(j) {
      max(asy_i[j], obs_i[j], RFD_m[j])
    })

    obs_i <- sapply(1:length(q), function(j) {
      max(RFD_m[j], obs_i[j])
    })

    beta0plus <- which(asy_i != obs_i)
    beta[beta0plus] <- (obs_i[beta0plus] - RFD_m[beta0plus])/(asy_i[beta0plus] - RFD_m[beta0plus])

    if (sum(m < nT) != 0) {
      int.m = sort(unique(c(floor(m[m < nT]), ceiling(m[m < nT]))))
      mRFD = rbind(int.m, sapply(int.m, function(k) iNEXT.3D:::RFD(av, nT, k, q, V_bar)))
    }

    sapply(m, function(mm) {
      if (mm < nT) {
        if (mm == round(mm)) {
          mRFD[-1, mRFD[1, ] == mm]
        } else {
          (ceiling(mm) - mm) * mRFD[-1, mRFD[1, ] == floor(mm)] + (mm - floor(mm)) * mRFD[-1, mRFD[1, ] == ceiling(mm)]
        }
      }
      else if (mm == nT) {
        obs_i
      }
      else if (mm == Inf) {
        asy_i
      }
      else {
        EFD(m = mm, qs = q, obs = obs_i, asy = asy_i, beta = beta, av = av)
      }
    }) %>% t %>% as.numeric
  })
  matrix(out, ncol = ncol(ai_vi$ai))
}
