library(geiger)
library(caper)
library(pez)
library(phytools)
library(hexbin)
library(ggplot2)


find.clade <- function(tree, tips){
  clade.mat <- clade.matrix(tree)$clade.matrix
  sums <- rowSums(clade.mat)
  ancestors <- apply(clade.mat, 1, function(x) all(x[tips]==1))
  sums[!ancestors] <- length(tree$tip.label)
  return(unname(which.min(sums)))
}


bind.replace <- function(backbone, donor, replacing.tip.label){
  size <- max(branching.times(donor))
  bind.point <- which(backbone$tip.label == replacing.tip.label)
  tree <- bind.tree(backbone, donor, where=bind.point)
  clade <- find.clade(tree, which(tree$tip.label %in% donor$tip.label))
  which.edge <- which(tree$edge[,2] == clade)
  if(tree$edge.length[which.edge] < size)
    stop("Donor phylogeny too large for backbone phylogeny")
  tree$edge.length[which.edge] <- tree$edge.length[which.edge] - size
  return(tree)
}


get.clades <- function(spp, size)
{
  size.met <- FALSE
  while (!size.met)
  {
    tree <- sim.bdtree(n = spp)
    clade.mat <- clade.matrix(tree)$clade.matrix
    mat.sums <- rowSums(clade.mat)
    size.met <- any(mat.sums == size)
  }
  selected.clades <- which(mat.sums == size)
  if (length(selected.clades) == 1){
    clade.number <- selected.clades
  }else{
    clade.number<-sample(selected.clades, 1)
  }
  random.clade <- extract.clade(tree, clade.number)
  original.ed<-ed.calc(tree)$spp[,2]
  node.name <- which(names(branching.times(tree)) == clade.number)
  node.age <- branching.times(tree)[node.name]
  nasty.bl <- dist.nodes(tree)[clade.number,][clade.number-1]
  gammatree<-ltt(tree, plot = FALSE, gamma = FALSE)
  gamma <- gammatest(gammatree)$gamma
  
  dropped.tree<-drop.tip(tree, sample(random.clade$tip.label, length(random.clade$tip.label)-1))

  r<-NULL
  attempt<-0
  while(is.null(r)){
    attempt<-attempt+1
    try({
      donor.clade<-sim.bdtree(n = size)
      donor.clade$tip.label <- letters[seq_along(donor.clade$tip.label)]
      for (i in dropped.tree$tip.label){
        if(any(random.clade$tip.label==i)){
          tip<-i
        }
      }
      r<-bind.replace(dropped.tree, donor.clade, tip)
      dropped.tree<-r})
}
  imputed.ed<-ed.calc(dropped.tree)$spp[,2]
  ed.corr<-cor(original.ed, imputed.ed)
  return(c(nasty.bl,node.age,gamma,ed.corr))
}


wrapper<-function(n.spp=seq(600, 1200, by = 200), clade.size=seq(2, 10,by=1)){
  params <- expand.grid(n.spp=n.spp, clade.size=clade.size)
  data <- matrix(nrow=nrow(params), ncol=4)
  for(i in seq_len(nrow(params))){
    data[i,] <- get.clades(params$n.spp[i], params$clade.size[i])
  }
  data <- cbind(params, data)
  colnames(data)<-c("n.spp","clade.size","nasty.bl","node.age","gamma","ed.corr")
  ageplot<-ggplot(data, aes(clade.size, node.age))
  ageplot+geom_hex()
  return(data)
}


get.full.ed <- function(spp, size)
{
  #Running the function in order to get the ED value correlations for the full phylogenies.
  size.met <- FALSE
  while (!size.met)
  {
    tree <- sim.bdtree(n = spp)
    clade.mat <- clade.matrix(tree)$clade.matrix
    mat.sums <- rowSums(clade.mat)
    size.met <- any(mat.sums == size)
  }
  selected.clades <- which(mat.sums == size)
  if (length(selected.clades) == 1){
    clade.number <- selected.clades
  }else{
    clade.number<-sample(selected.clades, 1)
  }
  random.clade <- extract.clade(tree, clade.number)
  focal.clade.spp <- unique(na.omit(as.numeric(unlist(strsplit(unlist(random.clade$tip.label), "[^0-9]+")))))
  #Getting the spp of the focal clade so that the ED values of these spp can be selected out of 
  #original.ed and imputed.ed.
  original.ed<-ed.calc(tree)$spp
  original.ed <- setNames(original.ed[,2], original.ed[,1])
  node.name <- which(names(branching.times(tree)) == clade.number)
  node.age <- branching.times(tree)[node.name]
  nasty.bl <- dist.nodes(tree)[clade.number,][clade.number-1]
  gammatree<-ltt(tree, plot = FALSE, gamma = FALSE)
  gamma <- gammatest(gammatree)$gamma
  
  dropped.tree<-drop.tip(tree, sample(random.clade$tip.label, length(random.clade$tip.label)-1))
  
  r<-NULL
  while(is.null(r)){
    try({
      donor.clade<-sim.bdtree(n = size)
      donor.clade$tip.label <- random.clade$tip.label#letters[seq_along(donor.clade$tip.label)]
      for (i in dropped.tree$tip.label){
        if(any(random.clade$tip.label==i)){
          tip<-i
        }
      }
      r<-bind.replace(dropped.tree, donor.clade, tip)
      dropped.tree<-r})
  }
  imputed.ed<-ed.calc(dropped.tree)$spp
  imputed.ed <- setNames(imputed.ed[,2], imputed.ed[,1])
  full.ed.corr <- cor(original.ed[tree$tip.label], imputed.ed[tree$tip.label])
  focal.ed.corr<-cor(original.ed[focal.clade.spp], imputed.ed[focal.clade.spp])
  return(full.ed.corr)
}

get.focal.ed <- function(spp, size)
{
  #Running the function in order to get the ED value correlations for the focal clades.
  
  size.met <- FALSE
  while (!size.met)
  {
    tree <- sim.bdtree(n = spp)
    clade.mat <- clade.matrix(tree)$clade.matrix
    mat.sums <- rowSums(clade.mat)
    size.met <- any(mat.sums == size)
  }
  selected.clades <- which(mat.sums == size)
  if (length(selected.clades) == 1){
    clade.number <- selected.clades
  }else{
    clade.number<-sample(selected.clades, 1)
  }
  random.clade <- extract.clade(tree, clade.number)
  focal.clade.spp <- unique(na.omit(as.numeric(unlist(strsplit(unlist(random.clade$tip.label), "[^0-9]+")))))
  #Getting the spp of the focal clade so that the ED values of these spp can be selected out of 
    #original.ed and imputed.ed.
  original.ed<-ed.calc(tree)$spp
  original.ed <- setNames(original.ed[,2], original.ed[,1])
  node.name <- which(names(branching.times(tree)) == clade.number)
  node.age <- branching.times(tree)[node.name]
  nasty.bl <- dist.nodes(tree)[clade.number,][clade.number-1]
  gammatree<-ltt(tree, plot = FALSE, gamma = FALSE)
  gamma <- gammatest(gammatree)$gamma
  
  dropped.tree<-drop.tip(tree, sample(random.clade$tip.label, length(random.clade$tip.label)-1))
  
  r<-NULL
  while(is.null(r)){
    try({
      donor.clade<-sim.bdtree(n = size)
      donor.clade$tip.label <- random.clade$tip.label#letters[seq_along(donor.clade$tip.label)]
      for (i in dropped.tree$tip.label){
        if(any(random.clade$tip.label==i)){
          tip<-i
        }
      }
      r<-bind.replace(dropped.tree, donor.clade, tip)
      dropped.tree<-r})
  }
  imputed.ed<-ed.calc(dropped.tree)$spp
  imputed.ed <- setNames(imputed.ed[,2], imputed.ed[,1])
  full.ed.corr <- cor(original.ed[tree$tip.label], imputed.ed[tree$tip.label])
  focal.ed.corr<-cor(original.ed[focal.clade.spp], imputed.ed[focal.clade.spp])
  return(focal.ed.corr)
}


multiple.wrapper<-function(n.spp= c(64, 128, 256, 512, 1024), clade.size=c(2, 3, 4, 8, 16, 32), reps = 100){
  
  params <- expand.grid(n.spp=n.spp, clade.size=clade.size)
  
  full.data <- matrix(nrow=nrow(params), ncol=reps)
  for(i in seq_len(nrow(params))){
    for (j in seq_len(ncol(full.data))){
      full.data[i,j] <- get.full.ed(params$n.spp[i], params$clade.size[i])
    }
  }
  
  focal.data <- matrix(nrow=nrow(params), ncol=reps)
  for(i in seq_len(nrow(params))){
    for (j in seq_len(ncol(focal.data))){
      focal.data[i,j] <- get.focal.ed(params$n.spp[i], params$clade.size[i])
    }
  }
  
  
  colnames(full.data)<-paste("trial.", 1:reps, sep = "")
  colnames(focal.data)<-paste("trial.", 1:reps, sep = "")
  
  full.data <- cbind(params, full.data)
  focal.data <- cbind(params, focal.data)
  
  write.table(full.data, file = "FullCorrelationData.txt")
  write.table(focal.data, file = "FocalCorrelationData.txt")
}



