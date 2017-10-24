library(geiger)
library(caper)
library(pez)
library(phytools)
library(moments)
library(apTreeshape)

find.clade <- function(tree, tips)
{
  clade.mat <- clade.matrix(tree)$clade.matrix
  sums <- rowSums(clade.mat)
  ancestors <- apply(clade.mat, 1, function(x) all(x[tips]==1))
  sums[!ancestors] <- length(tree$tip.label)
  return(unname(which.min(sums)))
}


bind.replace <- function(backbone, donor, replacing.tip.label)
{
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

get.ed <- function(spp, size)
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
  #Getting the spp of the focal clade so that the ED values of these spp can be selected out of 
  #original.ed and imputed.ed.
  random.clade <- extract.clade(tree, clade.number)
  focal.clade.spp <- unique(na.omit(as.numeric(unlist(strsplit(unlist(random.clade$tip.label), "[^0-9]+")))))
  
  original.ed<-ed.calc(tree)$spp
  original.ed.ranking <- original.ed[order(original.ed$ED),]
  original.ed <- setNames(original.ed[,2], original.ed[,1])
  excluding.clade.original <- original.ed[!names(original.ed) %in% random.clade$tip.label]
  node.name <- which(names(branching.times(tree)) == clade.number)
  
  #Calculating for original.tree
  original.lambda <- yule(tree)$lambda
  node.age <- branching.times(tree)[node.name]
  nasty.bl <- dist.nodes(tree)[clade.number,][clade.number-1]
  gammatree<-ltt(tree, plot = FALSE, gamma = FALSE)
  original.gamma <- gammatest(gammatree)$gamma
  colless.tree <- as.treeshape(tree, model = "yule")
  original.colless <- colless(colless.tree)
  original.kurtosis <- kurtosis(original.ed)
  original.skew <- skewness(original.ed)
  original.sd <- sd(original.ed)
  original.branches <- sum(tree$edge.length)

  dropped.tree<-drop.tip(tree, sample(random.clade$tip.label, length(random.clade$tip.label)-1))
  
  r<-NULL
  while(is.null(r)){
    try({
      donor.clade<-sim.bdtree(n = size)
      donor.clade$tip.label <- random.clade$tip.label
      for (i in dropped.tree$tip.label){
        if(any(random.clade$tip.label==i)){
          tip<-i
        }
      }
      r<-bind.replace(dropped.tree, donor.clade, tip)
      dropped.tree<-r})
  }

  #calculating for imputed.tree
  imputed.colless.tree <- as.treeshape(dropped.tree, model = "yule")
  imputed.colless <- colless(imputed.colless.tree)
  imputed.lambda <- yule(dropped.tree)$lambda
  imputed.gammatree<-ltt(dropped.tree, plot = FALSE, gamma = FALSE)
  imputed.gamma <- gammatest(imputed.gammatree)$gamma
  imputed.ed <- ed.calc(dropped.tree)$spp
  imputed.ed.ranking <- imputed.ed[order(imputed.ed$ED),]
  imputed.ed <- setNames(imputed.ed[,2], imputed.ed[,1])
  imputed.kurtosis <- kurtosis(imputed.ed)
  imputed.skew <- skewness(imputed.ed)
  imputed.sd <- sd(imputed.ed)
  imputed.branches <- sum(dropped.tree$edge.length)

  #calculating for imputed clade (donor.clade)
  imputed.clade.lambda <- yule(donor.clade)$lambda
  original.clade.lambda <- yule(random.clade)$lambda
  imputed.clade.branches <- sum(donor.clade$edge.length)
  original.clade.branches <- sum(random.clade$edge.length)

  full.ranking.compare <- sum(abs(original.ed.ranking$ED - imputed.ed.ranking$ED))/((spp-size)/size)
  error.params <- c(50,100,200,0.05*spp,0.1*spp,0.2*spp)
  error.rate <- NA
  for (i in error.params){
    error.rate[which(error.params == i)] <- sum(!(imputed.ed.ranking$species[1:i] %in% original.ed.ranking$species[1:i]))/i
  }

  excluding.clade.imputed <- imputed.ed[!names(imputed.ed) %in% random.clade$tip.label]
  full.ed.corr <- cor(excluding.clade.original, excluding.clade.imputed)
  focal.ed.corr <- cor(original.ed[random.clade$tip.label], imputed.ed[random.clade$tip.label])
  
  return(c(full.ed.corr, focal.ed.corr, original.gamma, imputed.gamma, original.lambda, imputed.lambda, 
            original.colless, imputed.colless, original.kurtosis, imputed.kurtosis, original.skew, imputed.skew, 
            original.sd, imputed.sd, imputed.clade.lambda, original.clade.lambda, original.branches, imputed.branches,
            original.clade.branches, imputed.clade.branches, error.rate))
}

wrapper<-function(n.spp=c(64, 128, 256, 512, 1024), clade.size=c(3, 4, 8, 16), reps = 100)
{
    data <- expand.grid(n.spp = n.spp, clade.size = clade.size, reps=1:reps, full.ed=NA, focal.ed=NA, 
    original.gamma = NA, imputed.gamma = NA, original.lambda = NA, imputed.lambda = NA,  original.colless = NA, 
    imputed.colless = NA, original.kurtosis = NA, imputed.kurtosis = NA, original.skew = NA, imputed.skew = NA, 
    original.sd = NA, imputed.sd = NA, original.branches = NA, imputed.branches = NA, error.rate.50 = NA, 
    error.rate.100 = NA ,error.rate.200 = NA, error.rate.5pct = NA, error.rate.10pct = NA, error.rate.20pct = NA)
    for(i in seq_len(nrow(data))){
        temp.vec <- get.ed(data$n.spp[i], data$clade.size[i])
        data$full.ed[i] <- temp.vec[1]
        data$focal.ed[i] <- temp.vec[2]
        data$original.gamma[i] <- temp.vec[3]
        data$imputed.gamma[i] <- temp.vec[4]
        data$original.lambda[i] <- temp.vec[5]
        data$imputed.lambda[i] <- temp.vec[6]
        data$original.colless[i] <- temp.vec[7]
        data$imputed.colless[i] <- temp.vec[8]
        data$original.kurtosis[i] <- temp.vec[9]
        data$imputed.kurtosis[i] <- temp.vec[10]
        data$original.skew[i] <- temp.vec[11]
        data$imputed.skew[i] <- temp.vec[12]
        data$original.sd[i] <- temp.vec[13]
        data$imputed.sd[i] <- temp.vec[14]
        data$imputed.clade.lambda[i] <- temp.vec[15]
        data$original.clade.lambda[i] <- temp.vec[16]
        data$original.branches[i] <- temp.vec[17]
        data$imputed.branches[i] <- temp.vec[18]
        data$original.clade.branches[i] <- temp.vec[19]
        data$imputed.clade.branches[i] <- temp.vec[20]
        data$error.rate.50[i] <- temp.vec[21]
        data$error.rate.100[i] <- temp.vec[22]
        data$error.rate.200[i] <- temp.vec[23]
        data$error.rate.5pct[i] <- temp.vec[24]
        data$error.rate.10pct[i] <- temp.vec[25]
        data$error.rate.20pct[i] <- temp.vec[26]
    }

    write.table(data, "correlations.txt")
}

wrapper()