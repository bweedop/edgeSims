library(caper)
library(geiger)
library(pez)
library(phytools)
library(moments)
library(apTreeshape)

drop_random_tips<-function(spp, dropped_fraction)
{
  #Simulate phylogeny.
  tree<-sim.bdtree(n = spp)
  
  #Calculate ED values for the original tree.
  ed <-ed.calc(tree)$spp
  ed <-setNames(ed[,2], ed[,1])

  original.lambda <- yule(tree)$lambda
  gammatree <-ltt(tree, plot = FALSE, gamma = FALSE)
  original.gamma <- gammatest(gammatree)$gamma
  colless.tree <- as.treeshape(tree, model = "yule")
  original.colless <- colless(colless.tree)
  original.kurtosis <- kurtosis(original.ed)
  original.skew <- skewness(original.ed)
  original.branches <- sum(tree$edge.length)
  
  #Use the original tree for the tree which will have a determined fraction of random tips dropped from it. 
  imputed_tree<-tree
  for (i in seq_len(spp*dropped_fraction))
  {
    spp_to_be_dropped<-sample(imputed_tree$tip.label, 1) 
    ed <- ed[!names(ed) %in% spp_to_be_dropped]
    imputed_tree<-drop.tip(imputed_tree, spp_to_be_dropped)
  }
  
  #Calculate ED values for the manipulated tree.
  imputed_ed<-ed.calc(imputed_tree)$spp
  imputed_ed<-setNames(imputed_ed[,2], imputed_ed[,1])
  
  #Calculate the correlation between the ED values of the original tree and the manipulated tree.
  ed_corr<-cor(ed, imputed_ed)
  return(c(ed_corr, original.gamma, original.lambda, original.colless, 
            original.kurtosis, original.skew, original.sd, original.branches))
}


drop_clustered_tips<-function(spp, dropped_fraction)
{
  #Simulate phylogeny.
  tree<-sim.bdtree(n = spp)
  
  #Calculate ED values for the original tree.
  ed<-ed.calc(tree)$spp
  ed <- setNames(ed[,2], ed[,1])

  original.lambda <- yule(tree)$lambda
  gammatree <-ltt(tree, plot = FALSE, gamma = FALSE)
  original.gamma <- gammatest(gammatree)$gamma
  colless.tree <- as.treeshape(tree, model = "yule")
  original.colless <- colless(colless.tree)
  original.kurtosis <- kurtosis(original.ed)
  original.skew <- skewness(original.ed)
  original.branches <- sum(tree$edge.length)
  
  #Use the original tree for the tree which will have a determined fraction of clustered tips dropped from it.
  imputed_tree<-tree
  
  #Simulate continuous trait data on the original tree using a constant rate Brownian-motion model.
  brown_evol<-sim.char(imputed_tree, 0.05, 1, model = "BM")[,,1]
  #Assessing which spp fall into the percentile which is to be dropped and the dropping them from the tree.
  quantile <- quantile(brown_evol, 1-dropped_fraction)
  to.drop <- which(brown_evol >= quantile)
  if (dropped_fraction != 0){
    imputed_tree <- drop.tip(tree, to.drop)
  }else{

  }
  ed <- ed[!names(ed) %in% names(to.drop)]
  #Calculate the ED values for the manipulated tree.
  imputed_ed<-ed.calc(imputed_tree)$spp
  imputed_ed<-setNames(imputed_ed[,2], imputed_ed[,1])
  
  #Calculate the correlation between the original ED values and the manipulated trees' ED values.
  ed_corr<-cor(ed, imputed_ed)
  return(c((ed_corr, original.gamma, original.lambda, original.colless, 
            original.kurtosis, original.skew, original.sd, original.branches))
}

data_wrapper<-function(n.spp = c(64, 128, 256, 512, 1024, 2048, 4096), fraction.dropped = c(seq(0.01, 0.99, 0.01)), reps = 100)
{
  #Wrapper that can run either of the methods (random or clustered) and then returns the data and graphs for the
  # amount of repititions.
  data <- expand.grid(n.spp = n.spp, fraction.dropped = fraction.dropped, reps = 0:reps, clustered.corr= NA, random.corr = NA)
  for(i in seq_len(nrow(data))){
        random.temp.vec <- drop_random_tips(data$n.spp[i], data$fraction.dropped[i])
	      clustered.temp.vec <- drop_clustered_tips(data$n.spp[i], data$fraction.dropped[i])
        
        data$clustered.corr[i] <- clustered.temp.vec
	      data$random.corr[i] <- random.temp.vec
    }
  write.table(data, "dropCorrelations.txt")
}