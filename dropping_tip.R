library(caper)
library(geiger)
library(plyr)
library(ggplot2)


drop_random_tips<-function(spp, dropped_fraction)
{
  #Simulate phylogeny and write phylogeny data to file.
  tree<-sim.bdtree(n = spp)
  write.tree(tree, file = "original_random_tree.tre")
  
  #Calculate ED values for the original tree and record to .txt file.
  ed<-ed.calc(tree)$spp
  ed<-ed[order(ed$ED, decreasing = TRUE),]
  rownames(ed)<-1:nrow(ed)
  write.table(ed, file = "original_random_tree_data.txt")
  
  #Use the original tree for the tree which will have a determined fraction of random tips dropped from it. 
  imputed_tree<-tree
  for (i in seq_len(spp*dropped_fraction))
  {
    spp_to_be_dropped<-sample(imputed_tree$tip.label, 1) 
    delete_row<-as.numeric(which(ed$species == spp_to_be_dropped))
    ed<-ed[-delete_row,]
    imputed_tree<-drop.tip(imputed_tree, spp_to_be_dropped)
  }
  #Write the manipulated tree to a file for documentation.
  write.tree(imputed_tree, file = "imputed_random_tree.tre")
  
  #Calculate ED values for the manipulated tree and write them to a .txt file for documentation.
  imputed_ed<-ed.calc(imputed_tree)$spp
  imputed_ed<-imputed_ed[order(imputed_ed$ED, decreasing = TRUE),]
  rownames(imputed_ed)<-1:nrow(imputed_ed)
  write.table(imputed_ed, file = "imputed_random_tree_data.txt")
  
  #Calculate the correlation between the ED values of the original tree and the manipulated tree.
  ed_corr<-cor(ed[2], imputed_ed[2])
  return(ed_corr)
}


drop_clustered_tips<-function(spp, dropped_fraction)
{
  #Simulate phylogeny and write phylogeny data to file.
  tree<-sim.bdtree(n = spp)
  write.tree(tree, file = "original_clustered_tree.tre")
  
  #Calculate ED values for the original tree and record to .txt file.
  ed<-ed.calc(tree)$spp
  ordered_ed<-ed[order(ed$ED, decreasing = TRUE),]
  rownames(ordered_ed)<-1:nrow(ordered_ed)
  write.table(ordered_ed, file = "original_clustered_tree_data.txt")
  
  #Use the original tree for the tree which will have a determined fraction of clustered tips dropped from it.
  imputed_tree<-tree
  
  #Simulate continuous trait data on the original tree using a constant rate Brownian-motion model.
  brown_evol<-sim.char(tree, 0.05, 1, model = "BM")[,,1]
  #Assessing which spp fall into the percentile which is to be dropped and the dropping them from the tree.
  quantile <- quantile(brown_evol, 1-dropped_fraction)
  to.drop <- which(brown_evol >= quantile)
  if (dropped_fraction != 0){
    ed<-ed[-to.drop,]
    imputed_tree <- drop.tip(tree, to.drop)
  }else{

  }
  write.tree(imputed_tree, file = "imputed_clustered_tree.tre")
  
  #Calculate the ED values for the manipulated tree and record the data.
  imputed_ed<-ed.calc(imputed_tree)$spp
  ordered_imputed_ed<-imputed_ed[order(imputed_ed$ED, decreasing = TRUE),]
  rownames(ordered_imputed_ed)<-1:nrow(ordered_imputed_ed)
  write.table(ordered_imputed_ed, file = "imputed_clustered_tree_data.txt")
  
  #Calculate the correlation between the original ED values and the manipulated trees' ED values.
  ed_corr<-cor(ed[2], imputed_ed[2])
  return(ed_corr)
}

data_wrapper<-function(method, reps)
{
  #Wrapper that can run either of the methods (random or clustered) and then returns the data and graphs for the
  # amount of repititions set by the user.
  data <- expand.grid(n.spp = c(64, 128, 256, 512, 1024), fraction.dropped = seq(0,0.2,by=0.01), edge.corr=0)
  if (method == "clustered")
  {
    for (i in seq_len(reps))
    {
      for(j in seq_len(nrow(data)))
      {
        data$edge.corr[j] <- data$edge.corr[j] + drop_clustered_tips(data$n.spp[j], data$fraction.dropped[j])
      }
    }
  }else if(method == "random"){
    for (i in seq_len(reps))
    {
      for(j in seq_len(nrow(data)))
      {
        data$edge.corr[j] <- data$edge.corr[j] + drop_random_tips(data$n.spp[j], data$fraction.dropped[j])
      }
    }
  }else{
    print("Invalid method entry")
  }
  
  data[,3]<-data[,3]/reps
  graph<-ggplot(data, aes(fraction.dropped, edge.corr, colour=factor(n.spp))) + geom_line()
  graph+labs(x = "Fraction Dropped", y = "Correlation of ED", colour = "No. of Spp")
  return(c(data, graph))
}




















