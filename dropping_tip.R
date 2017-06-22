library(caper)
library(geiger)


drop_random_tips<-function(spp, dropped_fraction)
{
  tree<-sim.bdtree(n = spp)
  write.tree(tree, file = "original_random_tree.tre")
  
  ed<-ed.calc(tree)$spp
  ed<-ed[order(ed$ED, decreasing = TRUE),]
  rownames(ed)<-1:nrow(ed)
  write.table(ed, file = "original_random_tree_data.txt")
  
  imputed_tree<-tree
  for (i in seq_len(spp*dropped_fraction))
  {
    spp_to_be_dropped<-sample(imputed_tree$tip.label, 1)
    delete_row<-as.numeric(which(ed$species == spp_to_be_dropped))
    ed<-ed[-delete_row,]
    imputed_tree<-drop.tip(imputed_tree, spp_to_be_dropped)
  }
  write.tree(imputed_tree, file = "imputed_random_tree.tre")
  
  imputed_ed<-ed.calc(imputed_tree)$spp
  imputed_ed<-imputed_ed[order(imputed_ed$ED, decreasing = TRUE),]
  rownames(imputed_ed)<-1:nrow(imputed_ed)
  write.table(imputed_ed, file = "imputed_random_tree_data.txt")
  
  ed_corr<-cor(ed[2], imputed_ed[2])
  return(ed_corr)
}


drop_clustered_tips<-function(spp, dropped_fraction)
{
  tree<-sim.bdtree(n = spp)
  write.tree(tree, file = "original_clustered_tree.tre")
  
  ed<-ed.calc(tree)$spp
  ordered_ed<-ed[order(ed$ED, decreasing = TRUE),]
  rownames(ordered_ed)<-1:nrow(ordered_ed)
  write.table(ordered_ed, file = "original_clustered_tree_data.txt")
  
  imputed_tree<-tree
  brown_evol<-sim.char(tree, 0.05, 1, model = "BM")[,,1]
  quantile <- quantile(brown_evol, 1-0.1)
  to.drop <- which(brown_evol >= quantile)
  ed<-ed[-to.drop,]
  imputed_tree <- drop.tip(tree, to.drop)
  write.tree(imputed_tree, file = "imputed_clustered_tree.tre")
  
  imputed_ed<-ed.calc(imputed_tree)$spp
  imputed_ed<-imputed_ed[order(imputed_ed$ED, decreasing = TRUE),]
  rownames(imputed_ed)<-1:nrow(imputed_ed)
  write.table(imputed_ed, file = "imputed_clustered_tree_data.txt")
  
  ed_corr<-cor(ed[2], imputed_ed[2])
  return(ed_corr)
}

data_wrapper<-function(method)
{
  data <- expand.grid(n.spp=c(200,400,800,1000,1200), fraction.dropped=c(.01,.05,.1,.2), edge.corr=NA)
  if (method == "clustered")
  {
    for(i in seq_len(nrow(data)))
    {
      data$edge.corr[i] <- drop_clustered_tips(data$n.spp[i], data$fraction.dropped[i])
    }
  }else if(method == "random"){
    for(i in seq_len(nrow(data)))
    {
      data$edge.corr[i] <- drop_random_tips(data$n.spp[i], data$fraction.dropped[i])
    }
  }else{
    print("Invalid method entry")
  }

  return(data)
}

