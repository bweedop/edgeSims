library(caper)
library(geiger)


drop_random_tips<-function(spp)
{
  tree<-sim.bdtree(n = spp)
  write.tree(tree, file = "original_random_tree.tre")
  
  ed<-ed.calc(tree)$spp
  ed<-ed[order(ed$ED, decreasing = TRUE),]
  rownames(ed)<-1:nrow(ed)
  write.table(ed, file = "original_random_tree_data.txt")
  
  imputed_tree<-tree
  for (i in seq_len(spp/10))
  {
    spp_to_be_dropped<-sample(imputed_tree$tip.label)
    imputed_tree<-drop.tip(imputed_tree, i)
  }
  write.tree(imputed_tree, file = "imputed_random_tree.tre")
  
  imputed_ed<-ed.calc(imputed_tree)$spp
  imputed_ed<-imputed_ed[order(imputed_ed$ED, decreasing = TRUE),]
  rownames(imputed_ed)<-1:nrow(imputed_ed)
  write.table(imputed_ed, file = "imputed_random_tree_data.txt")
}


drop_clustered_tips<-function(spp)
{
  tree<-sim.bdtree(n = spp)
  write.tree(tree, file = "original_clustered_tree.tre")
  
  ed<-ed.calc(tree)$spp
  ed<-ed[order(ed$ED, decreasing = TRUE),]
  rownames(ed)<-1:nrow(ed)
  write.table(ed, file = "original_clustered_tree_data.txt")
  
  imputed_tree<-tree
  brown_evol<-sim.char(tree, 0.05, 1, model = "BM")
  brown_evol_topQ<-summary(brown_evol)[5]
  i_list<-numeric()
  for(i in brown_evol[,,1])
  {
    if(i > brown_evol_topQ)
    {
      imputed_tree<-drop.tip(imputed_tree, paste("s",grep(i, brown_evol[,,1]), sep = ""))
    }
  }
  write.tree(imputed_tree, file = "imputed_clustered_tree.tre")
  
  imputed_ed<-ed.calc(imputed_tree)$spp
  imputed_ed<-imputed_ed[order(imputed_ed$ED, decreasing = TRUE),]
  rownames(imputed_ed)<-1:nrow(imputed_ed)
  write.table(imputed_ed, file = "imputed_clustered_tree_data.txt")
}







