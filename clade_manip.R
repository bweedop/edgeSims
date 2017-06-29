library(geiger)
library(caper)
library(pez)
library(phytools)
get.clades <- function(spp, size)
{
  tree <- sim.bdtree(n = spp)
  clade.mat <- clade.matrix(tree)$clade.matrix
  mat.sums <- rowSums(clade.mat)
  selected.clades <- which(mat.sums == size)
  clade.number<-sample(selected.clades, 1)
  random.clade <- extract.clade(tree, clade.number)
  node.name <- which(names(branching.times(tree)) == clade.number)
  node.age <- branching.times(tree)[node.name]
  nasty.bl <- dist.nodes(tree)[clade.number,][clade.number-1]
  gammatree<-ltt(tree, plot = FALSE, gamma = FALSE)
  gamma <- gammatest(gammatree)$gamma
  print(nasty.bl)
  print(node.age)
}

replace<-function(x, y)
{
  p<-y$tip.label
  for (i in p){
    pez:::bind.replace(x, y, i)
  }
}