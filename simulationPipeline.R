require(ape)
require(geiger)

simulate_trees<-function(species)
{
  tree <- sim.bdtree(n = species)
  tree$edge.length <- tree$edge.length / (max(branching.times(tree)) * 10)
  
  return(tree)
}

tree_files <- function(trees, species)
{
  number_of_trees <- c(rep(1:trees))
  for (i in number_of_trees)
  {
    tree <- simulate_trees(species)
    write.tree(tree, file = paste(i,"simulatedTree.tre", sep = "_"))
  }
}

run_dawg<-function(trees, species, seq_length, model, format)
{
  tree_files(trees, species)
  user_length <- paste("Length =", seq_length)
  alt_model<- paste('"',model,'"', sep = "")
  user_model <- paste("Model =", noquote(alt_model))
  user_format <- paste("Format =", paste('"', format, '"', sep = ""))
  number_of_trees <- c(rep(1:trees))
  for (i in number_of_trees)
  {
    user_file <- paste("File =", paste('"',i,"simulatedSeqs.fasta",'"',sep = ""))
    tree <- readLines(paste(i,"simulatedTree.tre", sep = "_"))
    tree <- paste("Tree =", noquote(tree))
    writeLines(paste(tree, user_length, user_model, user_file, user_format, sep = "\n\n"), 
               "simulationParams.dawg")
    system("dawg simulationParams.dawg")
  }
  
}

