#' Visualize: Generate a graph 
#' 
#' 
#' @param g 
#' @param result A function to find gene z
#' @param x Gene X
#' @param k A number of the length of step 
#' @param cutoff A specific number to filter gene w from gene z
#' @example visualize(,,,,0.8)
#' @return 
#' @export 

visualize <- function(g,result, x, k, cutoff=0.8)
{
  
  X = as.character(x)
  Y = V(g)$name[unlist(igraph::neighborhood(g, 2, nodes=X))]
  
  z = result[X,]
  W = names(z[z>cutoff])
  W1 = V(g)$name[unlist(igraph::neighborhood(g, 1, nodes=W))]
  subg = induced.subgraph(g, unique(c(X,Y,W,W1)))
  
  type <- vector(mode="character", length=length(V(subg)))
  type <- setNames(type, V(subg)$name)
  
  
  for(v in W1)
  {
    type[v] = "W neighbor"
  }
  for(v in W)
  {
    type[v] = "W"
  }
  
  for(v in Y)
  {
    type[v] = "Y"
  }
  type[X] = "X"
  print(type)
  network = asNetwork(subg)
  ggnet(network,node.group=type,segment.size=1,label.nodes=T,col="white")
}

#' Visualize without w1: Generate a graph contains gene x,y,w1
#' @param g 
#' @param result A function to find gene z
#' @param x Gene X
#' @param k A number of the length of step 
#' @param cutoff A specific number to filter gene w from gene z
#' @example visualize(,,,,0.8)
#' @return  
#' @export
#'
visualizewitoutw1 <- function(g,result, x, k, cutoff=0.8)
{
  
  X = as.character(x)
  Y = V(g)$name[unlist(igraph::neighborhood(g, 2, nodes=X))]
  
  z = result[X,]
  W = names(z[z>cutoff])
  #W1 = V(g)$name[unlist(igraph::neighborhood(g, 1, nodes=W))]
  
  #print(c(X,Y,W,W1))
  subg = induced.subgraph(g, unique(c(X,Y,W)))
  
  type <- vector(mode="character", length=length(V(subg)))
  type <- setNames(type, V(subg)$name)
  
  
  
  for(v in W)
  {
    type[v] = "W"
  }
  for(v in Y)
  {
    type[v] = "Y"
  }
  type[X] = "X"
  #print(type)
  network = asNetwork(subg)
  
  output =  ggnet(network,node.group=type,segment.size=1,label.nodes=T,col="black")
  ggsave(output, file=paste("test/",as.character(x),".jpg",sep = ""), w=8, h=6, scale=5,limitsize=FALSE)
}


