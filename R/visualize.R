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
visualize <- function(g,result, x, k, cutoff=1, path=NULL)
{
  
  X = as.character(x)
  Y = V(g)$name[unlist(igraph::neighborhood(g, 2, nodes=X))]
  
  z = result[X,]
  W = names(z[z>cutoff])

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

  network = asNetwork(subg)
  size = length(V(subg))
  scale = (17/961)*size + 1999/961
  
  output =  ggnet(network,node.group=type,segment.size=1,label.nodes=T,col="black",subset.threshold = 1)
  ggsave(output, file=paste(as.character(x),".jpg",sep = ""), path=path,w=4, h=3, scale=scale,limitsize=FALSE)
  return(output)
}

#' visualize with w community
#' @export
#'
visualize.with.community<-function(g,result, x, k=2, cutoff=1,community.min=5,path=NULL)
{
  X = as.character(x)
  Y = V(g)$name[unlist(igraph::neighborhood(g, k, nodes=X))]
  z = result[X,]

  W = names(z[z>cutoff])
  
  wc = getCommunity(z, g,cutoff,  community.min)
  member = membership(wc)
  community_index = names(sizes(wc)[sizes(wc)>community.min])
  if(length(community_index)>7)
  {
    print("comunity too large.")
    return()
  }
  
  subg = induced.subgraph(g, unique(c(X,Y,W)))
  
  type <- rep("other", length(V(subg)))
  type <- setNames(type, V(subg)$name)

  time = 1
  for(ci in community_index)
  {
    w = names(member[member==ci])
    type[w]=paste("w",time,sep =  "")
    time<-time+1
  }

  type[Y] = "Y"

  type[X] = "X"
  
  network = asNetwork(subg)
  size = length(V(subg))
  scale = (17/961)*size + 1999/961

  output =  ggnet(network,node.group=type,segment.size=1,label.nodes=T,col="black",subset.threshold = 1)
  
  ggsave(output, file=paste(as.character(x),".jpg",sep = ""), path=path,w=4, h=3, scale=scale,limitsize=FALSE)
  return(output)
}