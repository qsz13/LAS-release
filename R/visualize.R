#' Visualize: Generate a graph which vividly displays the gene x,y and w.
#' 
#' \code{visualize()} generates a graph. It is used to intuitively and vividly display the layout of gene x, y and w. 
#' 
#' @param graph The graph of gene network.
#' @param kernel.result The result of kernel.density which finds genes z of a gene x.
#' @seealso \code{\link{kernel.result}}
#' @param x The Gene the graph is generated for.
#' @param k A specific number stands for the length of step from gene y to gene x.
#' @param cutoff A specific number to filter gene w from gene z.
#' @example 
#' # Create sample data for examples.
#' library(stats)
#' library(igraph)
#' graph <- erdos.renyi.game(50,0.3)
#' relate_matrix <- matrix(data=rexp(200,rate=.1), nrow=50, ncol=5, byrow= TRUE, dimnames=NULL)
#' #use the first normalize method as an example
#' kernel.result <- kernel.density(relate_matrix, graph, smoothing.normalize=c("one"))
#' visualize(graph,kernel.result,x,k=2,cutoff=1,path=NULL)
#' @return a graph of gene x,y and w
#' @export 
visualize <- function(graph,kernel.result, x, k=2, cutoff=1, path=NULL)
{
  
  X = as.character(x)
  Y = V(graph)$name[unlist(igraph::neighborhood(graph, k, nodes=X))]
  
  z = kernel.result[X,]
  W = names(z[z>cutoff])

  subg = induced.subgraph(graph, unique(c(X,Y,W)))

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

#' visualize with w community: Generate a graph of genes w and their community in different colors.
#' @param graph The fraph of gene network.
#' @param kernel.result The result of kernel.density which finds genes z of a gene x.
#' @seealso \code{\link{kernel.result}} 
#' @param x The Gene the graph is generated for.
#' @param k A specific number stands for the length of step from gene y to gene x.
#' @param cutoff A specific number to filter gene w from gene z.
#' @param cummunity.min An Integer confines the least number of genes in a community of w shown in graph. 
#' @param path The path where the result graph is saved to.The default path is the original path of input graph.
#' @return a graph displays genes w and their correspongding community in different colors.
#' @example 
#' # Create sample data for examples.
#' library(stats)
#' library(igraph)
#' graph <- erdos.renyi.game(50,0.3)
#' relate_matrix <- matrix(date=rexp(200,rate=.1), nrow=50, ncol=5, byrow= TURE, dimnames=NULL)
#' #use the first normalize method as an example
#' kernel.result <- kernel.density(relate_matrix, graph, smoothing.normalize=c("one"))
#' visualize.with.community(graph,kernel.result,x,k=2,cutoff=1,community.min=5,path=NULL)
#'
#' @export
#'
visualize.with.community<-function(graph,kernel.result, x, k=2, cutoff=1,community.min=5,path=NULL)
{
  X = as.character(x)
  Y = V(graph)$name[unlist(igraph::neighborhood(graph, k, nodes=X))]
  z = kernel.result[X,]

  W = names(z[z>cutoff])
  
  wc = getCommunity(z, graph,cutoff,  community.min)
  member = membership(wc)
  community_index = names(sizes(wc)[sizes(wc)>community.min])
  if(length(community_index)>7)
  {
    print("comunity too large.")
    return()
  }
  
  subg = induced.subgraph(graph, unique(c(X,Y,W)))
  
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