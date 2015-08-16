#' Find the liquid association scouting gene
#' @useDynLib LAS
#' @param network.graph An igraph object representing the gene network.
#' @param express.matrix A matrix represeting the express matrix for the genes in gene network.Row names are the gene id in gene network.
#' @param k Integer giving the order of the network.
#' @param n.cores Core number used for parallel computing.
#' @return A logical matrix representing the LA-scouting genes for each gene. Rows represent the center gene id and columns represents the LA-scouting genes.
#' @example lascouting(,,2,4)
#' @export
#' 
lascouting <- function(network.graph, express.matrix, k=2, n.cores=4){
  
  network.node <- V(network.graph)$name
  matrix.node <- row.names(express.matrix)
  if(!identical(intersect(network.node,matrix.node),union(network.node,matrix.node))){
    common.node <- getCommonNode(network.graph, express.matrix)
    network.graph <- cleanGraph(network.graph, common.node)
    express.matrix <- cleanMatrix(express.matrix, common.node)
  }
  size <- length(common.node)
  express.matrix = normalizeInputMatrix(express.matrix)
  if(k!=1)
  {
    graph.connected <- connect.neighborhood(network.graph,k)
    connected.list <- as.matrix(get.edgelist(graph.connected))
  }
  else
  {
    connected.list <- as.matrix(get.edgelist(network.graph))
  }
  row.size <- nrow(connected.list)
  express.matrix.t <- t(express.matrix)/ncol(express.matrix)
  
  cl <- makeCluster(n.cores, outfile="")
  registerDoParallel(cl)
  
  ptm <- proc.time()
  
  cat("loop begin\n")
  
  result <- foreach(i=1:row.size) %dopar%
{
  xy <- express.matrix[connected.list[i,1],]*express.matrix[connected.list[i,2],]
  la.vector <- c(xy%*%express.matrix.t)
  lfdr <- fdrtool(la.vector, verbose=FALSE, plot = FALSE)$lfdr
  return(rownames(express.matrix)[which(lfdr<0.2)])
  
}
stopCluster(cl)
print(proc.time() - ptm)

node.z <- Matrix(0, nrow = size, ncol = size,dimnames=list(rownames(express.matrix),rownames(express.matrix)))
for(i in 1:row.size)
{
  if(length(result[[i]])!=0 )
  {
    x = connected.list[i,1]
    y = connected.list[i,2]
    node.z[x,c(result[[i]])] = 1
    node.z[y,c(result[[i]])] = 1
  }
}

result = node.z#z.kernel.density(node.z, network.graph)
return(result)

}

