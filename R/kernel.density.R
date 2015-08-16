#' Evaluate the result using kernel density estimation.
#' 
#' @param relate.matrix The matrix returned by lascouting.
#' @param network.graph The igraph object representing the gene network.
#' @param smoothing.normalize Different ways to normalize the result.
#' @return A matrix representing the kernel density of each gene. Each row is a gene, columns are the weights of scouting genes for the gene.
#' @examples
#'  kernel.density(,,one)
#'  kernel.density(,,squareM)
#'  kernel.density(,,none)
#' @export
#' 
kernel.density <- function(relate.matrix, network.graph, smoothing.normalize=c("one","squareM","none") ) {
  smoothing.normalize <- match.arg(smoothing.normalize)
  
  network.node <- V(network.graph)$name
  matrix.node <- row.names(relate.matrix)
  if(!identical(intersect(network.node,matrix.node),union(network.node,matrix.node))){
    common.node <- getCommonNode(network.graph, relate.matrix)
    network.graph <- cleanGraph(network.graph, common.node)
  }
  
  weight0 = dnorm(0)
  weight1 = dnorm(1)
  weight2 = dnorm(2)
  
  size = nrow(relate.matrix)
  
  
  relate.matrix = relate.matrix[order(rownames(relate.matrix)), ] 
  relate.matrix = relate.matrix[,order(colnames(relate.matrix)) ] 
  
  
  
  
  
  adjacency1 <- get.adjacency(network.graph, type="both")
  adjacency2 <- get.adjacency(connect.neighborhood(network.graph,2), type="both")-adjacency1
  
  adjacency1 = adjacency1[order(rownames(adjacency1)), ] 
  adjacency1 = adjacency1[,order(colnames(adjacency1)) ] 
  adjacency2 = adjacency2[order(rownames(adjacency2)), ] 
  adjacency2 = adjacency2[,order(colnames(adjacency2)) ] 
  
  temp = diag(size)*weight0 
  weight.matrix = temp+adjacency1*weight1+adjacency2*weight2
  
  if(smoothing.normalize=="one")
  {
    rsum = rowSums(as.matrix(weight.matrix))
    nmatrix = diag(1/rsum)
    colnames(nmatrix) <- rownames(weight.matrix)
    weight.matrix = nmatrix %*% weight.matrix
  }
  else if(smoothing.normalize=="squareM")
  {
    temp.weight.matrix = as.matrix(weight.matrix)
    rsum = rowSums(temp.weight.matrix)
    m = rowSums(temp.weight.matrix != 0)
    nmatrix = diag(sqrt(m)/rsum)
    colnames(nmatrix) <- rownames(weight.matrix)
    weight.matrix = nmatrix%*%weight.matrix
    
  }
  
  result <- relate.matrix%*%weight.matrix
  
  
  return(result)
  
}