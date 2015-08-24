#' Create a table to record Gene Ontology Biological Process mapping results.  Every gene W's community takes a row.
#' 
#' \code{getgobp.community()} generates a result file of ego gene X,  significant GO terms of X, significant GO terms 
#' of genes within k steps of X, gene W, significant GO terms  of W,
# ’ the similarity of gene W and genes within k steps of gene X, the average distance between gene X and gene
# W. ‘ A gene X may correspond with several W communities. Thus one community takes a row in the table.
#' @param graph The graph of gene network.
#' @param z.matrix A matrix representing gene Z (selected scouting genes). Row names are the gene id in gene network.
#' @param k An Integer giving the order of the network.
#' @param n.cores The number of cores used for parallel computing.
#' @param cutoff The threshold to find LA scouting genes.
#' @param community.min Integer. The minimum number of genes numbers in a community. 
#' @param term.limit The maximum number of GO terms to list in a row of the table.
#' @return A table containing the IDs of scouting center genes W, over-represented GO terms by
#'  W, semantic similarity on the Gene Ontology system between the X ego network and all 
#'  scouting center genes, average graph distance between gene X and W. W are grouped by 
#'  network community. Each W community occupies a row. 
#' @export
#' 
#' 
getgobp.community <- function(graph, z.matrix, k = 2, n.cores = 4, cutoff = 1, community.min = 5, term.limit = NA) {
  community <- apply(z.matrix, 1, getCommunity, graph, cutoff, community.min)
  community <- community[!sapply(community, is.null)]
  all.entrez <- colnames(z.matrix)

  cl <- makeCluster(n.cores, outfile = "")
  registerDoParallel(cl)
  cat("loop begin\n")
  
  resulttable <- foreach(i = 1:length(names(community)), .combine = "rbind") %dopar% {
    x <- names(community)[i]
    wc <- community[[x]]
    member <- membership(wc)
    
    community_index <- names(sizes(wc)[sizes(wc) > community.min])
    
    
    sel.entrez <- x
    xgo <- getGO(sel.entrez, all.entrez)
    
    if (is.null(xgo) || is.na(xgo$Pvalue) || length(xgo$Term) == 0) {
      return(NULL)
    } else {
      if (!is.na(term.limit)) {
        xgo <- xgo[1:term.limit, ]
      }
      xgo <- paste(xgo$Term, signif(xgo$Pvalue, digits = 5), sep = ": ", collapse = "\n")
    }
    xk <- V(graph)[unlist(igraph::neighborhood(graph, k, nodes = x))]$name
    
    sel.entrez <- xk
    xkgo <- getGO(sel.entrez, all.entrez)
    
    
    if (is.null(xkgo) || is.na(xkgo$Pvalue) || length(xkgo$Term) == 0) {
      return(NULL)
    } else {
      if (!is.na(term.limit)) {
        xkgo <- xkgo[1:term.limit, ]
      }
      xkgo <- paste(xkgo$Term, signif(xkgo$Pvalue, digits = 5), sep = ": ", collapse = "\n")
    }
    
    w.result <- do.call("rbind", lapply(community_index, get.W.GO, member, xk, x, graph, all.entrez, term.limit))
    if (is.null(w.result)) {
      return(NULL)
    } else {
      print(x)
      return(rbind(resulttable, cbind(x, xgo, xkgo, w.result)))
    }
  }
  stopCluster(cl)
  return(resulttable)
  
}

#' Create a table to record Gene Ontology Biological Process mapping results. Every ego node (X) occupies a row.
#' 
#' \code{getgobp.all()}generates a result file of ego gene X,  significant GO terms of X, significant GO terms
#'  of genes within k steps of X, gene W, significant GO terms  of W,
# ‘ the similarity of gene W and genes within k steps of gene X, the average distance between gene X and gene
# W. ’ A gene X takes a row in the table.
#' @param z.matrix A matrix representing gene Z. (selected scouting genes). Row names are the gene id in gene network.
#' @param k An Integer giving the order of the network.
#' @param n.cores The number of cores used for parallel computing.
#' @param cutoff The threshold to find LA scouting genes.
#' @param community.min An Integer confines the min gene numbers in community. 
#' @param term.limit The maximum number of GO terms to list in a row of the table.
#' @return A table containing the IDs of scouting center genes W, over-represented GO terms by 
#' W, semantic similarity on the Gene Ontology system between the X ego network and all scouting
#'  center genes, average graph distance between gene X and W.
#' @export
#' 
getgobp.all <- function(graph, z.matrix, k = 2, n.cores = 4, cutoff = 1, community.min = 5, term.limit = NA) {
  community <- apply(z.matrix, 1, getCommunity, graph, cutoff, community.min)
  community <- community[!sapply(community, is.null)]
  all.entrez <- colnames(z.matrix)

  cl <- makeCluster(n.cores, outfile = "")
  registerDoParallel(cl)
  cat("loop begin\n")
  
  resulttable <- foreach(i = 1:length(names(community)), .combine = "rbind") %dopar% {
    x <- names(community)[i]
    
    wc <- community[[x]]
    member <- membership(wc)
    
    community_index <- names(sizes(wc)[sizes(wc) > 5])
    
    
    sel.entrez <- x
    xgo <- getGO(sel.entrez, all.entrez)
    
    if (is.null(xgo) || is.na(xgo$Pvalue) || length(xgo$Term) == 0) {
      return(NULL)
    } else {
      if (!is.na(term.limit)) {
        xgo <- xgo[1:term.limit, ]
      }
      xgo <- paste(xgo$Term, signif(xgo$Pvalue, digits = 5), sep = ": ", collapse = "\n")
    }
    xk <- V(graph)[unlist(igraph::neighborhood(graph, k, nodes = x))]$name
    
    sel.entrez <- xk
    xkgo <- getGO(sel.entrez, all.entrez)
   
    if (is.null(xkgo) || is.na(xkgo$Pvalue) || length(xkgo$Term) == 0) {
      return(NULL)
    } else {
      if (!is.na(term.limit)) {
        xkgo = xkgo[1:term.limit, ]
      }
      xkgo <- paste(xkgo$Term, signif(xkgo$Pvalue, digits = 5), sep = ": ", collapse = "\n")
    }
    
    w.result <- do.call("rbind", lapply(community_index, get.W.GO, member, xk, x, graph, all.entrez, term.limit))
    if (is.null(w.result)) {
      return(NULL)
    } else {
      print(x)
      w <- paste(w.result[, 1], collapse = " ")
      wgo <- paste(w.result[, 2], collapse = "\n")
      xk.w.semantic.similarity <- paste(w.result[, 3], collapse = " ")
      x.w.avg.distance <- paste(w.result[, 4], collapse = " ")
      return(rbind(resulttable, cbind(x, xgo, xkgo, w, wgo, xk.w.semantic.similarity, x.w.avg.distance)))
    }
  }
  stopCluster(cl)
  return(resulttable)
}




get.W.GO <- function(ci, member, xk, x, graph, all.entrez, term.limit) {
  
  
  w <- names(member[member == ci])
  
  sel.entrez <- w
  wgo <- getGO(sel.entrez, all.entrez)
  if (is.null(wgo) || is.na(wgo$Pvalue) || length(wgo$Term) == 0) {
    return(NULL)
  } else {
    if (!is.na(term.limit)) {
      wgo <- wgo[1:term.limit, ]
    }
    wgo <- paste(wgo$Term, signif(wgo$Pvalue, digits = 5), sep = ": ", collapse = "\n")
    
  }
  
  
  xk.w.semantic.similarity <- clusterSim(c(w), c(xk), combine = "avg")
  x.w.avg.distance <- mean(shortest.paths(graph, v = w, to = x))
  w <- paste(w, collapse = " ")
  
  return(data.frame(w, wgo, xk.w.semantic.similarity, x.w.avg.distance))
}