#' Generate a result file of gene x,  gostats of x, gostats of genes within k steps of x, gene w, gostats of w, 
#' the similarity of gene w and genes within k steps of gene x, the average distance between gene x and gene w.
#' A gene x may correspond with several w communities. Thus one community takes a row in the table.
#' @param graph The graph of gene network.
#' @param z.matrix A matrix representing gene Z. Row names are the gene id in gene network.
#' @param k An Integer giving the order of the network.
#' @param n.cores A Core number used for parallel computing.
#' @param cutoff A number used to find LA scouting gene z.
#' @param community.min An Integer confines the min gene numbers in community. 
#' @param term.limit A parameter indicates there is no limit of content in a line of the table.
#' @return A form contains id of gene w, GO info of gene w , xk.w.semantic.similarity, x.w.avg.distance.
#' @example getgobp(graph, z_matrix, k=2, n.cores=4, cutoff=1, community.min=5, term.limit=NA)
#' @export
#' 
getgobp <- function(graph, z.matrix, k=2, n.cores=4, cutoff=1, community.min=5, term.limit=NA)
{
  
  community = apply(z.matrix, 1, getCommunity, graph,cutoff,community.min)
  community = community[!sapply(community, is.null)]
  all.entrez<-colnames(z.matrix)
  
  resulttable = NULL
  cl <- makeCluster(n.cores, outfile="")
  registerDoParallel(cl)
  cat('loop begin\n')
  
  resulttable <- foreach(i=1:length(names(community)), .combine='rbind') %dopar%
  {
    x = names(community)[i]
    wc = community[[x]]
    member = membership(wc)
    
    community_index = names(sizes(wc)[sizes(wc)>community.min])
    
    
    sel.entrez<-x
    xgo = getGO(sel.entrez, all.entrez)
    
    if(is.null(xgo)||is.na(xgo$Pvalue)||length(xgo$Term)==0)
    {
      return(NULL)
    }  
    else
    {
      if(!is.na(term.limit))
      {
        xgo = xgo[1:term.limit,]
      }
      xgo <- paste(xgo$Term, signif(xgo$Pvalue,digits = 5), sep=": ", collapse = '\n')
    }
    xk = V(graph)[unlist(igraph::neighborhood(graph,k,nodes=x))]$name

    sel.entrez = xk
    xkgo = getGO(sel.entrez, all.entrez)
    
    
    if(is.null(xkgo)||is.na(xkgo$Pvalue)||length(xkgo$Term)==0)
    {
      return(NULL)
    }
    else
    {
      if(!is.na(term.limit))
      {
        xkgo = xkgo[1:term.limit,]
      }
      xkgo <- paste(xkgo$Term, signif(xkgo$Pvalue,digits = 5), sep=": ", collapse = '\n')
    }
    
    w.result = do.call("rbind",lapply(community_index, get.W.GO, member, xk,x, graph, all.entrez, term.limit))
    if(is.null(w.result))
    {
      return(NULL)
    }
    else
    {
      print(x)
      return(rbind(resulttable,cbind(x, xgo, xkgo,w.result)))
    }
  }
  stopCluster(cl)
  return(resulttable)

}

#' Generate a result file of gene x,  gostats of x, gostats of genes within k steps of x, gene w, gostats of w, 
#' the similarity of gene w and genes within k steps of gene x, the average distance between gene x and gene w.
#' Regardless of a gene X may correspond with multiple w communities. A gene X only takes a row in the table.
#' @param graph The graph of gene network.
#' @param z.matrix A matrix representing gene Z. Row names are the gene id in gene network.
#' @param k An Integer giving the order of the network.
#' @param n.cores A Core number used for parallel computing.
#' @param cutoff A number used to find LA scouting gene z.
#' @param community.min An Integer confines the min gene numbers in community. 
#' @param term.limit A parameter indicates there is no limit of content in a line of the table.
#' @return A form contains id of gene w, GO info of gene w , semantic similarity of xk and gene w, average distance between gene x and w.
#' @example getgobp.x.in.one.line(graph, z.matrix, k=2, n.cores=4, cutoff=1, community.min=5, term.limit=NA)
#' @export
#' 
getgobp.x.in.one.line <- function(graph, z.matrix, k=2, n.cores=4, cutoff=1, community.min=5, term.limit=NA)
{
  community = apply(z.matrix, 1, getCommunity, graph,cutoff,community.min)
  community = community[!sapply(community, is.null)]
  all.entrez<-colnames(z.matrix)
  
  resulttable = NULL
  cl <- makeCluster(n.cores, outfile="")
  registerDoParallel(cl)
  cat('loop begin\n')
  
  resulttable <- foreach(i=1:length(names(community)), .combine='rbind') %dopar%
  {
    x = names(community)[i]
    
    wc = community[[x]]
    member = membership(wc)
    
    community_index = names(sizes(wc)[sizes(wc)>5])
    
    
    sel.entrez<-x
    xgo = getGO(sel.entrez, all.entrez)
    
    if(is.null(xgo)||is.na(xgo$Pvalue)||length(xgo$Term)==0)
    {
      return(NULL)
    }  
    else
    {
      if(!is.na(term.limit))
      {
        xgo = xgo[1:term.limit,]
      }
      xgo <- paste(xgo$Term, signif(xgo$Pvalue,digits = 5), sep=": ", collapse = '\n')
    }
    xk = V(graph)[unlist(igraph::neighborhood(graph,k,nodes=x))]$name
    
    sel.entrez = xk
    xkgo = getGO(sel.entrez, all.entrez)
    
    
    if(is.null(xkgo)||is.na(xkgo$Pvalue)||length(xkgo$Term)==0)
    {
      return(NULL)
    }
    else
    {
      if(!is.na(term.limit))
      {
        xkgo = xkgo[1:term.limit,]
      }
      xkgo <- paste(xkgo$Term, signif(xkgo$Pvalue,digits = 5), sep=": ", collapse = '\n')
    }
    
    w.result = do.call("rbind",lapply(community_index, get.W.GO, member, xk,x, graph, all.entrez, term.limit))
    if(is.null(w.result))
    {
      return(NULL)
    }
    else
    {
      print(x)
      w = paste(w.result[,1], collapse =" ")
      wgo = paste(w.result[,2], collapse ="\n")
      xk.w.semantic.similarity = paste(w.result[,3], collapse =" ")
      x.w.avg.distance = paste(w.result[,4], collapse =" ")
      return(rbind(resulttable,cbind(x, xgo, xkgo,w, wgo,xk.w.semantic.similarity,x.w.avg.distance )))
    }
  }
  stopCluster(cl)
  return(resulttable)
}




get.W.GO <- function(ci, member, xk,x,graph, all.entrez, term.limit)
{
  
  
  w = names(member[member==ci])
  
  sel.entrez<-w
  wgo = getGO(sel.entrez, all.entrez)
  if(is.null(wgo)||is.na(wgo$Pvalue)||length(wgo$Term)==0)
  {
    return(NULL)
  }
  else
  {
    if(!is.na(term.limit))
    {
      wgo = wgo[1:term.limit,]
    }
    wgo <- paste(wgo$Term, signif(wgo$Pvalue,digits = 5), sep=": ", collapse = '\n')
    
  }
  
  
  xk.w.semantic.similarity = clusterSim(c(w),c(xk),combine = "avg")
  x.w.avg.distance = mean(shortest.paths(graph,v=w,to=x))
  w<-paste(w, collapse = ' ')
  
  return(data.frame(w, wgo , xk.w.semantic.similarity, x.w.avg.distance))
}