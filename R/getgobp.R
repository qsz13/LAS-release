#' Generate a result file
#' @param graph
#' @param z.matrix A matrix representing gene Z. Row names are the gene id in gene network.
#' @param k Integer giving the order of the network.
#' @param n.cores Core number used for parallel computing.
#' @param cutoff 
#' @param community.min An Integer giving the min gene numbers in community. 
#' @param term.limit 
#' @return A form contains id of gene w, GO info of gene w , xk.w.semantic.similarity, x.w.avg.distance
#' @example getgobp(,,2,4,0.8,5,NA)
#' @export
#' 
getgobp <- function(graph, z.matrix, k=2, n.cores=4, cutoff=0.8, community.min=5, term.limit=NA)
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
  #   for(x in names(community)){
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
  
  #     xk = neighborhood(graph,k,nodes=x)
  xk = V(graph)[unlist(neighborhood(graph,k,nodes=x))]$name
  #print(xk)
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
  
  w.result = do.call("rbind",lapply(community_index, gen.data, member, xk,x, graph, all.entrez, term.limit))
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





gen.data <- function(ci, member, xk,x,graph, all.entrez, term.limit)
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