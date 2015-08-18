library(igraph)
library(microbenchmark)
library(fdrtool)
network_list = as.matrix(read.table("~/LAS/HumanBinaryHQ_HINT.txt"))

load("/Users/danielqiu/Workspace/Bio/GSE18864_entrez_norm.bin")
load("/Users/danielqiu/Workspace/iLab/Bio/GSE10255_entrez.bin")
g = graph.data.frame(as.matrix(read.table("/Users/danielqiu/Workspace/iLab/Bio/HumanBinaryHQ_HINT.txt")), directed=FALSE)

las(g,b)

e = E(g)


lasr = function(x,y,z)
{
  sum(x*y*z)
}

x = sample(1:100,100,replace=TRUE)
y = sample(1:100,100,replace=TRUE)
z = sample(1:100,100,replace=TRUE)

microbenchmark(
     lascore(x,y,z),
     lasr(x,y,z)
   )


x = b['1',]
y = b['310']
z = b['780']




for(i in 1:500)
{
  if(length(h[[i]])!=0 )
  {
    print(h[[i]][[1]])
  }
  
}



library(GGally)
library(intergraph)
library(network)

testdata =relate.matrix['5700',]
target = names(testdata[testdata==1])

terget.graph = induced.subgraph(graph=g,vids=target)
network = asNetwork(terget.graph)
ggnet(network,size=4,segment.size=1)

network = asNetwork(g)
istarget = rownames(result) %in% target
ggnet(network, size=3,node.group=istarget,color="blue" ,)



testdata =result['5700',]
target = names(testdata[testdata>0.6])
terget.graph = induced.subgraph(graph=g,vids=target)


target = names(testdata)[index.top.N(testdata,N=40)]

target

index.top.N = function(xs, N=10){
  if(length(xs) > 0) {
    o = order(xs, na.last=FALSE)
    o.length = length(o)
    if (N > o.length) N = o.length
    o[((o.length-N+1):o.length)]
  }
  else {
    0
  }
}


testz = result['1000',]
cutoff=0.8
w = names(testz[testz>cutoff])
subg <- induced.subgraph(graph=g,vids=w)
wc = walktrap.community(subg)
network = asNetwork(subg)
ggnet(network,size=4,segment.size=1)


plot(subg)

library(GOstats)


library(GOstats)

sel.entrez<-rownames(b)[1:5]
all.entrez<-rownames(b)
params <- new("GOHyperGParams", geneIds=sel.entrez, universeGeneIds=all.entrez, ontology="BP", pvalueCutoff=0.01,conditional=F, testDirection="over", annotation="hgu133a.db")
Over.pres<-hyperGTest(params)
ov<-summary(Over.pres)
ov$Term





for(c in com)
{
  print(names(com))
}

  sel.entrez = as.character(unlist(neighborhood(g,2,c('10148'))))
  all.entrez<-rownames(result)
  params <- new("GOHyperGParams", geneIds=sel.entrez, universeGeneIds=all.entrez, ontology="BP", pvalueCutoff=0.01,conditional=F, testDirection="over", annotation="hgu133a.db")
  Over.pres<-hyperGTest(params)
  ov<-summary(Over.pres)
ov$Term

typeof(sel.entrez)
unlist(sel.entrez)
sel.entrez
all.entrez



for(ci in community_index)
{
  print(ci)
}
cluster1 <- c("835", "5261","241", "994")
cluster2 <- c("307", "308", "317", "321", "506", "540", "378", "388", "396")
clusterSim(cluster1, cluster2, ont="MF", organism="human", measure="Wang")
temp <- geneSim("5921", "9046", ont = "BP", organism = "human", measure = "Wang", combine = "max")













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
  ggnet(network,node.group=type,size=4,segment.size=1)
}



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
  print(type)
  network = asNetwork(subg)
  ggnet(network,node.group=type,size=4,segment.size=1)
}

xcandidate = c(5696,23212,5696,5696,9683,9683,55739,55739,6634,2892,54903,10480,5685,2892,1534,1337,55739,1534,10480,23075,9312,10480,9683
,3693
,5685
,9040
,23212
,9683
,54903
,65005
,5885
,23075
,1308
,23212
,54903
,6222
,6636
,8541
,5984
,7004
,10746
,5193
,5230
,9040
,10142
,10480
,65005
,7818
,79735
,5885
,8789
,6205
,6222
,11157
,64699
,10480
,54903
,9158
,10142
,10758
,10914
,2815
,3693
,4163
,6205
,23522
,11171
,65005
,9219
,1846
,51065
,4352
,8789
,1008
,6135
,8666
,9219
,6222
,6634
,7135
,10746
,9158
,11171
,4625
,2892
,23204
,5885
,23204
,23204
,56145
,6158
,4355
,8541
,1308
,27258
,5885
,6159
,11157
,5432
,65005
,5984
,5089
,5660
,1337
,10445
,26058
,26762
,2815
,51631
,5660
,5685
,58529
,6636
,7818
,9040
,23204
,11157
,7039
,9158
,1345
,5885
,1058
,5682
,2192
,5885
,6222
,8541
,8065
,54929
,6135
,6222
,6205
,9219
,64699
,6205
,51065
,51371
,57111
,3185
,54929
,1308
,5880
,10746
,1475
,25788
,93974
,9518
,9656
,2971
,3069
,8789
,54929
,10142
,1846
,10382
,10758
,10445
,8668
,2334
,3185
,5230
,54929
,6132
,6430
,6634
,7004
,7316
,6171
,9551
,1442
,56145
,10142
,23204
,2971
,4130
,4355
,51690
,5291
,55739
,6230
,7311
,4215
,1160
,56955
,1160
,1337
,23522
,4969
,29911
,4772
,56145
,10580
,12
,2065
,2192
,27231
,5725
,6137
,80831)

output = visualizewitoutw1(g,result,10181,2,cutoff=2)


for( x in xcandidate)
{
  print(x)
  visualizewitoutw1(g,result10255,x,2,cutoff=1)
}


graph <- erdos.renyi.game(50,0.3)
matrix <- matrix(data=rexp(200,rate=.1), nrow=50, ncol=5, byrow= TRUE, dimnames=NULL)
seq <- seq(from=1,to=50)
row.names(matrix) <-seq
V(graph)$name <- seq
lascouting(graph,matrix,k=2,n.cores=4)

