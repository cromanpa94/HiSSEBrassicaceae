setwd("~/Desktop/Brassicaceae/ClusterJan12/PIP_example")
library(data.table)
library(ape)
library(geiger)
library(hisse)
library(qpcR)
library(ape)
library(rotl)
library(dismo)
library(letsR)
library(ape)
library(plyr)
library(geiger)
library(phangorn)
library(splitstackshape)
library(raster)
library(rgdal)
library(data.table)
library(dplyr)
library(ggplot2)
library(raster)
library(grid)
library(taxize)
library(treebase)
library(caper)
library(spocc)
library("rgeos")
library(maptools)
library(raster)
library(rgdal)
library(sdm)
library(seegSDM)
library(virtualspecies)
library(rangeBuilder)
require(gridExtra)
library(data.table)
library(parallel)
library(maptools)
library(randomForest)
library(phylopath)
library(picante)
library(e1071)   
library(diversitree)

Ploidy<-list.files(".", ".tre$", recursive = T)


Ploidy<-list.files(".", "ploidy.txt", recursive = T)
Ploidy2<-Ploidy
Ploidy2<-Ploidy2[-c(43)]
Ploidy2<-Ploidy2[-c(66)]
Ploidy2<-Ploidy2[-c(69)]
Ploidy2<-Ploidy2[-c(74)]
Ploidy2<-Ploidy2[-c(82)]
Ploidy2<-Ploidy2[-c(115)]


files<-lapply(1:length(Ploidy2), function(x){
  print(x)
  read.delim(Ploidy2[x], header=FALSE, comment.char="#")
  
})

names(files)<-Ploidy2
DF<-rbindlist(files, idcol   = T, fill = T)

spp<-as.character(unique(DF$V1))

incongruence<-lapply(seq_along(spp), function(x){
  length(unique(na.omit(DF[DF$V1 == spp[x], "V5"])))
  
})
which(unlist(incongruence)==2)

DF<-DF[!is.na(DF$V5),]

dupli<-DF[!duplicated(DF$V1),]


###Running HiSSE
setwd("~/Desktop/Brassicaceae/ClusterJan12")
load("Script_ploidy_Feb14-ChromEvol.RData")

tree<-read.tree("Brassicaceae.tree")

data<-as.matrix(dupli$V5)
row.names(data)<-dupli$V1
tree2<-drop.tip(tree, name.check(tree, data)$tree_not_data)
data2<-cbind.data.frame(dupli$V1,dupli$V5)


write.csv(ntree$tip.label, "tips.csv")

tips_corr<-read.csv("tips_corrected.csv")

spp_drop<-as.character(tips_corr[which(tips_corr$Drop == "Delete"),"Tip"])


tree2<-drop.tip(tree2,spp_drop)  
tree<-drop.tip(tree,spp_drop)  
data2<-data2[-which(data2$`dupli$V1` %in% spp_drop),]





##Run HiSSE and BiSSE
  
  ntree<-tree2
  nnls<-nnls.tree(cophenetic(ntree),ntree,rooted=TRUE)
  
  hdbin_first<-data2  
  
  h<-length(which(hdbin_first[,2] == 0))
  c<-length(which(hdbin_first[,2] == 1))
  sap<-c(h,c)/length(hdbin_first[,2])
  
  
##BiSSE analyses
  
  trans.rates.bisse <- ParEqual(TransMatMaker(hidden.states = FALSE), c(1, 2))
  ##Bisse 1; Tau=Equal, Epsilon=Equal,q=Equal
  trans.rates.bisse_eq<-trans.rates.bisse
  trans.rates.bisse_eq[!is.na(trans.rates.bisse_eq) & !trans.rates.bisse_eq == 0] = 1
  
  b1 <- hisse(ntree, hdbin_first, f=sap, hidden.states = FALSE,
              turnover.anc = c(1,1,0,0), eps.anc = c(1,1,0,0),
              trans.rate = trans.rates.bisse_eq, output.type="raw")
  
  ##Bisse 2; Tau=Free, Epislon=Free, q=Equal
  b2 <- hisse(ntree, hdbin_first, f=sap, hidden.states = FALSE,
              turnover.anc = c(1,1,0,0), eps.anc = c(1,2,0,0),
              trans.rate = trans.rates.bisse_eq, output.type="raw")
  
  
  ##Bisse 3; Tau=equal, Epislon=equal, q=free
  b3 <- hisse(ntree, hdbin_first, f=sap, hidden.states = FALSE,
              turnover.anc = c(1,2,0,0), eps.anc = c(1,1,0,0),
              trans.rate = trans.rates.bisse, output.type="raw")
  
  ##Bisse 4; Tau=free, Epislon=free, q=free
  b4 <- hisse(ntree, hdbin_first, f=sap, hidden.states = FALSE,
                         turnover.anc = c(1,2,0,0), eps.anc = c(1,2,0,0),
                         trans.rate = trans.rates.bisse, output.type="raw")
  
##HiSSE analyses  

  #hisse 1 Hidden state present for both 0 & 1; Tau=Free, Epsilon=Free,q=Equal
  trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
  trans.rates.hisse[!is.na(trans.rates.hisse) & !trans.rates.hisse == 0] = 1
  trans.rates.hisse
  h1 <- hisse(ntree, hdbin_first, f=sap ,hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4),
               trans.rate=trans.rates.hisse, output.type="raw")

  #hisse 2 Hidden state present for 0; Tau=Free, Epsilon=Free,q=Equal
  trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
  trans.rates.hisse <- ParDrop(trans.rates.hisse, c(10,11,12,3,6,9,8,5))
  trans.rates.hisse[!is.na(trans.rates.hisse) & !trans.rates.hisse == 0] = 1
  trans.rates.hisse
  h2 <- hisse(ntree, hdbin_first,f=sap , hidden.states=TRUE, turnover.anc=c(1,2,3,0), eps.anc=c(1,2,3,0),
              trans.rate=trans.rates.hisse, output.type="raw")
  
  #hisse 3 Hidden state present for 1; Tau=Free, Epsilon=Free,q=Equal
  trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
  trans.rates.hisse <- ParDrop(trans.rates.hisse, c(2,3,5,7,8,9,10,12)) 
  trans.rates.hisse[!is.na(trans.rates.hisse) & !trans.rates.hisse == 0] = 1
  trans.rates.hisse
  h3 <- hisse(ntree, hdbin_first,f=sap , hidden.states=TRUE, turnover.anc=c(1,2,0,3), eps.anc=c(1,2,0,3),
              trans.rate=trans.rates.hisse, output.type="raw")
  
  #hisse 4 Hidden state present for both 0 & 1; Tau=Free, Epsilon=Free,q=Free
  trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
  h4 <- hisse(ntree, hdbin_first, f=sap ,hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4),
              trans.rate=trans.rates.hisse, output.type="raw")
  
  #hisse 5 Hidden state present for 0; Tau=Free, Epsilon=Free,q=Free
  trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
  trans.rates.hisse <- ParDrop(trans.rates.hisse, c(10,11,12,3,6,9,8,5))
  h5 <- hisse(ntree, hdbin_first,f=sap , hidden.states=TRUE, turnover.anc=c(1,2,3,0), eps.anc=c(1,2,3,0),
              trans.rate=trans.rates.hisse, output.type="raw")
  
  #hisse 6 Hidden state present for 1; Tau=Free, Epsilon=Free,q=Free
  trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
  trans.rates.hisse <- ParDrop(trans.rates.hisse, c(2,3,5,7,8,9,10,12))
  trans.rates.hisse
  h6 <- hisse(ntree, hdbin_first,f=sap , hidden.states=TRUE, turnover.anc=c(1,2,0,3), eps.anc=c(1,2,0,3),
              trans.rate=trans.rates.hisse, output.type="raw")
  
  #hisse CID-2
  trans.rates = TransMatMaker(hidden.states=TRUE)
  trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
  cid2_1 = hisse(ntree, hdbin_first, f=sap, hidden.states=TRUE, turnover.anc=c(1,1,2,2),
             eps.anc=c(1,1,2,2), trans.rate=trans.rates.nodual)
  
  trans.rates = TransMatMaker(hidden.states=TRUE)
  trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
  trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
  trans.rates.nodual.allequal
  cid2_2 = hisse(ntree, hdbin_first, f=sap, hidden.states=TRUE, turnover.anc=c(1,1,2,2),
                 eps.anc=c(1,1,2,2), trans.rate=trans.rates.nodual)
  

  trans.rates.nodual.threerates <- trans.rates.nodual
  to.change <- cbind(c(1,3), c(2,4))
  trans.rates.nodual.threerates[to.change] = 1
  to.change <- cbind(c(2,4), c(1,3))
  trans.rates.nodual.threerates[to.change] = 2
  to.change <- cbind(c(1,3,2,4), c(3,1,4,2))
  trans.rates.nodual.threerates[to.change] = 3
  trans.rates.nodual.threerates
  
  cid2_3 = hisse(ntree, hdbin_first, f=sap, hidden.states=TRUE, turnover.anc=c(1,1,2,2),
                 eps.anc=c(1,1,2,2), trans.rate=trans.rates.nodual.threerates)
  
  #hisse CID-4
  cid4 <- hisse.null4(ntree, hdbin_first, f=sap, turnover.anc=rep(c(1,2,3,4),2), 
                               eps.anc=rep(c(1,2,3,4),2), trans.type="equal")
  
  
  models<-c(paste0("BiSSE_", 1:4), paste0("HiSSE_", 1:6), "HiSSE-CID2_1","HiSSE-CID2_2","HiSSE-CID2_3", "HiSSE-CID4")
  
 Table<- cbind.data.frame(Model=models, 
        lnLik=c(b1$loglik,b2$loglik,b3$loglik,b4$loglik,
                h1$loglik,h2$loglik,h3$loglik,h4$loglik,h5$loglik,h6$loglik,
                cid2_1$loglik,cid2_2$loglik,cid2_3$loglik,cid4$loglik),
  AIC=c(b1$AIC,b2$AIC,b3$AIC,b4$AIC,h1$AIC,h2$AIC,h3$AIC,h4$AIC,h5$AIC,h6$AIC,cid2_1$AIC,
        cid2_2$AIC,cid2_3$AIC,cid4$AIC))
 Table<-Table[-c(12:13),]
 
 Table$deltaAIC<- akaike.weights(Table$AIC)$deltaAIC
 Table$weights<- akaike.weights(Table$AIC)$weights
 
 write.csv(Table, "Hisse.table.csv")
 
 ##Marginal reconstruction
 h4_rec <- MarginRecon(ntree, hdbin_first, f=sap, hidden.states=TRUE, pars=h4$solution,
                        aic=h4$aic, n.cores=4)

 ##Support region
 h4_sr <- SupportRegion(h4, scale.int=0.05)
 h4_sr2 <- SupportRegion(h4, scale.int=0.05, output.type = "raw")
 
 
 
boxplot(h4_sr2$points.within.region[,c(2:5)]) #Speciation
boxplot(h4_sr2$points.within.region[,c(6:9)]) #Extinction
boxplot(h4_sr2$points.within.region[,c(2:5)]-h4_sr2$points.within.region[,c(6:9)]) #Diversification
boxplot(h4_sr2$points.within.region[,c(6:9)]/h4_sr2$points.within.region[,c(2:5)]) #Turnover


colMeans(h4_sr2$points.within.region[,c(2:5)]-h4_sr2$points.within.region[,c(6:9)])


 mr<-GetModelAveRates(h4_rec, type = "tips")

write.csv(mr, "Table S4.csv")
  

 mean(mr[which(mr$state ==0),"net.div"])
 mean(mr[which(mr$state ==1),"net.div"])
 
 mean(mr[which(mr$state ==0),"speciation"])
 mean(mr[which(mr$state ==1),"speciation"])
 
 mean(mr[which(mr$state ==0),"extinction"])
 mean(mr[which(mr$state ==1),"extinction"])
 
 mean(mr[which(mr$state ==0),"turnover"])
 mean(mr[which(mr$state ==1),"extinction"])
 
###Clades----
 Phylo_Gen_Summary<-function(tree=NULL){
   genera<-  unique(sub("_.*", "", tree$tip.label))
   
   spp_vec<-c()
   spp_vec_con<-c()
   for (i in 1:length(genera)){
     spp_vec<-tree$tip.label[sub("_.*", "", tree$tip.label) %in% genera[i]]
     spp_vec_con<-rbind(spp_vec_con,paste(spp_vec, collapse = ', '))
   }
   
   
   spp_vec_csin<-cSplit(spp_vec_con, 'V1', sep=", ", type.convert=FALSE)
   spp_vec_csin<-as.matrix(spp_vec_csin)
   
   
   ##Obtain node number
   
   noden<-c()
   node<-c()
   for (i in 1:length(genera)){
     node<-getMRCA(tree, na.omit(spp_vec_csin[i,]))
     node_1<-replace(node,which(is.null(node)),0) 
     noden<-c(noden,node_1)
   }
   
   ##Organize genus and node number
   organiza<-cbind(genera, noden)
   
   ###extract and obtain Nnumber and Ntips
   tree$node.label<-(length(tree$tip.label)+1: c(length(tree$tip.label)+tree$Nnode))
   
   nnodes<-c()
   nntips<-c()
   formanodes<-c()
   formatips<-c()
   
   for (i in 1:dim(organiza)[1]){
     stree<-extract.clade(tree, replace(organiza[,2][i], organiza[,2][i]==0, tree$node.label[1]))
     nnodes<-stree$Nnode
     nntips<-length(stree$tip.label)
     formanodes<-c(formanodes,nnodes)
     formatips<-c(formatips,nntips)
   }
   
   newdata<-cbind(organiza, formanodes, formatips)
   newdata<-newdata[,-1]
   
   if(is.null(dim(newdata) )== TRUE){
     newdata= t(as.data.frame(newdata))
     rownames(newdata)<-genera
     
   }else {
     
     newdata[formanodes ==(tree$node.label[1]-2)]<-0 
     newdata[formatips ==(tree$node.label[1]+1)] <-0 
     rownames(newdata)<-genera
     
   }  
   ##Number of species in each genus
   
   sppnum<-c()
   for (i in 1:dim(spp_vec_csin)[1]){
     ri<-length(na.omit(spp_vec_csin[i,]))
     sppnum<-c(sppnum,ri)
   }
   
   newdata<-cbind(newdata, sppnum)
   newdata<-as.data.frame(newdata)
   nsppintree<- replace(as.vector(newdata$formatips), as.vector(newdata$formatips)==0, 1) 
   newdata$formatips<-nsppintree
   
   
   ##IF ELSE
   
   mon<-c()
   monoph<-c()
   for (i in 1:dim(newdata)[1]){
     mon<-if(newdata[i,3] == newdata[i,4]) "Mono" else "Non-Mono"
     monoph<-rbind(monoph,mon)
   }
   
   newdata<-cbind(newdata, monoph)
   
   newdataf<-as.data.frame(newdata)
   #newdata_2 <-newdataf[as.numeric(newdataf[,3]) < c(2*as.numeric(as.vector(newdataf[,4]))), ] ##Include paraphyletic
   newdata_2 <- subset(newdata , monoph == "Mono") ##Without paraphyletic
   colnames(newdata_2)<-c("node number", "# nodes", "# tips (clade)","# tips (species)", "Monophyly")
   
   newdata_2[newdata_2==0] <- NA
   
   print(newdata_2)
   return(newdata_2)
   
 }
 Phylo_Get_Richness<-function(newdata_2=NULL, database="col"){
   
   rich<-c()
   richness_1<-c()
   
   pudf<-row.names(newdata_2)
   
   for (i in 1: length(pudf) ){
     rich<-dim(downstream(pudf[i], db = database, downto = 'Species', verbose = FALSE, rows=1)[[1]])[1]
     rich<-replace(rich,which(is.null(rich[1])),0) 
     richness_1<-c(richness_1,rich)
     cat("Getting richness for genus ", i, "of", length(pudf), "\n") 
     flush.console()
   }
   
   #richness<- pmax(richness_1, richness_2, na.rm = TRUE)
   richness<- richness_1
   newdata_2<-cbind(newdata_2, richness) ##Remember this change
   
   richnessr<-c()
   for (i in 1:length(newdata_2$`node number`)){
     richnessr[i]<- ifelse(as.numeric(as.character(newdata_2$`# tips (species)`[i]))>newdata_2$richness[i],as.numeric(as.character(newdata_2$`# tips (species)`[i])) ,newdata_2$richness[i] )
   }
   
   newdata_2$richness<-richnessr
   
   print(newdata_2)
   return(newdata_2)
 }
 Phylo_Get_Div<-function(tree,data){
   
   tree$node.label<-NULL
   
   cladeage<-function(phylo, cladelist){
     
     library(ape)
     library(phytools)
     
     uniquefam<-unique(cladelist[,1]) # this generates a vector of unique clade names
     
     # This loop generates a set of vectors with the family names of Acanthomorphs in the tree and assigns species names to them which I will use for estimating MRCA of all families in all the trees.
     
     familynames<-list()
     for(i in 1:length(uniquefam)) {
       familynames[[i]]<-row.names(cladelist)[cladelist ==as.character(uniquefam[i])]
     }
     
     
     maxcrowncladeage<-vector(length=length(familynames))
     mincrowncladeage <-vector(length=length(familynames))
     mediancrowncladeage <-vector(length=length(familynames))
     maxstemcladeage <-vector(length=length(familynames))
     minstemcladeage <-vector(length=length(familynames))
     medianstemcladeage <-vector(length=length(familynames))
     percentparaphyly<-vector(length=length(familynames))
     
     # this loops over the clade MRCA function in phytools and calculates the age of the crown and stem if you have more than 1 species in a clade or just the stem age if you only have 1 species
     
     allcrownages<-matrix(nrow=length(phylo), ncol=length(uniquefam))
     colnames(allcrownages)<-uniquefam
     
     allstemages<-matrix(nrow=length(phylo), ncol=length(uniquefam))
     colnames(allstemages)<-uniquefam
     
     paraphyly<-matrix(nrow=length(phylo), ncol=length(uniquefam))
     colnames(paraphyly)<-uniquefam
     
     
     for(i in 1:length(familynames)){
       
       print(i)
       
       if (length(familynames[[i]])>1){
         
         for (j in 1:length(phylo)) {
           print(j)
           
           bt<-branching.times(phylo[[j]])
           foocrown<-findMRCA(phylo[[j]], familynames[[i]])			
           allcrownages[j,i]<-bt[names(bt)==foocrown]
           stemnodemulti<-phylo[[j]]$edge[,1][phylo[[j]]$edge[,2]==foocrown]
           allstemages[j,i]<-bt[names(bt)==stemnodemulti]
           len<-length(phylo[[j]]$tip.label) # finds the length of the tip labels so can identify tip nodes when trying to identify paraphyly
           foodesc<-getDescendants(phylo[[j]], foocrown)
           cladetips<-foodesc[foodesc<=len] # this identifies just the descendants that are tips
           paraphyly[j,i]<-if(length(familynames[[i]])<length(cladetips)) 1 else 0 # then compares the length of the all the descendant tips to that of the number of species in the clade if there are fewer species in the clade then it is paraphyletic (1) and if not then the clade is monophyletic (0) 
         }
         
         
       } else {
         
         for (j in 1:length(phylo)) {
           bt<-branching.times(phylo[[j]])
           footip<-match(familynames[[i]], phylo[[j]]$tip.label)
           stemnode<-phylo[[j]]$edge[,1][phylo[[j]]$edge[,2]==footip]
           allstemages[j,i]<-bt[names(bt)== stemnode]
         }	
         
         
       }# this will summarize the m age for each tree across trees as a summary
       
       maxcrowncladeage[i]<-max(allcrownages[,i])	
       mincrowncladeage[i]<-min(allcrownages[,i]) 
       mediancrowncladeage[i]<-median(allcrownages[,i])
       
       maxstemcladeage[i]<-max(allstemages[,i])	
       minstemcladeage[i]<-min(allstemages[,i]) 
       medianstemcladeage[i]<-median(allstemages[,i])
       
       percentparaphyly[i]<-(sum(paraphyly[,i])/length(paraphyly[,i]))*100
       
       summaryages<-cbind(mediancrowncladeage,mincrowncladeage,maxcrowncladeage, medianstemcladeage, minstemcladeage, maxstemcladeage, percentparaphyly)
       row.names(summaryages)<-as.character(uniquefam)
       
       res<-list(summaryages, allcrownages, allstemages, paraphyly)
     }
     
     res
     
   }
   genera<-as.data.frame(sapply(strsplit(tree$tip.label,"_"),function(x) x[1]))
   row.names(genera)<-tree$tip.label
   
   ca<-cladeage(list(tree, tree), genera)
   crown<-ca[[2]][1,]
   stem<-ca[[3]][1,]
   dates<-cbind.data.frame(stem,crown)
   
   new_data<-merge(dates, b, by=0)
   
   rates<-do.call("cbind.data.frame",list(sapply(1:dim(new_data)[1], function(x){
     bd.ms(time=new_data$stem[x], n=new_data$richness[x], crown = F, epsilon = 0)
   }),
   
   sapply(1:dim(new_data)[1], function(x){
     bd.ms(time=new_data$stem[x], n=new_data$richness[x], crown = F, epsilon = 0.5)
     
   }),
   sapply(1:dim(new_data)[1], function(x){
     bd.ms(time=new_data$stem[x], n=new_data$richness[x], crown = F, epsilon = 0.9)
     
   }),
   sapply(1:dim(new_data)[1], function(x){
     bd.ms(time=new_data$stem[x], n=new_data$richness[x], crown = F, epsilon = 0.99)
     
   }),
   sapply(1:dim(new_data)[1], function(x){
     bd.ms(time=new_data$crown[x], n=new_data$richness[x], crown = T, epsilon = 0)
     
   }),
   sapply(1:dim(new_data)[1], function(x){
     bd.ms(time=new_data$crown[x], n=new_data$richness[x], crown = T, epsilon = 0.5)
     
   }),
   sapply(1:dim(new_data)[1], function(x){
     bd.ms(time=new_data$crown[x], n=new_data$richness[x], crown = T, epsilon = 0.9)
     
   }),
   
   sapply(1:dim(new_data)[1], function(x){
     bd.ms(time=new_data$crown[x], n=new_data$richness[x], crown = T, epsilon = 0.99)
     
   })
   
   
   ))
   
   names(rates)<-c("S0","S5","S9","S99","C0","C5","C9","C99")
   q<-cbind.data.frame(new_data,rates)
   row.names(q)<-q[,1]
   return(q)
 }
 
 a<-Phylo_Gen_Summary(tree)
 b<-Phylo_Get_Richness(a)
 c<-Phylo_Get_Div(tree,b)
 
 
 jetzDivRates <- function(tree) {
   
   spRate <- function(sp, tree) {
     #get branch lengths from root to tip
     edges <- vector()
     daughterNode <- match(sp, tree$tip.label)
     while (daughterNode != (length(tree$tip.label) + 1)) {
       parentNode <- tree$edge[which(tree$edge[,2] == daughterNode), 1]
       edges <- c(edges, tree$edge.length[which(tree$edge[,1] == parentNode & tree$edge[,2] == daughterNode)])
       daughterNode <- parentNode
     }
     
     res <- sum(sapply(1:length(edges), function(x) edges[x] * (1/(2 ^ (x-1)))))
     res <- res ^ (-1)
     
     return(res)
   }
   
   rates <- unlist(lapply(tree$tip.label, function(x) spRate(x, tree)))
   names(rates) <- tree$tip.label
   
   return(rates)
   
 }
 DR<-jetzDivRates(tree)
 DR<-as.data.frame(DR)
 
 write.csv(DR, "DR.csv")
 
 
 Add<-lapply(1:length(c[,1]), function(x){
  s<- DR[ which(gsub("\\_.*","",row.names(DR)) == row.names(c)[x] ),]
   cbind.data.frame( DR_mean=mean(s), DR_skewness=
                     skewness(s)  )
 })
 Add<-do.call("rbind",Add)
 row.names(Add)<-row.names(c)
 
 c<-merge(c, Add, by=0)
 row.names(c)<-row.names(Add)
 
 df<-hdbin_first
 df$genus<-gsub("_.*","\\1",hdbin_first$`dupli$V1`)
 

 fp<-lapply( 1:length(c[,1]), function(i){
   ric<- c[i,"richness"]
  su<-df[which( df$genus ==  row.names(c)[i]  ),]
  prop<-if(dim(su)[1]==0){0}else{length(which(su$`dupli$V5` == 1))/dim(su)[1]}
  propr<-if(dim(su)[1]==0){0}else{length(which(su$`dupli$V5` == 1))/ric}
  Absolute<-if(dim(su)[1]==0){0}else{length(which(su$`dupli$V5` == 1))}
  
  cbind.data.frame(Tip_proportion=prop,
  Richness_proportion=propr,
  Absolute=Absolute)
 })
 
 
fp<-do.call("rbind",fp)
row.names(fp)<-row.names(c)

data_ppa <- merge(c, fp, by=0)[,-1]
colnames(data_ppa)
row.names(data_ppa)<-data_ppa$Row.names


models_stem <- define_model_set(
MS3_1   = c(richnessL ~ S9, richnessL ~ stem),
MS3_2   = c(richnessL ~ S9, richnessL ~ stem, S9 ~ Absolute),
MS3_3   = c(richnessL ~ S5, richnessL ~ stem),
MS3_4   = c(richnessL ~ S5, richnessL ~ stem, S5 ~ Absolute),
MS3_5   = c(richnessL ~ S0, richnessL ~ stem),
MS3_6   = c(richnessL ~ S0, richnessL ~ stem, S0 ~ Absolute)

)


models_stemDi <- define_model_set(
  MS3_1   = c(richnessL ~ S9, richnessL ~ stem),
  MS3_2   = c(richnessL ~ S9, richnessL ~ stem, S9 ~ Absolute),
  MS3_22   = c(richnessL ~ S9, richnessL ~ stem, S9 ~ AbsoluteDi),
  MS3_222   = c(richnessL ~ S9, richnessL ~ stem, S9 ~ AbsoluteDi, S9 ~ Absolute),
  MSAR   = c( Absolute ~richnessL),
  MSRA   = c(richnessL~Absolute ),
  MS3_3   = c(richnessL ~ S5, richnessL ~ stem),
  MS3_4   = c(richnessL ~ S5, richnessL ~ stem, S5 ~ Absolute),
  MS3_42   = c(richnessL ~ S5, richnessL ~ stem, S5 ~ AbsoluteDi),
  MS3_422   = c(richnessL ~ S5, richnessL ~ stem, S5 ~ AbsoluteDi, S5 ~ Absolute),
  MS3_5   = c(richnessL ~ S0, richnessL ~ stem),
  MS3_6   = c(richnessL ~ S0, richnessL ~ stem, S0 ~ Absolute),
  MS3_62   = c(richnessL ~ S0, richnessL ~ stem, S0 ~ AbsoluteDi),
  MS3_622   = c(richnessL ~ S0, richnessL ~ stem, S0 ~ AbsoluteDi, S0 ~ Absolute)
)


pdf("models_stem.pdf", 10,20)
plot_model_set(models_stemDi)
dev.off()


models_crown_MS <- define_model_set(
MS3_1   = c(richnessL ~ C9, richnessL ~ crown),
MS3_2   = c(richnessL ~ C9, richnessL ~ crown, C9 ~ Absolute),
MS3_22   = c(richnessL ~ C9, richnessL ~ crown, C9 ~ AbsoluteDi),
MS3_222   = c(richnessL ~ C9, richnessL ~ crown, C9 ~ AbsoluteDi, C9 ~ Absolute),
MS3_3   = c(richnessL ~ C5, richnessL ~ crown),
MS3_4   = c(richnessL ~ C5, richnessL ~ crown, C5 ~ Absolute),
MS3_42   = c(richnessL ~ C5, richnessL ~ crown, C5 ~ AbsoluteDi),
MS3_422   = c(richnessL ~ C5, richnessL ~ crown, C5 ~ AbsoluteDi, C5 ~ Absolute),
MS3_5   = c(richnessL ~ C0, richnessL ~ crown),
MS3_6   = c(richnessL ~ C0, richnessL ~ crown, C0 ~ Absolute),
MS3_62   = c(richnessL ~ C0, richnessL ~ crown, C0 ~ AbsoluteDi),
MS3_622   = c(richnessL ~ C0, richnessL ~ crown, C0 ~ AbsoluteDi, C0 ~ Absolute),
MSAR   = c( Absolute ~richnessL),
MSRA   = c(richnessL~Absolute )
)


pdf("models_crown.pdf", 10,20)
plot_model_set(models_crown_MS)
dev.off()


models_crown_DR <- define_model_set(
DR3_1   = c(richnessL ~ DR_mean, richnessL ~ crown),
DR3_2   = c(richnessL ~ DR_mean, richnessL ~ crown, DR_mean ~ Absolute),
DR3_22   = c(richnessL ~ DR_mean, richnessL ~ crown, DR_mean ~ AbsoluteDi),
DR3_222   = c(richnessL ~ DR_mean, richnessL ~ crown, DR_mean ~ AbsoluteDi, DR_mean ~ Absolute),
MSAR   = c( Absolute ~richnessL),
MSRA   = c(richnessL~Absolute )
)

pdf("models_DR.pdf", 10,10)
plot_model_set(models_crown_DR)
dev.off()



 ##Generic-level tree
 
 tips<-tree$tip.label
 genera<-unique(sapply(strsplit(tips,"_"),function(x) x[1]))
 ii<-sapply(genera,function(x,y) grep(x,y)[1],y=tips)
 ingroup_genera<-drop.tip(tree,setdiff(tree$tip.label,tips[ii])[-126]) ##Eruca was being dropped
 ingroup_genera$tip.label<-gsub("_.*","\\1",ingroup_genera$tip.label)
 
 name.check(ingroup_genera, data_ppa) ##Ones that we dropped
 
 new_ingroup_genera<-drop.tip(ingroup_genera, name.check(ingroup_genera, data_ppa)$tree_not_data)
 new_dataset<-data_ppa
 
 new_dataset$Absolute<-log( new_dataset$Absolute+ 1)
 new_dataset$AbsoluteDi<-log( (as.numeric(as.character(new_dataset$`# tips (species)`)) - new_dataset$Absolute) + 1)
 
 new_dataset$richnessL<-log( new_dataset$richness+ 1)
 
 ##Fit for stem
 result_stem <- phylo_path(models_stem, data = new_dataset, tree = new_ingroup_genera, model = 'lambda')
 result_stemDi <- phylo_path(models_stemDi, data = new_dataset, tree = new_ingroup_genera, model = 'lambda')
 
 result_crown_MS <- phylo_path(models_crown_MS, data = new_dataset, tree = new_ingroup_genera, model = 'lambda')
 result_crown_DR <- phylo_path(models_crown_DR, data = new_dataset, tree = new_ingroup_genera, model = 'lambda')
 
 write.csv(new_dataset, "Phylogenetic.path.analysis.table.csv")
 
 (S1 <- summary(result_stem))
 write.csv(S1, "Stem.csv")
 (S1Di <- summary(result_stemDi))
 
 
 (S2 <- summary(result_crown_MS))
 write.csv(S2, "CrownMS.csv")
 
 (S3 <- summary(result_crown_DR))
 write.csv(S3, "CrownDR.csv")
 
 S11 <- average(result_stem)
 plot(S11, algorithm = 'mds', curvature = 0.1) 

 S11Di <- average(result_stemDi)
 plot(S11Di, algorithm = 'mds', curvature = 0.1) 
 
 S22 <- average(result_crown_MS)
 plot(S22, algorithm = 'mds', curvature = 0.1) 
 
 S33 <- average(result_crown_DR)
 plot(S33, algorithm = 'mds', curvature = 0.1) 
 
 
 best_cichlids <- best(result_crown_DR)
 coef_plot(best_cichlids, error_bar = "se", reverse_order = TRUE) + ggplot2::coord_flip()
  plot(best_cichlids)
 
  my_model <- choice(result_stem, "MS3_6")
  coef_plot(my_model, error_bar = "se", reverse_order = TRUE) + ggplot2::coord_flip()
  plot(my_model)
 
  
  
 ##Phylogenetic anova

 DR2<-DR[-which( row.names(DR) %in% name.check(ntree, DR)$data_not_tree ),]
names(DR2)<-  row.names(DR)[-which( row.names(DR) %in% name.check(ntree, DR)$data_not_tree )]
 
 
 y1<-hdbin_first$`dupli$V5`
 names(y1)<-hdbin_first$`dupli$V1`
 
 phylANOVA(ntree,DR2,y1)
 
 phyl.pairedttest(ntree,DR2,y1)
 
 
 library(doBy)

 mean(DR2[y1==1]) 
 mean(DR2[y1==0]) 
 
 ##Phylosignal
 
library(phytools) 
x<-hdbin_first$`dupli$V5`
 names(x)<-hdbin_first$`dupli$V1`
phylosig(ntree,x,method="lambda",test=TRUE)


write.tree(new_ingroup_genera, "generic.tree")
write.tree(ntree, "hisse.tree")



##New PGLS Diversification and polyploid richness
library(caper)
CD <- comparative.data(new_ingroup_genera, new_dataset, Row.names, vcv=TRUE, vcv.dim=3, na.omit = F)

summary(pgls(C9 ~log(Absolute+1), CD, lambda='ML'))


write.csv(hdbin_first, "TableS2.csv")




###Plot Figure 1--------

pdf("hisse_plot.pdf",10,10)
hisse::plot.hisse.states(h4_rec, show.tip.label=F, rate.param="net.div",
                         state.colors=c("corn flower blue","#FF6347"),rate.colors=c("white", "black"),
                         edge.width.rate=6,
                         edge.width.state=2)
dev.off()
library(phytools)
pdf("time_scale_plot.pdf",10,10)
plotTree(ntree,type="fan",part=0.93,fsize=0.0001)
h<-max(nodeHeights(ntree))
axis(1,pos=-0.05*h,at=round(seq(0,h,by=h/5)),lwd=2)
dev.off()



##Plot
trait<-as.data.frame(data2[,2])
row.names(trait)<-data2[,1]
colnames(trait)<-c("Ploidy")
pdf("trait.plot.pdf",10,10)
trait.plot(tree2, trait, cols = list(ploidy = c("corn flower blue","#FF6347")), cex.lab = 0.3, cex.legend = 2,
           legend = T)
dev.off()


tips<-tree2$tip.label
genera2<-sapply(strsplit(tips,"_"),function(x) x[1])

pdf("trait2.plot.pdf",10,10)
trait.plot(tree2, trait, cols = list(ploidy = c("corn flower blue","#FF6347")), cex.lab = 0.3, cex.legend = 2,
           legend = T, class = genera2)
dev.off()




###Plot figure 3-------
library(ape)
library(adephylo)
library(phylobase)
library(phylosignal)

new_ingroup_genera
df3<-new_dataset
df3[is.na(df3)] <- 0
head(df3)

p4d <- phylo4d(new_ingroup_genera, df3[,c("richnessL","stem","S9", "Richness_proportion", "Absolute")])


pdf("Generic_level_tree.pdf", 10,10)
barplot(p4d, tree.type = "fan", tip.cex = 0.001, tree.open.angle = 180,
        trait.cex = 0.6, scale=F, tree.ladderize=T, center=F,
        trait.bg.col = c("#F6CED8", "#CED8F6", "#f2f6ce","#CEF6CE", "#ebcef6"), 
        bar.col = "black", bar.lwd=6, tree.ratio=0.2)
dev.off()


pdf("Generic_level_tree2.pdf", 10,10)
barplot(p4d, tip.cex = .5,
        trait.cex = 0.6, scale=F, tree.ladderize=T, center=F,
        trait.bg.col = c("#F6CED8", "#CED8F6", "#f2f6ce","#CEF6CE", "#ebcef6"), 
        bar.col = "black", bar.lwd=4, tree.ratio=0.35)
dev.off()




pdf("Generic_level_tree3.pdf", 10,30)
barplot(p4d, tip.cex = 1,
        trait.cex = 0.6, scale=F, tree.ladderize=T, center=F,
        trait.bg.col = c("#F6CED8", "#CED8F6", "#f2f6ce","#CEF6CE", "#ebcef6"), 
        bar.col = "black", bar.lwd=6, tree.ratio=0.35, tree.type = "fan")
dev.off()



###Plot bpxplots

##sub.type <- Speciation, extinction, diversification and turnover
##task.type   <- 0A, 1A, 0B, 1B
##data.value   <- val

colnames(h4_sr2$points.within.region)

data.tbl1<- rbind.fill(lapply(2:9, function(x){
  cbind.data.frame(sub.type=strsplit(colnames(h4_sr2$points.within.region)[x], ".", fixed = T)[[1]][1], 
                   task.type=strsplit(colnames(h4_sr2$points.within.region)[x], ".", fixed = T)[[1]][2],
                   data.value=h4_sr2$all.points[,x])
}))
data.tbl2<- rbind.fill(lapply(2:5, function(x){
  cbind.data.frame(sub.type="net.div", 
                   task.type=strsplit(colnames(h4_sr2$points.within.region)[x], ".", fixed = T)[[1]][2],
                   data.value=h4_sr2$all.points[,x] - h4_sr2$all.points[,x+4]   )
}))
data.tbl3<- rbind.fill(lapply(2:5, function(x){
  cbind.data.frame(sub.type="turnover", 
                   task.type=strsplit(colnames(h4_sr2$points.within.region)[x], ".", fixed = T)[[1]][2],
                   data.value=h4_sr2$all.points[,x+4]/h4_sr2$all.points[,x])
}))

data.tbl<-rbind.fill(list(data.tbl1,data.tbl2,data.tbl3))

sub.types <- as.character(unique(data.tbl$sub.type))
task.types <- as.character(unique(data.tbl$task.type))




# now have a fake dataset, so can shift to plotting. The code first sets some parameters used for both plots, then generates each plot separately.  
#windows(width=7, height=3); # specify plot window size. quartz() is the mac equivalent. the code will run without this, but the plot size (and so needed boxplot  
# widths, font sizes, etc will vary.  
#layout(matrix(c(1,1,2), c(1,3)));  # put two plots side-by-side, first twice as wide as the second  
# mar: c(bottom, left, top, right), number of lines of margin. default is c(5, 4, 4, 2) + 0.1.  
# mgp: margin line (in mex units) for the axis title, axis labels and axis line. mgp[1] is title, mgp[2:3] is axis. default is c(3, 1, 0).  
# The length of tick marks as a fraction of the height of a line of text. The default value is -0.5; setting tcl = NA sets tck = -0.01 which is S' default.  


# plotting-related variables used to simplify plotting code below. I generally fiddle with the cx, shifts, boxwex variables to make the plots readable.  
cx <- 1;  # cex for axis label text; default is 1. bigger numbers make the text bigger, smaller, smaller.  
x.ttl <- ""; # blank x-axis label  
y.ttl <- "Rate (events / My)";  
ttl <- "";  
y.lim <- c(min(data.tbl$data.value), max(data.tbl$data.value)+1.5); # y-axis limits. +2 at top to give room for legend  


# 1st: plot four boxes together: one for each task.type, grouped by sub.type  
x.lim <- c(0.5, (length(unique(data.tbl$sub.type))+0.5));  # x-axis limits; one for each of the two sub.types  
cols <- c('gold3', 'chartreuse3', 'corn flower blue', '#FF6347');  # color to plot each of the task.types (in task.types order)  
shifts <- c(-0.3, -0.1, 0.1, 0.3); # x offset for the bars (see boxplot() call below).  
# I find it easier to hard-code where the boxplots will fall (via shifts) rather than calculating the plotting locations on the fly;   

plot(x=0, y=0, xlim=x.lim, ylim=y.lim, xaxt='n', col='white', xlab=x.ttl, ylab=y.ttl, main=ttl, cex.main=cx, cex.axis=cx, cex.lab=cx);  
axis(side=1, at=1:length(unique(data.tbl$sub.type)), labels=unique(data.tbl$sub.type), cex.axis=cx, cex.lab=cx);  # put on the x-axis labels  
grid(nx=NA, ny=NULL, col='darkgrey');  
lines(x=c(-1,100), y=c(0,0), col='darkgrey');  
for (t.id in 1:length(task.types)) {  
  for (i in 1:length(sub.types)) {  # t.id <- 1; i <- 1;   
    inds <- which(data.tbl$sub.type == sub.types[i] & data.tbl$task.type == task.types[t.id]); # rows for this sub.type and task.type; what we'll plot.  
    #if (length(inds) != length(c.ids)) { stop("length(inds) != length(c.ids)"); } # error-checking code; works since length(c.ids) == length(p.ids)  
    boxplot(data.tbl$data.value[inds], at=i+shifts[t.id], col=cols[t.id], add=TRUE, xaxt='n', yaxt='n', bty='n', boxwex=0.4, cex=1); # draw the boxplot!  
  }  
}  
legend(x='top', legend=task.types, fill=cols, horiz=TRUE, cex=0.8, bg='white', bty='n'); # set font size a bit smaller to fit  
box(); # redraw the box around the outside to be nice and neat  



pdf("boxplot_ChromEvol.pdf",5,5)
library(ggplot2)
bp <- ggplot(data.tbl, aes(x=task.type, y=data.value, group=task.type))+
  geom_hline(yintercept=0, linetype="dotted", color = "black") + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(aes(fill=task.type), outlier.size = NA) + facet_wrap(.~sub.type, scales="free")
print(bp+theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
        scale_fill_manual(values=c('gold3', 'chartreuse3', 'corn flower blue', '#FF6347'))+
        xlab("")+ylab(""))
dev.off()
