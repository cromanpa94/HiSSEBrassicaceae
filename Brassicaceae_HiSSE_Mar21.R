library(ape)
library(geiger)
library(hisse)
library(diversitree)
library(phytools)
library(phangorn)

###Running HiSSE
setwd("/home/u6/cromanpa94/HiSSEBrassicaceae")

##Sampling fractions
sfs<-list('Hohmann'=c(1-.433,.433), 'Wood'=c(1-0.3358,0.3358), 'dataset_specific'=NULL)

##Function
fit_SSE_Brassicaceae<-function(tree=tree2, data=data2, f=NULL, n.cores=28){

  tree2<-tree; ntree<-tree2
  nnls<-nnls.tree(cophenetic(ntree),ntree,rooted=TRUE)
  
  data2<-data; hdbin_first<-data2  
  
  
  if(is.null(f)){
    h<-length(which(hdbin_first[,2] == 0))
    c<-length(which(hdbin_first[,2] == 1))
    sap<-c(h,c)/length(hdbin_first[,2])
  }else{ sap<- f}
  
  
  
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
  
  ##Diploidizations set to 0
  trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
  trans.rates.hisse <- ParDrop(trans.rates.hisse, c(1,8,3,9)) 
  h4_ND <- hisse(ntree, hdbin_first, f=sap ,hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4),
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

  cid2 = hisse(ntree, hdbin_first, f=sap, hidden.states=TRUE, turnover.anc=c(1,1,2,2),
                 eps.anc=c(1,1,2,2), trans.rate=trans.rates.nodual)
  
  #hisse CID-4
  cid4 <- hisse.null4(ntree, hdbin_first, f=sap, turnover.anc=rep(c(1,2,3,4),2), 
                               eps.anc=rep(c(1,2,3,4),2))
  
  modelnames<-c(paste0("BiSSE_", 1:4), paste0("HiSSE_", 1:3), "HiSSE_4","HiSSE_4Ndiploidi", paste0("HiSSE_", 5:6) , "HiSSE-CID2","HiSSE-CID4")
  
  models<-list(b1,b2,b3,b4,h1,h2,h3,h4,h4_ND,h5,h6, cid2,cid4)
  
  Table<- cbind.data.frame(Model=modelnames, 
        lnLik=c(b1$loglik,b2$loglik,b3$loglik,b4$loglik,
                h1$loglik,h2$loglik,h3$loglik,h4$loglik,h4_ND$loglik,h5$loglik,h6$loglik,
                cid2$loglik,cid4$loglik),
  AIC=c(b1$AIC,b2$AIC,b3$AIC,b4$AIC,h1$AIC,h2$AIC,h3$AIC,h4$AIC,h4_ND$AIC,h5$AIC,h6$AIC,cid2$AIC,
        cid4$AIC))

 Table$deltaAIC<- akaike.weights(Table$AIC)$deltaAIC
 Table$weights<- akaike.weights(Table$AIC)$weights

 best_model<-models[[which.min(Table$deltaAIC)]]
 
 bestrec <- MarginRecon(ntree, hdbin_first, f=sap, hidden.states=TRUE, pars=best_model$solution, n.cores=n.cores)
 
 res<-list('sap'=sap, 'Models'=models, 'Data'= Table, 'Best'=best_model, 'ASR'=bestrec)
 
 return(res)
}

##Chromevol
load("Script_ploidy_Feb14-ChromEvol.RData")
tree2
data2

CE<-lapply(seq_along(sfs), function(x){
  fit_SSE_Brassicaceae(tree=tree2, data=data2, f=sfs[[x]])
})
names(CE)<-names(sfs)

##Stebbins
load("Script_ploidy_Feb13-Stebbins.RData")
tree2
hdbin_first<-as.matrix(data_new[,c("Species","Ploidy")] ) 

SF<-lapply(seq_along(sfs), function(x){
  fit_SSE_Brassicaceae(tree=tree2, data=hdbin_first, f=sfs[[x]])
})
names(SF)<-names(sfs)

