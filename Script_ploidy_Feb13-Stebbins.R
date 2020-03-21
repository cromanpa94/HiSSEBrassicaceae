library(ape)
library(geiger)
library(hisse)
library(qpcR)


setwd("~/Desktop/Brassicaceae/ClusterJan12")
load("~/Desktop/Brassicaceae/ClusterJan12/Script_ploidy_Jan28.RData")

Species_Ploidy_Final_Class <- read.csv("~/Desktop/Plants_CCDB/Counts/Species_Ploidy_Final_Class.csv")
row.names(Species_Ploidy_Final_Class)<-Species_Ploidy_Final_Class$Species

data<-Species_Ploidy_Final_Class[which(Species_Ploidy_Final_Class$family == "Brassicaceae"),]


data_new<-data[-which(data$Species %in% name.check(tree2, data)$data_not_tree),]


spp_drop<-name.check(tree2, Species_Ploidy_Final_Class)$tree_not_data

tree2<-drop.tip(tree2,spp_drop)  


##Run HiSSE and BiSSE

ntree<-tree2

hdbin_first<-as.matrix(data_new[,c("Species","Ploidy")] ) 

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


##Marginal reconstruction
h4_rec <- MarginRecon(ntree, hdbin_first, f=sap, hidden.states=TRUE, pars=h4$solution,
                      aic=h4$aic, n.cores=4)

##Support region
h4_sr <- SupportRegion(h4)


mr<-GetModelAveRates(h4_rec, type = "tips")

mean(mr[which(mr$state ==0),"net.div"])
mean(mr[which(mr$state ==1),"net.div"])

mean(mr[which(mr$state ==0),"speciation"])
mean(mr[which(mr$state ==1),"speciation"])

mean(mr[which(mr$state ==0),"extinction"])
mean(mr[which(mr$state ==1),"extinction"])

mean(mr[which(mr$state ==0),"turnover"])
mean(mr[which(mr$state ==1),"turnover"])
