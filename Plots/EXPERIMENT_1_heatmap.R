library(reshape2)
library(ggplot2)
library(gplots)

samples<- read.csv("../EXPERIMENT1/SAMPLES_EXPERIMENT1.csv",sep="\t", header=T)
otus<- read.csv("../EXPERIMENT1/Results/otu_table_expanded.txt",sep="\t",header=T,row.names=1)


# EMPTY 
#otus<-otus[,! colnames(otus) %in%  c("EXPERIMENT_1_T_Neg_PCRI_1","EXPERIMENT_1_T_Neg_PCRI_1","EXPERIMENT_1_T_Neg_PCRI_2")]


#TO RELATIVE
	X<-otus
	info<-X[1:3]
	n<-ncol(X)
	X<-X[-c(1:3,n)]
	X<-X[order(as.numeric(gsub("[^0-9.-]", "", rownames(info)))),]
	X<-apply(X, 2, function(x){x/sum(x)} )
	out<-cbind(info,X)
otus<-X



otus<-otus[1:10,]

#otusMEAN<-otus[rev(order(rowSums(otus,na.rm=T))),]
#otusMAX<-otus[rev(order(apply(otus,1,max))),]
#otus<-unique(rbind(otusMEAN[1:5,], otusMAX[1:25,] ))

library(RColorBrewer)
kolory<-colorRampPalette(brewer.pal(9,"Blues"))(100)
    
#kolory<- colorRampPalette(c("black","red","yellow","green"))(n = 100)

svg("EXPERIMENT_1_heatmap.svg",width=40,height=20)
heatmap.2(as.matrix(otus),Rowv=F,Colv=F,dendrogram="none",trace="none",density.info="none",col=kolory,margins=c(14,6),key=T,scale="none",breaks = seq(0,1, length.out = 101) )
dev.off()

