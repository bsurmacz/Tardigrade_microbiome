library(reshape2)
library(ggplot2)

samples<- read.csv("../EXPERIMENT1/SAMPLES_EXPERIMENT1.csv",sep="\t", header=T)
otus<-read.csv("../EXPERIMENT1/Results/Decontaminated_OTU_Table.txt",sep="\t", header=T,row.names=1)


#TO RELATIVE
	X<-otus
	info<-X[1:5]
	n<-ncol(X)
	X<-X[-c(1:5,n)]
	X<-apply(X, 2, function(x){100*x/sum(x)} )
	out<-cbind(info,X)
otus_relative<-out


data<-t(otus_relative[-c(1:5)] )

types<-data.frame(sample_type=t(samples[5, rownames(data)]))
colnames(types)<-"sample_type"

data<-cbind(data, types)

Tardigrade_samples<-subset(data,sample_type=="1")
Medium_samples<-subset(data,sample_type=="medium")

Tardigrade_samples<-subset(Tardigrade_samples, select=-sample_type) 
Medium_samples<-subset(Medium_samples, select=-sample_type) 

Tardigrade_samples<-apply(Tardigrade_samples,2,as.numeric)
Medium_samples<-apply(Medium_samples,2,as.numeric)

Tardigrade_means<-apply(Tardigrade_samples,2,max)
Medium_means<-apply(Medium_samples,2,max)

plot(Tardigrade_means,Medium_means)








as.data.frame(t(samples[5, rownames(data)]),colnames=c("sample_type")) 

 data.frame(t(samples[5, rownames(data)]),colnames=c("sample_type"))






# EMPTY 
otus<-otus[,! colnames(otus) %in%  c("EXPERIMENT_1_T_Neg_PCRI_1","EXPERIMENT_1_T_Neg_PCRI_1","EXPERIMENT_1_T_Neg_PCRI_2")]


#TO RELATIVE
	X<-otus
	info<-X[1:5]
	n<-ncol(X)
	X<-X[-c(1:5,n)]
	X<-apply(X, 2, function(x){100*x/sum(x)} )
	out<-cbind(info,X)
otus<-out

otus<-otus[,-c(1,2,3,4,5)]

apply(otus,1,max)

otusMEAN<-otus[rev(order(rowSums(otus,na.rm=T))),]
otusMAX<-otus[rev(order(apply(otus,1,max))),]

otus<-unique(rbind(otusMEAN[1:5,], otusMAX[1:9,] ))



otus<-rbind(100-colSums(otus),otus)
rownames(otus)[1]<-"Others"

melted_otus<-melt(t(otus))

melted_otus$REPLICATE<-t( samples[1, paste0(melted_otus[,1])])
melted_otus$EXTRACTION<-t( samples[2, paste0(melted_otus[,1])])

melted_otus$EXTRACTION[melted_otus$EXTRACTION=="Chelex"]<-"1.Chelex"
melted_otus$EXTRACTION[melted_otus$EXTRACTION=="beads"]<-"2.SPRIbeads"


melted_otus$SPECIES<-t( samples[3, paste0(melted_otus[,1])])
melted_otus$BLANK_TYPE<-t( samples[4, paste0(melted_otus[,1])])
melted_otus$N<-t( samples[5, paste0(melted_otus[,1])])
melted_otus$WASHED<-t( samples[6, paste0(melted_otus[,1])])
melted_otus$EXPERIMENT<-t( samples[9, paste0(melted_otus[,1])])
melted_otus$GENUS<-t( samples[10, paste0(melted_otus[,1])])
melted_otus$SAMPLE_NAME<-t( samples[11, paste0(melted_otus[,1])])
melted_otus$BLANK_TYPE[melted_otus$BLANK_TYPE=="blank"]<-"0.blank"
melted_otus$BLANK_TYPE[melted_otus$BLANK_TYPE=="ludwik"]<-"1.ludwik"
melted_otus$BLANK_TYPE[melted_otus$BLANK_TYPE=="zywiec"]<-"1.zywiec"
melted_otus$BLANK_TYPE[melted_otus$BLANK_TYPE=="algae"]<-"2.algae"
melted_otus$BLANK_TYPE[melted_otus$BLANK_TYPE=="rotifers"]<-"3.algae"
melted_otus$BLANK_TYPE[melted_otus$BLANK_TYPE==""]<-"4.not_blank"
melted_otus$SAMPLE<-melted_otus[,1]
melted_otus$WITH_REPLICATES<-melted_otus$REPLICATE %in% c("A","B","C")
melted_otus$NAME<-paste(melted_otus$EXTRACTION, melted_otus$SPECIES,melted_otus$REPLICATE,melted_otus$N)

melted_otus<-subset(melted_otus,BLANK_TYPE=="4.not_blank")


THEME<-theme(panel.margin=unit(c(0), "lines"),  axis.ticks.x=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_rect(colour = "#ffffff", fill=NA, size=0.5), strip.background = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),strip.text   = element_text(size=0),legend.title=element_blank() )


library(RColorBrewer)
col_vector<-c("#dadada",brewer.pal(11, "Paired") )


dev.new()
svg("EXPERIMENT_1_decontaminated_otus.svg",width=20,height=9)
ggplot( melted_otus, aes(x=SAMPLE_NAME,y=value,fill=Var2) )+geom_bar(stat="identity",position="stack", width=2)+ theme(panel.border = element_blank() ,panel.background = element_blank() )+facet_grid(~EXTRACTION+WITH_REPLICATES+BLANK_TYPE+GENUS+SPECIES+REPLICATE+N+SAMPLE, scales = "free", space = "free")+THEME+ scale_y_continuous(limits = c(0,101), expand = c(0, 0))+xlab("")+ylab("Relative abundance") +scale_fill_manual(values=col_vector)

dev.off()


