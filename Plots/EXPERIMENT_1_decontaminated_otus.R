library(reshape2)
library(ggplot2)

samples_experiment_1<- read.csv("../EXPERIMENT1/SAMPLES_EXPERIMENT1.csv",sep="\t", header=T)
otus_experiment_1<-read.csv("../EXPERIMENT1/RENAMED/Results/Decontaminated_OTU_Table.txt",sep="\t", header=T,row.names=1)



otus_full_experiment_1<-read.csv("../EXPERIMENT1/RENAMED/Results/otu_table_experiment_1.txt",sep="\t", header=T,row.names=1)



samples_experiment_2<- read.csv("../EXPERIMENT2/SAMPLES_EXPERIMENT2.csv",sep="\t", header=T)
otus_experiment_2<-read.csv("../EXPERIMENT2/RENAMED/Results/Decontaminated_OTU_Table.txt",sep="\t", header=T,row.names=1)
otus_full_experiment_2<-read.csv("../EXPERIMENT2/RENAMED/Results/otu_table_experiment_2.txt",sep="\t", header=T,row.names=1)




samples_experiment_3<- read.csv("../EXPERIMENT3/SAMPLES_EXPERIMENT3.csv",sep="\t", header=T)
otus_experiment_3<-read.csv("../EXPERIMENT3/RENAMED/Results/Decontaminated_OTU_Table.txt",sep="\t", header=T,row.names=1)
otus_full_experiment_3<-read.csv("../EXPERIMENT3/RENAMED/Results/otu_table_experiment_3.txt",sep="\t", header=T,row.names=1)








OTU_PLOT<-function(otus,N_col_to_remove,nonblanks_only,samples){


rownames(otus)<-tolower(rownames(otus))





#TO RELATIVE
	X<-otus
	info<-X[1:N_col_to_remove]
	n<-ncol(X)
	X<-X[-c(1:N_col_to_remove)]
	X<-X[,colSums(X)>0]  #removing empty libraries
	X<-apply(X, 2, function(x){100*x/sum(x)} )
	out<-cbind(info,X)
otus<-out

otus<-otus[,-c(1:4)]

otus<-otus[ order( as.numeric(gsub('\\D+','', rownames(otus))) ) , ]
	#otusMEAN<-otus[rev(order(rowSums(otus,na.rm=T))),]
	#otusMAX<-otus[rev(order(apply(otus,1,max))),]
	#otus<-unique(rbind(otusMEAN[1:5,], otusMAX[1:9,] ))
otus<-otus[1:13, ]



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

if(nonblanks_only==T){
	melted_otus<-subset(melted_otus,BLANK_TYPE=="4.not_blank")
	}

THEME<-theme(panel.margin=unit(c(0), "lines"),  axis.ticks.x=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_rect(colour = "#ffffff", fill=NA, size=0.5), strip.background = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),strip.text   = element_text(size=0),legend.title=element_blank() )


library(RColorBrewer)
col_vector<-c("#dadada",brewer.pal(12, "Paired"),"pink" )


color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

set.seed(113+10003)
col_vector<-c("#dadada",sample(color, 30))
col_vector<-setNames(col_vector, c("Others",paste0("otu",1:30)))


output_plot<-ggplot( melted_otus, aes(x=SAMPLE_NAME,y=value,fill=Var2) )+geom_bar(stat="identity",position="stack", width=2)+ theme(panel.border = element_blank() ,panel.background = element_blank() )+facet_grid(~EXTRACTION+WITH_REPLICATES+BLANK_TYPE+GENUS+SPECIES+REPLICATE+N+SAMPLE, scales = "free", space = "free")+THEME+ scale_y_continuous(limits = c(0,101), expand = c(0, 0))+xlab("")+ylab("Relative abundance") +scale_fill_manual(values=col_vector)


	return(output_plot)
	}
library(ggpubr)




FULL_PLOT_experiment_1<-OTU_PLOT(otus_full_experiment_1,0,FALSE,samples_experiment_1)
DECONTAMINATED_PLOT_experiment_1<-OTU_PLOT(otus_experiment_1,4,TRUE,samples_experiment_1)

FULL_PLOT_experiment_2<-OTU_PLOT(otus_full_experiment_2,0,FALSE,samples_experiment_2)
DECONTAMINATED_PLOT_experiment_2<-OTU_PLOT(otus_experiment_2,4,TRUE,samples_experiment_2)

FULL_PLOT_experiment_3<-OTU_PLOT(otus_full_experiment_3,0,FALSE,samples_experiment_3)
DECONTAMINATED_PLOT_experiment_3<-OTU_PLOT(otus_experiment_3,4,TRUE,samples_experiment_3)







dev.new()
svg("EXPERIMENT_1_otus.svg",width=20, height=20)
ggarrange(FULL_PLOT_experiment_1,DECONTAMINATED_PLOT_experiment_1, ncol=1,common.legend = FALSE, legend="right")
dev.off()


dev.new()
svg("EXPERIMENT_2_otus.svg",width=20, height=20)
ggarrange(FULL_PLOT_experiment_2,DECONTAMINATED_PLOT_experiment_2, ncol=1,common.legend = FALSE, legend="right")
dev.off()


dev.new()
svg("EXPERIMENT_3_otus.svg",width=20, height=20)
ggarrange(FULL_PLOT_experiment_3,DECONTAMINATED_PLOT_experiment_3, ncol=1,common.legend = FALSE, legend="right")
dev.off()


