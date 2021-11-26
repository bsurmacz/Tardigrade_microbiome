library(reshape2)
library(ggplot2)

samples<- read.csv("../EXPERIMENT1/SAMPLES_EXPERIMENT1.csv",sep="\t", header=T)
Statistic_table<- read.csv("../EXPERIMENT1/Results/Statistics_table.txt",sep="\t",header=T,row.names=1)

percentages<-Statistic_table[c(12,22,21,13,14,17),  ] 
melted<-melt(t(percentages))

melted$REPLICATE<-t( samples[1, paste0(melted[,1])])


melted$EXTRACTION<-t( samples[2, paste0(melted[,1])])
melted$EXTRACTION[melted$EXTRACTION=="Chelex"]<-"1.Chelex"
melted$EXTRACTION[melted$EXTRACTION=="beads"]<-"2.SPRIbeads"
melted$SPECIES<-t( samples[3, paste0(melted[,1])])
melted$BLANK_TYPE<-t( samples[4, paste0(melted[,1])])
melted$N<-t( samples[5, paste0(melted[,1])])
melted$WASHED<-t( samples[6, paste0(melted[,1])])
melted$EXPERIMENT<-t( samples[9, paste0(melted[,1])])
melted$GENUS<-t( samples[10, paste0(melted[,1])])
melted$SAMPLE_NAME<-t( samples[11, paste0(melted[,1])])
melted$BLANK_TYPE[melted$BLANK_TYPE=="blank"]<-"0.blank"
melted$BLANK_TYPE[melted$BLANK_TYPE=="ludwik"]<-"1.ludwik"
melted$BLANK_TYPE[melted$BLANK_TYPE=="zywiec"]<-"1.zywiec"
melted$BLANK_TYPE[melted$BLANK_TYPE=="algae"]<-"2.algae"
melted$BLANK_TYPE[melted$BLANK_TYPE=="rotifers"]<-"3.algae"
melted$BLANK_TYPE[melted$BLANK_TYPE==""]<-"4.not_blank"
melted$SAMPLE<-melted[,1]
melted$WITH_REPLICATES<-melted$REPLICATE %in% c("A","B","C")
melted$NAME<-paste(melted$EXTRACTION, melted$SPECIES,melted$REPLICATE,melted$N)



COLORS<-c("darkolivegreen1","darkorange3","chartreuse4","#777777","#555555","#ff00ff")

THEME<-theme(panel.margin=unit(c(0), "lines"),  axis.ticks.x=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_rect(colour = "#ffffff", fill=NA, size=0.5), strip.background = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),strip.text   = element_text(size=0),legend.title=element_blank() )



dev.new()
svg("EXPERIMENT_1_contamination_plot.svg",width=20,height=9)

ggplot( melted, aes(x=SAMPLE_NAME,y=value,fill=Var2) )+geom_bar(stat="identity",position="stack", width=2)+ theme(panel.border = element_blank() ,panel.background = element_blank() )+facet_grid(~EXTRACTION+WITH_REPLICATES+BLANK_TYPE+GENUS+SPECIES+REPLICATE+N+SAMPLE, scales = "free", space = "free")+THEME+ scale_y_continuous(limits = c(0,101), expand = c(0, 0))+scale_fill_manual(values=COLORS)+xlab("")+ylab("Relative abundance")
dev.off()


