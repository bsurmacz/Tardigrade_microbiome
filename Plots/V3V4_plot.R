library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

V3V4<-read.table("../V3V4/OTU_PLOT.csv",sep="\t", header=T,row.names=1,dec=",")
V3V4<-V3V4[ c(nrow(V3V4), 1:(nrow(V3V4)-1) ),]
V3V4

V3V4_contamination<-read.table("../V3V4/V3V4_contamination.csv",sep="\t", header=T,row.names=1,dec=",")
V3V4_contamination<-V3V4_contamination[ nrow(V3V4_contamination):1,]


samples<-read.table("../V3V4/SAMPLES_V3V4.csv",sep="\t", header=T,dec=",")

melted_v3v4<-melt(t(V3V4))

melted_v3v4_contamination<-melt(t(V3V4_contamination))







melted_v3v4$SPECIES<-t( samples[3, paste0(melted_v3v4[,1])])
melted_v3v4$PROJECT<-t( samples[2, paste0(melted_v3v4[,1])])
melted_v3v4$SAMPLE_NAME<-t( samples[11, paste0(melted_v3v4[,1])])
melted_v3v4$SAMPLE<-melted_v3v4[,1]



melted_v3v4_contamination$SPECIES<-t( samples[3, paste0(melted_v3v4_contamination[,1])])
melted_v3v4_contamination$PROJECT<-t( samples[2, paste0(melted_v3v4_contamination[,1])])
melted_v3v4_contamination$SAMPLE_NAME<-t( samples[11, paste0(melted_v3v4_contamination[,1])])
melted_v3v4_contamination$SAMPLE<-melted_v3v4_contamination[,1]







THEME<-theme(panel.margin=unit(c(0), "lines"),  axis.ticks.x=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_rect(colour = "#ffffff", fill=NA, size=0.5), strip.background = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),strip.text   = element_text(size=0),legend.title=element_blank() )



n <- 40
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

set.seed(18)
#col_vector<-sample(color, n)

col_vector<-c("#dadada",brewer.pal(9, "Paired")) #Paired


KOLORY_V3V4<-col_vector

PLOT_OTUS <- ggplot( melted_v3v4, aes(x=SAMPLE_NAME,y=100*value,fill=Var2) )+geom_bar(stat="identity",position="stack", width=2)+ theme(panel.border = element_blank() ,panel.background = element_blank() )+facet_grid(~PROJECT+SPECIES+SAMPLE, scales = "free", space = "free")+THEME+ scale_y_continuous(limits = c(0,101), expand = c(0, 0))+xlab("")+ylab("Relative abundance")+scale_fill_manual(values=KOLORY_V3V4)+ scale_x_discrete(labels= samples[11,])+ theme(legend.position="top")


KOLORY_V3V4<-c( "orange","lightgreen", "#a2a2a2")

PLOT_COMMON<-ggplot( melted_v3v4_contamination, aes(x=SAMPLE_NAME,y=value,fill=Var2) )+geom_bar(stat="identity",position="stack", width=2)+ theme(panel.border = element_blank() ,panel.background = element_blank() )+facet_grid(~PROJECT+SPECIES+SAMPLE, scales = "free", space = "free")+THEME+ scale_y_continuous(limits = c(0,101), expand = c(0, 0))+xlab("")+ylab("Relative abundance")+scale_fill_manual(values=KOLORY_V3V4)+ scale_x_discrete(labels= samples[11,])+ theme(legend.position="top")

ggarrange( PLOT_COMMON,PLOT_OTUS, ncol=1)


dev.new()
svg("V3V4.svg",width=8,height=14)
ggarrange( PLOT_COMMON,PLOT_OTUS, ncol=1)
dev.off()
dev.off()






