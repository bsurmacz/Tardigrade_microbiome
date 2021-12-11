library(reshape2)
library(ggplot2)
	library(RColorBrewer)
	
	
samples_experiment_1<- read.csv("../EXPERIMENT1/SAMPLES_EXPERIMENT1.csv",sep="\t", header=T)
otus_experiment_1<-read.csv("../EXPERIMENT1/RENAMED/Results/Decontaminated_OTU_Table.txt",sep="\t", header=T,row.names=1)
Statistic_table_experiment_1<- read.csv("../EXPERIMENT1/RENAMED/Results/Statistics_table.txt",sep="\t",header=T,row.names=1)
otus_full_experiment_1<-read.csv("../EXPERIMENT1/RENAMED/Results/otu_table_experiment_1.txt",sep="\t", header=T,row.names=1)


samples_experiment_1<-samples_experiment_1[,!colnames(samples_experiment_1) %in% c("EXPERIMENT_1_TZ_Neg_PCRI","EXPERIMENT_1_T_AU_031_1")  ]
otus_experiment_1<-otus_experiment_1[,!colnames(otus_experiment_1) %in% c("EXPERIMENT_1_TZ_Neg_PCRI","EXPERIMENT_1_T_AU_031_1")  ]
otus_full_experiment_1<-otus_full_experiment_1[,!colnames(otus_full_experiment_1) %in% c("EXPERIMENT_1_TZ_Neg_PCRI","EXPERIMENT_1_T_AU_031_1")  ]
Statistic_table_experiment_1<-Statistic_table_experiment_1[,!colnames(Statistic_table_experiment_1) %in% c("EXPERIMENT_1_TZ_Neg_PCRI","EXPERIMENT_1_T_AU_031_1")  ]


samples_experiment_2<- read.csv("../EXPERIMENT2/SAMPLES_EXPERIMENT2.csv",sep="\t", header=T)
otus_experiment_2<-read.csv("../EXPERIMENT2/RENAMED/Results/Decontaminated_OTU_Table.txt",sep="\t", header=T,row.names=1)
otus_full_experiment_2<-read.csv("../EXPERIMENT2/RENAMED/Results/otu_table_experiment_2.txt",sep="\t", header=T,row.names=1)
Statistic_table_experiment_2<- read.csv("../EXPERIMENT2/RENAMED/Results/Statistics_table.txt",sep="\t",header=T,row.names=1)


samples_experiment_3<- read.csv("../EXPERIMENT3/SAMPLES_EXPERIMENT3.csv",sep="\t", header=T)
otus_experiment_3<-read.csv("../EXPERIMENT3/RENAMED/Results/Decontaminated_OTU_Table.txt",sep="\t", header=T,row.names=1)
otus_full_experiment_3<-read.csv("../EXPERIMENT3/RENAMED/Results/otu_table_experiment_3.txt",sep="\t", header=T,row.names=1)
Statistic_table_experiment_3<- read.csv("../EXPERIMENT3/RENAMED/Results/Statistics_table.txt",sep="\t",header=T,row.names=1)




####
### otus to show
####
OTUS_TO_SHOW<-function(otus,otus_full,N_mean_otus,N_mean_full,N_max_otus,N_max_full){
	otus_data<-otus[-c(1:4)]
	
	otus_full<-apply(otus_full, 2, function(x){100*x/max(1,sum(x))} )
	otus_data<-apply(otus_data, 2, function(x){100*x/max(1,sum(x))} )
		
	otus_MEAN<-otus_data[rev(order(rowSums(otus_data,na.rm=T))),]
	otus_MAX<-otus_data[rev(order(apply(otus_data,1,max))),]
	otus_full_MEAN<-otus_full[rev(order(rowSums(otus_full,na.rm=T))),]
	otus_full_MAX<-otus_full[rev(order(apply(otus_full,1,max))),]

	max_full<-rownames(otus_full_MAX)[1:N_max_full]
	max_decontaminated<-rownames(otus_full_MAX)[1:N_max_otus]
	mean_full<-rownames(otus_full_MEAN)[1:N_mean_full]
	mean_decontaminated<-rownames(otus_MEAN)[1:N_mean_otus]

	otus_to_show_1<-unique(tolower(c(max_full,max_decontaminated,mean_full,mean_decontaminated)))
	otus_to_show_1<-otus_to_show_1[ order( as.numeric(gsub('\\D+','', otus_to_show_1)) )  ]

	return(otus_to_show_1)
}






OTU_PLOT<-function(otus,N_col_to_remove,nonblanks_only,samples,otus_to_show){

	rownames(otus)<-tolower(rownames(otus))
#TO RELATIVE
	X<-otus
	if(N_col_to_remove>0){info<-X[1:N_col_to_remove] }
	if(N_col_to_remove>0){X<-X[-c(1:N_col_to_remove)] }
	X<-apply(X, 2, function(x){100*x/max(1,sum(x))} )
	#out<-cbind(info,X)
	otus<-X


	otus<-otus[ rownames(otus) %in% otus_to_show,  ]
	otus<-otus[ order( as.numeric(gsub('\\D+','', rownames(otus))) ) , ]
	

	otus<-rbind(100-colSums(otus),otus)
	print(otus)
	print(rownames(otus))
	rownames(otus)[1]<-"Others"

	print(otus)
	if(nonblanks_only==T){

		nothing<-as.data.frame(t(c(otus[1,]*0)))
		rownames(nothing)<-"nothing"
		nothing[,samples[ 4,colnames(nothing)]!=""  | colSums(otus)==0  ]<-100
		otus[,samples[ 4,colnames(nothing)]!=""] <-0
		otus<-rbind(otus,nothing)
		}


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

#	if(nonblanks_only==T){
#		melted_otus<-subset(melted_otus,BLANK_TYPE=="4.not_blank")
#		}



	THEME<-theme(panel.margin=unit(c(0), "lines"),  axis.ticks.x=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_rect(colour = "#ffffff", fill=NA, size=0.5), strip.background = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),strip.text   = element_text(size=0),legend.title=element_blank() )





	color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

	set.seed(113+1009+333)
	col_vector<-c("#dadada",brewer.pal(n = 12, name = "Paired"),sample(color, 408),"white")
	col_vector<-setNames(col_vector, c("Others",paste0("otu",1:420),"nothing"))
	col_vector[["otu14"]]<-"cyan3"
	col_vector[["otu72"]]<-"bisque1"
	col_vector[["otu133"]]<-"yellow2"


	output_plot<-ggplot( melted_otus, aes(x=SAMPLE_NAME,y=value,fill=Var2) )+geom_bar(stat="identity",position="stack", width=2)+ theme(panel.border = element_blank() ,panel.background = element_blank() )+facet_grid(~EXTRACTION+WITH_REPLICATES+BLANK_TYPE+GENUS+SPECIES+REPLICATE+N+SAMPLE, scales = "free", space = "free")+THEME+ scale_y_continuous(limits = c(0,101), expand = c(0, 0))+xlab("")+ylab("Relative abundance") +scale_fill_manual(values=col_vector)

	return(output_plot)
	}



###################################################
#   Contamination plots
##################################################
CONTAMINATION_PLOT<-function(Statistic_table,samples){

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

OUTPUT_PLOT<-ggplot( melted, aes(x=SAMPLE_NAME,y=value,fill=Var2) )+geom_bar(stat="identity",position="stack", width=2)+ theme(panel.border = element_blank() ,panel.background = element_blank() )+facet_grid(~EXTRACTION+WITH_REPLICATES+BLANK_TYPE+GENUS+SPECIES+REPLICATE+N+SAMPLE, scales = "free", space = "free")+THEME+ scale_y_continuous(limits = c(0,101), expand = c(0, 0))+scale_fill_manual(values=COLORS)+xlab("")+ylab("Relative abundance")
return (OUTPUT_PLOT)
}











####################################################

library(ggpubr)

#		otus,otus_full,N_mean_otus,N_mean_full,N_max_otus,N_max_full

otus_to_show_1<-OTUS_TO_SHOW( otus_experiment_1, otus_full_experiment_1,6,6,10,8)
otus_to_show_2<-OTUS_TO_SHOW( otus_experiment_2, otus_full_experiment_2,6,6,10,8)
otus_to_show_3<-OTUS_TO_SHOW( otus_experiment_3, otus_full_experiment_3,6,6,10,8)

otus_to_show_1<-c(otus_to_show_1,"otu26")
otus_to_show_2<-c(otus_to_show_2,"otu26")
otus_to_show_3<-c(otus_to_show_3,"otu26","otu17")

FULL_PLOT_experiment_1<-OTU_PLOT(otus_full_experiment_1,0,FALSE,samples_experiment_1,otus_to_show_1)
DECONTAMINATED_PLOT_experiment_1<-OTU_PLOT(otus_experiment_1,4,TRUE,samples_experiment_1,otus_to_show_1)

FULL_PLOT_experiment_2<-OTU_PLOT(otus_full_experiment_2,0,FALSE,samples_experiment_2,otus_to_show_2)
DECONTAMINATED_PLOT_experiment_2<-OTU_PLOT(otus_experiment_2,4,TRUE,samples_experiment_2,otus_to_show_2)


FULL_PLOT_experiment_3<-OTU_PLOT(otus_full_experiment_3,0,FALSE,samples_experiment_3,otus_to_show_3)
DECONTAMINATED_PLOT_experiment_3<-OTU_PLOT(otus_experiment_3,4,TRUE,samples_experiment_3,otus_to_show_3)



CONTAMINATION_PLOT_experiment_1<-CONTAMINATION_PLOT(Statistic_table_experiment_1,samples_experiment_1)
CONTAMINATION_PLOT_experiment_2<-CONTAMINATION_PLOT(Statistic_table_experiment_2,samples_experiment_2)
CONTAMINATION_PLOT_experiment_3<-CONTAMINATION_PLOT(Statistic_table_experiment_3,samples_experiment_3)


print(otus_to_show_1)
print(otus_to_show_2)
print(otus_to_show_3)

dev.new()
svg("EXPERIMENT_1_otus.svg",width=20, height=26)
ggarrange(CONTAMINATION_PLOT_experiment_1,FULL_PLOT_experiment_1,DECONTAMINATED_PLOT_experiment_1, ncol=1,common.legend = FALSE, legend="bottom")
dev.off()

dev.new()
svg("EXPERIMENT_2_otus.svg",width=20, height=26)
ggarrange(CONTAMINATION_PLOT_experiment_2,FULL_PLOT_experiment_2,DECONTAMINATED_PLOT_experiment_2, ncol=1,common.legend = FALSE, legend="bottom")
dev.off()

dev.new()
svg("EXPERIMENT_3_otus.svg",width=20, height=26)
ggarrange(CONTAMINATION_PLOT_experiment_3,FULL_PLOT_experiment_3,DECONTAMINATED_PLOT_experiment_3, ncol=1,common.legend = FALSE, legend="bottom")
dev.off()


