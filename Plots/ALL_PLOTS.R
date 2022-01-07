library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(gplots)

THEME<-theme(panel.margin=unit(c(0), "lines"),  axis.ticks.x=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_rect(colour = "#ffffff", fill=NA, size=0.5), strip.background = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),strip.text   = element_text(size=0),legend.title=element_blank() )

samples_experiment_2<- read.csv("../EXPERIMENT2/SAMPLES_EXPERIMENT2.csv",sep="\t", header=T)

samples<-samples_experiment_2
otus_in<-otus_experiment_2


HEATMAP<-function(samples,otus_in)
	{
	N_col_to_remove=4
	X<-otus_in
	if(N_col_to_remove>0){info<-X[1:N_col_to_remove] }
	if(N_col_to_remove>0){X<-X[-c(1:N_col_to_remove)] }
	X<-apply(X, 2, function(x){100*x/max(1,sum(x))} )
	#out<-cbind(info,X)
	otus<-X
	otus<-otus[, samples[4,colnames(otus)]==""]
	samples_ch<-apply(samples,2,as.character)
	rownames(samples_ch)<- c( "REPLICATE", "EXTRACTION","SPECIES","BLANK_TYPE", "N","WASHED","lab","xx","experiment","SAMPLE","NAME"   )
	with_replicates<-samples_ch["REPLICATE",colnames(otus)] %in% c("A","B","C")

	ORDER<-order(samples_ch["EXTRACTION",colnames(otus)], t(with_replicates),
	samples_ch["BLANK_TYPE",colnames(otus)],
	samples_ch["SPECIES",colnames(otus)],
	samples_ch["REPLICATE",colnames(otus)],
	samples_ch["N",colnames(otus)],
	samples_ch["SAMPLE",colnames(otus)],
	samples_ch["NAME",colnames(otus)]  )

	otus<-otus[ ,ORDER ]
	otus<-otus[ rev(order(rowSums(otus)) ), ]
	otus<-otus[1:20,]
	otus<-otus[order(as.numeric(gsub("[^0-9.-]", "", rownames(otus)))),]

	#kolory<-colorRampPalette(brewer.pal(9,"Blues"))(100)
	kolory<- colorRampPalette(c("#ffffff","#0000a0","#000090","#000060"))(n = 100)
	
	HEATMAP<-heatmap.2(as.matrix(otus),Rowv=F,Colv=F,dendrogram="none",trace="none",density.info="none",col=kolory,margins=c(14,6),key=T,scale="none",breaks = seq(0,100, length.out = 101) )
	return(HEATMAP)
	}



		
samples_experiment_1<- read.csv("../EXPERIMENT1/SAMPLES_EXPERIMENT1.csv",sep="\t", header=T)
otus_experiment_1<-read.csv("../EXPERIMENT1/RENAMED/Results/Decontaminated_OTU_Table.txt",sep="\t", header=T,row.names=1)
Statistic_table_experiment_1<- read.csv("../EXPERIMENT1/RENAMED/Results/Statistics_table.txt",sep="\t",header=T,row.names=1)
otus_full_experiment_1<-read.csv("../EXPERIMENT1/RENAMED/Results/otu_table_experiment_1.txt",sep="\t", header=T,row.names=1)


excluded_experiment1<-c("EXPERIMENT_1_TZ_Neg_PCRI","EXPERIMENT_1_T_AU_031_1","EXPERIMENT_1_T_AU_031_2")
samples_experiment_1<-samples_experiment_1[,!colnames(samples_experiment_1) %in% excluded_experiment1 ]
otus_experiment_1<-otus_experiment_1[,!colnames(otus_experiment_1) %in% excluded_experiment1 ]
otus_full_experiment_1<-otus_full_experiment_1[,!colnames(otus_full_experiment_1) %in% excluded_experiment1 ]
Statistic_table_experiment_1<-Statistic_table_experiment_1[,!colnames(Statistic_table_experiment_1) %in% excluded_experiment1 ]


samples_experiment_2<- read.csv("../EXPERIMENT2/SAMPLES_EXPERIMENT2.csv",sep="\t", header=T)
otus_experiment_2<-read.csv("../EXPERIMENT2/RENAMED/Results/Decontaminated_OTU_Table.txt",sep="\t", header=T,row.names=1)
otus_full_experiment_2<-read.csv("../EXPERIMENT2/RENAMED/Results/otu_table_experiment_2.txt",sep="\t", header=T,row.names=1)
Statistic_table_experiment_2<- read.csv("../EXPERIMENT2/RENAMED/Results/Statistics_table.txt",sep="\t",header=T,row.names=1)


samples_experiment_3<- read.csv("../EXPERIMENT3/SAMPLES_EXPERIMENT3.csv",sep="\t", header=T)
otus_experiment_3<-read.csv("../EXPERIMENT3/RENAMED/Results/Decontaminated_OTU_Table.txt",sep="\t", header=T,row.names=1)
otus_full_experiment_3<-read.csv("../EXPERIMENT3/RENAMED/Results/otu_table_experiment_3.txt",sep="\t", header=T,row.names=1)
Statistic_table_experiment_3<- read.csv("../EXPERIMENT3/RENAMED/Results/Statistics_table.txt",sep="\t",header=T,row.names=1)



############################# QUANTIFICATION PLOTS

QUANTIFICATION_PLOT<-function(samples,Statistic_table){
	number_of_spikein_copies<- 1000
	quant <- t(Statistic_table["Symbionts_reads",]/(Statistic_table["Extraction_Spikeins_reads",]))
	colnames(quant)<-"quant"
	rownames(samples)<- c( "REPLICATE", "EXTRACTION","SPECIES","BLANK_TYPE", "N","WASHED","lab","xx","experiment","SAMPLE","NAME"   )
	samples<-t(samples)
	samples[,"SAMPLE"]<- rownames(samples)
	q<-as.data.frame(cbind(samples,quant))
	q$NAME <- paste0(q$NAME,q$SAMPLE)
	q$WITH_REPLICATES<-q$REPLICATE %in% c("A","B","C")
	q$BLANK_TYPE<-as.character(q$BLANK_TYPE)
	q$BLANK_TYPE[q$BLANK_TYPE=="blank"]<-"0.blank"
	q$BLANK_TYPE[q$BLANK_TYPE=="ludwik"]<-"1.ludwik"
	q$BLANK_TYPE[q$BLANK_TYPE=="zywiec"]<-"1.zywiec"
	q$BLANK_TYPE[q$BLANK_TYPE=="algae"]<-"2.algae"
	q$BLANK_TYPE[q$BLANK_TYPE=="rotifers"]<-"3.algae"
	q$BLANK_TYPE[q$BLANK_TYPE==""]<-"4.not_blank"
	q$quant<-as.numeric(as.character(q$quant))
	q$quant[q$quant==Inf]<-4000
	q$quant[q$quant>3000]<-0
	q[q$EXTRACTION=="beads",]$quant <-  5* q[q$EXTRACTION=="beads",]$quant  # in SPRIbeads extraction - spikein was added to 0.2 of sample
	q$EXTRACTION<-q$EXTRACTION*number_of_spikein_copies
	PLOT<-ggplot( q,aes(x=NAME, y=quant ))+geom_bar(stat="identity")+facet_grid(~EXTRACTION+WITH_REPLICATES+BLANK_TYPE+SPECIES+REPLICATE+N+SAMPLE+NAME, scales = "free", space = "free")+THEME
	return(PLOT)
}


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







	color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

	set.seed(113+1009+333)
	col_vector<-c("#dadada",brewer.pal(n = 12, name = "Paired"),sample(color, 408),"white")
	col_vector<-setNames(col_vector, c("Others",paste0("otu",1:420),"nothing"))
	col_vector[["otu14"]]<-"cyan3"
	col_vector[["otu18"]]<-"sienna1"
		col_vector[["otu25"]]<-"palevioletred"
	col_vector[["otu72"]]<-"bisque1"
	col_vector[["otu88"]]<-"honeydew2"
	col_vector[["otu133"]]<-"yellow2"
	col_vector[["otu1273"]]<-"aquamarine2"

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


#		otus,otus_full,N_mean_otus,N_mean_full,N_max_otus,N_max_full

otus_to_show_1<-OTUS_TO_SHOW( otus_experiment_1, otus_full_experiment_1,6,6,10,8)
otus_to_show_2<-OTUS_TO_SHOW( otus_experiment_2, otus_full_experiment_2,6,6,10,8)
otus_to_show_3<-OTUS_TO_SHOW( otus_experiment_3, otus_full_experiment_3,6,6,10,8)

otus_to_show_1<-c(otus_to_show_1,"otu26","otu1","otu4","otu7")
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

Q_PLOT_2<-QUANTIFICATION_PLOT(samples_experiment_2, Statistic_table_experiment_2)
Q_PLOT_3<-QUANTIFICATION_PLOT(samples_experiment_3, Statistic_table_experiment_3)

HEATMAP_EXP2<-HEATMAP(samples_experiment_2,otus_experiment_2)
HEATMAP_EXP3<-HEATMAP(samples_experiment_3,otus_experiment_3)


dev.new()
svg("HEATMAP_EXP2.svg",width=40, height=18)
HEATMAP(samples_experiment_2,otus_experiment_2)
dev.off()

dev.new()
svg("HEATMAP_EXP3.svg",width=40, height=18)
HEATMAP(samples_experiment_3,otus_experiment_3)
dev.off()


dev.new()
svg("EXPERIMENT_1_otus.svg",width=20, height=26)
ggarrange(CONTAMINATION_PLOT_experiment_1,FULL_PLOT_experiment_1,DECONTAMINATED_PLOT_experiment_1, ncol=1,common.legend = FALSE, legend="bottom")
dev.off()



dev.new()
svg("EXPERIMENT_2_otus.svg",width=20, height=30)
ggarrange(CONTAMINATION_PLOT_experiment_2,FULL_PLOT_experiment_2,DECONTAMINATED_PLOT_experiment_2,Q_PLOT_2, ncol=1,common.legend = FALSE, legend="bottom")
dev.off()

dev.new()
svg("EXPERIMENT_3_otus.svg",width=20, height=30)
ggarrange(CONTAMINATION_PLOT_experiment_3,FULL_PLOT_experiment_3,DECONTAMINATED_PLOT_experiment_3,Q_PLOT_3, ncol=1,common.legend = FALSE, legend="bottom")
dev.off()


