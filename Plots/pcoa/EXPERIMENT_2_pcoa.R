library(ggpubr)

library("vegan")
library("ggplot2")


samples<- read.csv("../../EXPERIMENT2/SAMPLES_EXPERIMENT2.csv",sep="\t", header=T)
full_data<-read.table("../../EXPERIMENT2/Results/zotu_table_expanded.txt",sep="\t",header=T)
decontaminated_data<-read.table("../../EXPERIMENT2/Results/Decontaminated_zOTU_table.txt",sep="\t",header=T)


decontaminated_data<-decontaminated_data[-ncol(decontaminated_data)]
#
full_data<-full_data[-ncol(full_data)]

decontaminated_data<-decontaminated_data[-c(1:5)]
full_data<-full_data[-c(1:5)]


relative_decontaminated<- apply(decontaminated_data, 2, function(x){x/sum(x)} )
relative_full<- apply(full_data, 2, function(x){x/sum(x)} )


EXPERIMENT_decontaminated<-t( samples[9, paste0(colnames(decontaminated_data))])
SAMPLE_TYPE<-t( samples[5, paste0(colnames(decontaminated_data))])
SPECIES_decontaminated<-t( samples[3, paste0(colnames(decontaminated_data))])


EXPERIMENT_full<-t( samples[9, paste0(colnames(full_data))])
SAMPLE_TYPE_full<-t( samples[5, paste0(colnames(full_data))])
SPECIES_full<-t( samples[3, paste0(colnames(full_data))])

EXP2_decontaminated<-relative_decontaminated[, EXPERIMENT_decontaminated=="exp3" & !(SAMPLE_TYPE  %in% c( "blank","rotifers","algae","ludwik","zywiec"))]
EXP2_full<-relative_full[, EXPERIMENT_full=="exp3"]

EXP2_decontaminated<-t(EXP2_decontaminated)
EXP2_full<-t(EXP2_full)


max<-1000
EXP2_MDS_decontaminated<-metaMDS( EXP2_decontaminated, distance="bray", k=2, trymax=max)
EXP2_MDS_full<-metaMDS( EXP2_full, distance="bray", k=2, trymax=max)

EXP2_scores_decontaminated <- as.data.frame(scores(EXP2_MDS_decontaminated)) 
EXP2_scores_full <- as.data.frame(scores(EXP2_MDS_full)) 


#######333   KOLORY
n <- 52
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
set.seed(117)
col_vector<-sample(color, n)
col_vector<-c("orange","green","#ca8b44",col_vector)

#########

EXP2_scores_decontaminated$species<-t( samples[3, paste0(rownames(EXP2_scores_decontaminated))])
EXP2_scores_decontaminated$extraction<-t( samples[2, paste0(rownames(EXP2_scores_decontaminated))])
EXP2_scores_decontaminated$type<-t( samples[5, paste0(rownames(EXP2_scores_decontaminated))])

EXP2_scores_full$species<-t( samples[3, paste0(rownames(EXP2_scores_full))])
EXP2_scores_full$extraction<-t( samples[2, paste0(rownames(EXP2_scores_full))])
EXP2_scores_full$type<-t( samples[5, paste0(rownames(EXP2_scores_full))])



library(dplyr)
 library(ggnewscale)

#by species
pcoa_plot<-function(data,Title,scale_shape,colors,extraction){
	hulls<-data %>%group_by(species) %>% slice(chull(NMDS1,NMDS2))

	plot_data<-ggplot( aes(x=NMDS1,y=NMDS2,colour=species,shape=type),data=data)+ geom_point(data=data,size=6)+ geom_polygon(data =hulls, aes(group=species,fill=species,color=species), alpha = 0.1)+scale_fill_manual(values=colors)+scale_color_manual(values=colors)+theme_bw()+scale_shape_manual(values=scale_shape)+labs(title=Title)
	if(extraction){
	plot_data<-plot_data+new_scale_fill() +geom_point(shape = 21, aes(x=NMDS1,y=NMDS2,fill=extraction),size=9,alpha=0.08)+scale_fill_manual(values=c("red","blue","gray"))
	}
return(plot_data)
}




EXP2_scores_full$species[EXP2_scores_full$species=="algae"]<-"_algae"
EXP2_scores_full$species[EXP2_scores_full$species=="blank"]<-"_blank"
EXP2_scores_full$species[EXP2_scores_full$species=="rotifers"]<-"_rotifers"
EXP2_scores_full$species[EXP2_scores_full$species=="ludwik"]<-"_blank"
EXP2_scores_full$species[EXP2_scores_full$species=="zywiec"]<-"_blank"


EXP2_scores_full$type[EXP2_scores_full$type=="zywiec"]<-"blank"
EXP2_scores_full$type[EXP2_scores_full$type=="ludwik"]<-"blank"
EXP2_scores_full$type[EXP2_scores_full$type=="ludwik"]<-"blank"


#EXP2_scores_decontaminated$species[EXP2_scores_decontaminated$species=="algae"]<-"_algae"
#EXP2_scores_decontaminated$species[EXP2_scores_decontaminated$species=="blank"]<-"_blank"
#EXP2_scores_decontaminated$species[EXP2_scores_decontaminated$species=="rotifers"]<-"_rotifers"

EXP2_scores_full$extraction[EXP2_scores_full$extraction==""]<-"3_"
EXP2_scores_full$extraction[EXP2_scores_full$extraction=="Chelex"]<-"1.Chelex"
EXP2_scores_full$extraction[EXP2_scores_full$extraction=="beads"]<-"2.beads"


EXP2_scores_decontaminated$extraction[EXP2_scores_decontaminated$extraction==""]<-"3_"
EXP2_scores_decontaminated$extraction[EXP2_scores_decontaminated$extraction=="Chelex"]<-"1.Chelex"
EXP2_scores_decontaminated$extraction[EXP2_scores_decontaminated$extraction=="beads"]<-"2.beads"



exp2_full<-pcoa_plot(EXP2_scores_full,"Experiment 2: raw data",c(20,19,15,4,1,17),c("chartreuse4","777777","darkorange3",col_vector),FALSE)
exp2_decontaminated<-pcoa_plot(EXP2_scores_decontaminated,"Experiment 2: decontaminated",c(20,19, 1),col_vector,FALSE)




ggarrange(exp2_full, exp2_decontaminated, ncol=2,common.legend = TRUE, legend="bottom")



dev.new()
svg("PCoA_exp2.svg",width=20,height=9)
ggarrange(exp2_full, exp2_decontaminated, ncol=2,common.legend = TRUE, legend="bottom")
dev.off()



