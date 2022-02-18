library(ggpubr)

library("vegan")
library("ggplot2")


samples<- read.csv("../../EXPERIMENT1/SAMPLES_EXPERIMENT1.csv",sep="\t", header=T)
full_data<-read.table("../../EXPERIMENT1/Results/zotu_table_expanded.txt",sep="\t",header=T)
decontaminated_data<-read.table("../../EXPERIMENT1/Results/Decontaminated_zOTU_table.txt",sep="\t",header=T)

excluded<-c("EXPERIMENT_1_T_AU_031_1","EXPERIMENT_1_T_AU_031_2")

full_data<-full_data[, ! colnames(full_data) %in% excluded ]
decontaminated_data<-decontaminated_data[, ! colnames(decontaminated_data) %in% excluded ]


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

EXP1_decontaminated<-relative_decontaminated[, EXPERIMENT_decontaminated=="exp1" & !(SAMPLE_TYPE  %in% c( "blank","rotifers","algae","ludwik","zywiec"))]
detaied_EXP1_decontaminated<-relative_decontaminated[, EXPERIMENT_decontaminated=="exp1" & SPECIES_decontaminated %in% c("Meb.KE.008","Mac.PL.010","Mac.PL.015") ]
cursory_EXP1_decontaminated<-relative_decontaminated[, EXPERIMENT_decontaminated=="exp1" & !(SPECIES_decontaminated %in% c("Meb.KE.008","Mac.PL.010","Mac.PL.015","algae","rotifers","blank")) ]


EXP1_full<-relative_full[, EXPERIMENT_full=="exp1" & !(SAMPLE_TYPE_full  %in% c( "blank","rotifers","algae","ludwik","zywiec"))]
detaied_EXP1_full<-relative_full[, EXPERIMENT_full=="exp1" & SPECIES_full %in% c("Meb.KE.008","Mac.PL.010","Mac.PL.015","blank","algae","rotifers") ]
cursory_EXP1_full<-relative_full[, EXPERIMENT_full=="exp1" & !(SPECIES_full %in% c("Meb.KE.008","Mac.PL.010","Mac.PL.015")) ]


EXP1_decontaminated<-t(EXP1_decontaminated)
detaied_EXP1_decontaminated<-t(detaied_EXP1_decontaminated)
detaied_EXP1_full<-t(detaied_EXP1_full)
cursory_EXP1_decontaminated<-t(cursory_EXP1_decontaminated)
cursory_EXP1_full<-t(cursory_EXP1_full)


EXP1_full<-t(EXP1_full)


max<-1000
EXP1_MDS_decontaminated<-metaMDS( EXP1_decontaminated, distance="bray", k=2, trymax=max)

detailed_EXP1_MDS_decontaminated<-metaMDS(detaied_EXP1_decontaminated, distance="bray", k=2, trymax=max)
detailed_EXP1_MDS_full<-metaMDS(detaied_EXP1_full, distance="bray", k=2, trymax=max)

cursory_EXP1_MDS_decontaminated<-metaMDS(cursory_EXP1_decontaminated, distance="bray", k=2, trymax=max)
cursory_EXP1_MDS_full<-metaMDS(cursory_EXP1_full, distance="bray", k=2, trymax=max)


EXP1_MDS_full<-metaMDS( EXP1_full, distance="bray", k=2, trymax=max)

EXP1_scores_decontaminated <- as.data.frame(scores(EXP1_MDS_decontaminated)) 

detailed_EXP1_scores_decontaminated <- as.data.frame(scores(detailed_EXP1_MDS_decontaminated)) 
detailed_EXP1_scores_full <- as.data.frame(scores(detailed_EXP1_MDS_full)) 
cursory_EXP1_scores_decontaminated <- as.data.frame(scores(cursory_EXP1_MDS_decontaminated)) 
cursory_EXP1_scores_full <- as.data.frame(scores(cursory_EXP1_MDS_full)) 


EXP1_scores_full <- as.data.frame(scores(EXP1_MDS_full)) 


#######333   KOLORY
n <- 52
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
set.seed(117)
col_vector<-sample(color, n)
col_vector<-c("orange","green","#ca8b44",col_vector)

#########

EXP1_scores_decontaminated$species<-t( samples[3, paste0(rownames(EXP1_scores_decontaminated))])

detailed_EXP1_scores_decontaminated$species<-t( samples[3, paste0(rownames(detailed_EXP1_scores_decontaminated))])
detailed_EXP1_scores_full$species<-t( samples[3, paste0(rownames(detailed_EXP1_scores_full))])

detailed_EXP1_scores_decontaminated$extraction<-t( samples[2, paste0(rownames(detailed_EXP1_scores_decontaminated))])
detailed_EXP1_scores_full$extraction<-t( samples[2, paste0(rownames(detailed_EXP1_scores_full))])



cursory_EXP1_scores_decontaminated$species<-t( samples[3, paste0(rownames(cursory_EXP1_scores_decontaminated))])
cursory_EXP1_scores_full$species<-t( samples[3, paste0(rownames(cursory_EXP1_scores_full))])


EXP1_scores_decontaminated$type<-t( samples[5, paste0(rownames(EXP1_scores_decontaminated))])
detailed_EXP1_scores_decontaminated$type<-t( samples[5, paste0(rownames(detailed_EXP1_scores_decontaminated))])
detailed_EXP1_scores_full$type<-t( samples[5, paste0(rownames(detailed_EXP1_scores_full))])

cursory_EXP1_scores_decontaminated$type<-t( samples[5, paste0(rownames(cursory_EXP1_scores_decontaminated))])
cursory_EXP1_scores_full$type<-t( samples[5, paste0(rownames(cursory_EXP1_scores_full))])

EXP1_scores_full$species<-t( samples[3, paste0(rownames(EXP1_scores_full))])
EXP1_scores_full$type<-t( samples[5, paste0(rownames(EXP1_scores_full))])

library(dplyr)
 library(ggnewscale)

#by species
pcoa_plot<-function(data,Title,scale_shape,colors,extraction){
	hulls<-data %>%group_by(species) %>% slice(chull(NMDS1,NMDS2))

	plot_data<-ggplot( aes(x=NMDS1,y=NMDS2,colour=species,shape=type),data=data)+ geom_point(data=data,size=5)+ geom_polygon(data =hulls, aes(group=species,fill=species,color=species), alpha = 0.1)+scale_fill_manual(values=colors)+scale_color_manual(values=colors)+theme_bw()+scale_shape_manual(values=scale_shape)+labs(title=Title)
	if(extraction){
#	plot_data<-plot_data+new_scale_color() +geom_point(shape = 1, aes(x=NMDS1,y=NMDS2,color=extraction),size=7)+scale_color_manual(values=c("red","blue","gray"))
	}
return(plot_data)
}



copy_cursory_EXP1_scores_full<-cursory_EXP1_scores_full



cursory_EXP1_scores_full$species[cursory_EXP1_scores_full$species=="algae"]<-"_algae"
cursory_EXP1_scores_full$species[cursory_EXP1_scores_full$species=="blank"]<-"_blank"
cursory_EXP1_scores_full$species[cursory_EXP1_scores_full$species=="rotifers"]<-"_rotifers"

detailed_EXP1_scores_full$species[detailed_EXP1_scores_full$species=="algae"]<-"_algae"
detailed_EXP1_scores_full$species[detailed_EXP1_scores_full$species=="blank"]<-"_blank"
detailed_EXP1_scores_full$species[detailed_EXP1_scores_full$species=="rotifers"]<-"_rotifers"



detailed_EXP1_scores_decontaminated$extraction[detailed_EXP1_scores_decontaminated$extraction==""]<-"3_"
detailed_EXP1_scores_decontaminated$extraction[detailed_EXP1_scores_decontaminated$extraction=="Chelex"]<-"1.Chelex"
detailed_EXP1_scores_decontaminated$extraction[detailed_EXP1_scores_decontaminated$extraction=="beads"]<-"2.beads"


detailed_EXP1_scores_full$extraction[detailed_EXP1_scores_full$extraction==""]<-"3_"
detailed_EXP1_scores_full$extraction[detailed_EXP1_scores_full$extraction=="Chelex"]<-"1.Chelex"
detailed_EXP1_scores_full$extraction[detailed_EXP1_scores_full$extraction=="beads"]<-"2.beads"


col_vector[16]<-"yellow4"
col_vector[17]<-"khaki"
col_vector[23]<-"tomato2"
col_vector[25]<-"maroon"
col_vector[27]<-"rosybrown1"


colorscale<-col_vector
names(colorscale)<-unique(cursory_EXP1_scores_full$species)

colorscale["_algae"]<-"#458a00"
colorscale["_rotifers"]<-"darkorange3"
colorscale["_blank"]<-"#777777"
colorscale["Macrobiotus cf. recens PT.056"]<-"#ad3b9c"
colorscale["Mesobiotus sp. VN.036"]<-"#458a00"
colorscale["Macrobiotus canaricus ES.004"]<-"#9c8936"
colorscale["Milnesium tardigradum FR.072"]<-"#5c3a1a"
colorscale["Paramacrobiotus sp.  FR.123"]<-"#dcc5e6"
colorscale["Milnesium inceptum NZ.030"]<-"#780128"
colorscale["Paramacrobiotus sp. PT.035"]<-"linen"
colorscale["Mesobiotus sp. RU.017"]<-"#7adecd"
colorscale["Paramacrobiotus sp. PT.006"]<-"#c6f7ef"
colorscale["Paramacrobiotus fairbanksi PL.018"]<-"#e3bfcd"
colorscale["Hypsibius sp.  PL.016"]<-"#c4620c"
colorscale["Macrobiotus polonicus SK.003"]<-"#fce4a2"
colorscale["Pseudobiotus sp.  PL.318"]<-"#7b856d"
colorscale["Macrobiotus polonicus PT.008"]<-"#d45500"


exp1_a_full<-pcoa_plot(cursory_EXP1_scores_full,"Experiment 1B: raw data",c(16, 18,4,1,17),colorscale,FALSE)
exp1_a_decontaminated<-pcoa_plot(cursory_EXP1_scores_decontaminated,"Experiment 1b: decontaminated",c(16, 1),colorscale,FALSE)
exp1_a_full



ggarrange(exp1_a_full, exp1_a_decontaminated, ncol=2,common.legend = TRUE, legend="bottom")

##

detailed_EXP1_scores_full$type[  detailed_EXP1_scores_full$type %in% c("1","10","medium")]<-paste( detailed_EXP1_scores_full$type[  detailed_EXP1_scores_full$type %in% c("1","10","medium")] , detailed_EXP1_scores_full$extraction[  detailed_EXP1_scores_full$type %in% c("1","10","medium")] )

detailed_EXP1_scores_decontaminated$type[  detailed_EXP1_scores_decontaminated$type %in% c("1","10","medium")]<-paste( detailed_EXP1_scores_decontaminated$type[  detailed_EXP1_scores_decontaminated$type %in% c("1","10","medium")] , detailed_EXP1_scores_decontaminated$extraction[  detailed_EXP1_scores_decontaminated$type %in% c("1","10","medium")] )



shapescale<-c(4, 18,15,0,17,16,1)
names(shapescale)<-unique(detailed_EXP1_scores_full$type)

colorscale<-c("#777777","#458a00", "#9fb87f","#f5b433","#c48eed","darkorange3" )
names(colorscale)<-unique(detailed_EXP1_scores_full$species)



exp1_b_full<-pcoa_plot(detailed_EXP1_scores_full,"Experiment 1a: raw data",shapescale,colorscale,F)
exp1_b_decontaminated<-pcoa_plot(detailed_EXP1_scores_decontaminated,"Experiment 1a: decontaminated",shapescale,colorscale,F)
ggarrange(exp1_b_full, exp1_b_decontaminated, ncol=2,common.legend = TRUE, legend="bottom")


detailed_EXP1_scores_full$extraction


exp1_a_full

dev.new()
svg("PCoA_exp1a_new.svg",width=20,height=9)
ggarrange(exp1_a_full, exp1_a_decontaminated, ncol=2,common.legend = TRUE, legend="bottom")
dev.off()

dev.new()
svg("PCoA_exp1b_new.svg",width=20,height=9)
ggarrange(exp1_b_full, exp1_b_decontaminated, ncol=2,common.legend = TRUE, legend="bottom")
dev.off()



