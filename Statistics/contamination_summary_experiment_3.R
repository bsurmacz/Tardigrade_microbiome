library(reshape2)
library(ggplot2)
library(dplyr)

samples<- read.csv("../EXPERIMENT3/SAMPLES_EXPERIMENT3.csv",sep="\t", header=T)
Statistic_table<- read.csv("../EXPERIMENT3/Results/Statistics_table.txt",sep="\t",header=T,row.names=1,dec=".")


percentages<-Statistic_table[c(12,22,21,13,14,17),  ] 
percentages


sample_type<-t(samples[5, paste0(colnames(Statistic_table))])

percentages<-as.data.frame(t(percentages))

percentages$sample_type<-sample_type



colnames(percentages)[colnames(percentages)=="%_of_all_Symbionts_reads_in_lib"]<-"Symbionts"
colnames(percentages)[colnames(percentages)=="%_of_all_Rotifers_reads_in_lib"]<-"Rotifers"
colnames(percentages)[colnames(percentages)=="%_of_Extraction_contamination_reads_in_lib"]<-"EXT_contam"
colnames(percentages)[colnames(percentages)=="%_of_PCR_contamination_reads_in_lib"]<-"PCR_contam"
colnames(percentages)[colnames(percentages)=="%_of_Extraction_Spikein_reads_in_lib"]<-"EXT_spikein"
colnames(percentages)[colnames(percentages)=="%_of_all_Algae_reads_in_lib"]<-"Algae"


OUT<-percentages %>% group_by(sample_type) %>% summarise (N=n(), PCR_median=median(PCR_contam), PCR_min=min(PCR_contam),PCR_max=max(PCR_contam),EXT_median=median(EXT_contam), EXT_min=min(EXT_contam),EXT_max=max(EXT_contam),ROT_median=median(Rotifers),ROT_min=min(Rotifers),ROT_max=max(Rotifers),ALG_median=median(Algae),ALG_min=min(Algae),ALG_max=max(Algae))
OUT



write.table(OUT,file="summary_EXPERIMENT3.csv")
