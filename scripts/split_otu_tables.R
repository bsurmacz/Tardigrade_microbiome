
library(dplyr)

zotus<-read.csv2("zotu_table.txt",sep="\t",header=TRUE)
info<-zotus[1]
colnames(info)="#OTU ID"
zotus<-zotus[,-1]

EXP1<- zotus %>% select(,matches("^EXPERIMENT_1"))
EXP2<- zotus %>% select(,matches("^EXPERIMENT_2"))
EXP3<- zotus %>% select(,matches("^EXPERIMENT_3"))
EXP4<- zotus %>% select(,matches("^EXPERIMENT_4"))

EXP1<-cbind(info,EXP1)
EXP2<-cbind(info,EXP2)
EXP3<-cbind(info,EXP3)
EXP4<-cbind(info,EXP4)

write.table(EXP1,file="zotu_table_experiment_1.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(EXP2,file="zotu_table_experiment_2.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(EXP3,file="zotu_table_experiment_3.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(EXP4,file="zotu_table_experiment_4.txt",sep="\t",col.names=T,row.names=F,quote=F)


otus<-read.csv2("otu_table.txt",sep="\t",header=TRUE)
info<-otus[1]
colnames(info)="#OTU ID"
otus<-otus[,-1]

EXP1<- otus %>% select(,matches("^EXPERIMENT_1"))
EXP2<- otus %>% select(,matches("^EXPERIMENT_2"))
EXP3<- otus %>% select(,matches("^EXPERIMENT_3"))
EXP4<- otus %>% select(,matches("^EXPERIMENT_4"))

EXP1<-cbind(info,EXP1)
EXP2<-cbind(info,EXP2)
EXP3<-cbind(info,EXP3)
EXP4<-cbind(info,EXP4)

write.table(EXP1,file="otu_table_experiment_1.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(EXP2,file="otu_table_experiment_2.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(EXP3,file="otu_table_experiment_3.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(EXP4,file="otu_table_experiment_4.txt",sep="\t",col.names=T,row.names=F,quote=F)




