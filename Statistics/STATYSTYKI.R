
samples<- read.csv("../EXPERIMENT1/SAMPLES_EXPERIMENT1.csv",sep="\t", header=T)
percentages<- read.table("../EXPERIMENT1/Results/Statistics_table.txt",sep="\t", header=T,row.names=1)



## 4 - PCR CONTAMINATION 
## 5 extraction contamination
nr<-5

extraction_contaminatinon<- t(percentages[12,])
pcr_contaminatinon<- t(percentages[13,])



X<-data.frame( t(percentages[12,]), t(percentages[13,]) ,t(percentages[21,]),t(percentages[20,]) ,t(percentages[16,]),   t( samples[c(1:9), paste0(rownames(pcr_contaminatinon))])  )

colnames(X)<- c(
"PCR_contaminants",
"EXT_contaminants",
"Rotifers",
"Algae",
"Spikein",
"replicate",
"extraction",
"species",
"blank_type",
"pool_size",
"treatment",
"lab",
"x",
"experiment"
)


X$replicate<-as.factor(paste(X$species,X$replicate))

X$PCR_contaminants<-X$PCR_contaminants/100
X$EXT_contaminants<-X$EXT_contaminants/100
X$Spikein<-X$Spikein/100
X$Algae<-X$Algae/100
X$Rotifers<-X$Rotifers/100

X$ALL_contaminants<-X$PCR_contaminants+X$EXT_contaminants


exp1 <-subset(X,experiment=="exp1")
exp1 <-subset(exp1,blank_type=="")


exp1_detailed<-subset(exp1,species %in% c("Meb.KE.008","Mac.PL.015","Mac.PL.010"))

exp1_minimalistic<-subset(exp1,!species %in% c("Meb.KE.008","Mac.PL.015","Mac.PL.010"))


library(glmmTMB)
#####################
#EXP1



exp1_detailed$pool_size<-as.factor(exp1_detailed$pool_size)
exp1_minimalistic$pool_size<-as.factor(exp1_minimalistic$pool_size)

M<-subset(exp1_minimalistic, pool_size=="medium")
I<-subset(exp1_minimalistic, pool_size=="1")
					
boxplot(ALL_contaminants~pool_size, data=exp1_minimalistic)

wilcox.test(M$ALL_contaminants, I$ALL_contaminants, paired = TRUE, alternative = "two.sided")

glmmTMB
+(1|replicate)
species+pool_size:species
model_ext_exp1_detailed<-glmmTMB(ALL_contaminants~extraction+pool_size+species+(1|replicate),data=exp1_detailed, family=list(family = "beta", link = "logit") )

hist(t(simulate(model_ext_exp1_detailed,nsim=1000)))


library(glmmTMB)
#library(DHARMa)
#	simulateResiduals(model_ext_exp1_detailed, plot = T)


#
ggplot(exp1_minimalistic, aes(y = PCR_contaminants, col = pool_size))+geom_bar(stat="identity")

model_pcr_exp1_detailed<-glmmTMB(PCR_contaminants~extraction+pool_size+species+(1|replicate) ,data=exp1_detailed)
summary(model_pcr_exp1_detailed)

model_all_exp1_detailed<-glmmTMB(ALL_contaminants~extraction+pool_size+species+(1|replicate) ,data=exp1_detailed)
summary(model_all_exp1_detailed)


#####################
#EXP3 (2)

exp3$treatment<-as.character(exp3$treatment)


exp3$sample_type<-""
exp3[exp3$pool_size=="medium",]$sample_type<-"medium"
exp3[exp3$pool_size!="medium",]$sample_type<-"individuals"


)+(1|pool_size/sample_type)
(1+sample_type|treatment)++treatment
+(1+sample_type|treatment)

model_z_medium<-glmmTMB(ALL_contaminants~extraction+species+pool_size+treatment+(1|pool_size/treatment)+(1|replicate),data=exp3,family=list(family = "beta", link = "logit"))
summary(model_z_medium)

+(1|pool_size/treatment)
SPIK<-glmmTMB(Spikein~extraction+species+pool_size+treatment+(1|replicate),data=exp3,family=list(family = "beta", link = "logit"))
summary(SPIK)


model_z_medium$drop.unused.levels <- TRUE

model<-glmmTMB(contaminants~extraction+pool_size+species+treatment+(1|replicate) ,data=subset(exp3,pool_size!="medium"))

summary(model)






