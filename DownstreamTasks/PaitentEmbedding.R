library(data.table)
library(stringr)
library(dplyr)
library(ggsurvfit)
library(ggplot2)

#Read the patient visit level occurrence data 
UMPC = fread('ADPatients.csv')
UMPC$patient_num = as.character(UMPC$patient_num)

#Load the knowledge graph obtained by PMI testing and LASSO regression 
load("cosine.Rdata")
load("CodeNetwork.sym.weight.Rdata")
ADcode = CodeNetwork.sym.weight[which(rownames(CodeNetwork.sym.weight) == "PheCode:290.11"),]
ADcode = ADcode[ADcode!=0]
ADcode = c(1,ADcode)
names(ADcode)[1] = 'PheCode:290.11'

#Load the embeddings and frequency 
load("embedding.Rdata")
load("dictfull.Rdata")

d = 1500 
PE = matrix(0, nrow = length(Patients), ncol = d)
rownames(PE) = Patients

N = sum(dict$rowsum)
embed = embed[rownames(embed)%in%names(ADcode),]
dict = dict[match(rownames(embed),dict$code),]

UMPCsub = UMPC %>% filter(feature_id %in% intersect(names(ADcode),dict$code))
rm(CodeNetwork.sym.weight,omop,omop.sub)

UMPCcounts = UMPCsub %>% group_by(feature_id) %>% summarize(counts = sum(counts))


save(N, d, dict, embed, Patients, PE, UMPCsub, UMPC, ADcode, UMPCcounts,
     file = "PreData.Rdata")

# after pre-process data, start from here ######
load("PreData.Rdata")
# Generate embeddings for Patients
Pstat = NULL
for(i in 1:length(Patients)){
  cat(i,"; ")
  if(i%%10==0) cat("\n")
  id = Patients[i]
  patientCode <- UMPCsub[UMPCsub$patient_num == id,]
  if(nrow(patientCode) == 0){
    Pstat = rbind(Pstat,c(id,0,0,0))
  }else{
    tmpid = match(patientCode$feature_id,rownames(embed))
    w = log(patientCode$counts+1) * log(N/dict$rowsum[tmpid])
    if(length(w)>1){
      PE[i,] <- colSums(w*embed[tmpid,])
    }else{
      PE[i,] <- w*embed[tmpid,] 
    }
    Pstat = rbind(Pstat,c(id,length(w),sum(patientCode$counts),sum(w)))
  }
}

save(PE,Pstat,
     file='PE_UPMC_AD.Rdata')

###Demographic
load("PE_UPMC_AD.Rdata") 
load("demographic.Rdata")
file = "patient_dates_info.csv"
demographic1 = fread(file)

demographic1$date_first_enc = as.Date(demographic1$date_first_enc,format="%Y-%m-%d") 
demographic1$date_last_followup = as.Date(demographic1$date_last_followup,format="%Y-%m-%d")
demographic1$date_first_phecode = as.Date(demographic1$date_first_phecode,format="%Y-%m-%d")
demographic1$date_last_followup_phecode = as.Date(demographic1$date_last_followup_phecode,format="%Y-%m-%d")
demographic1$death_date = as.Date(demographic1$death_date,format="%Y-%m-%d")
demographic1$followup = demographic1$date_last_followup - demographic1$date_first_enc
demographic1$survival = demographic1$death_date - demographic1$date_first_enc 

demographic1$event = 1
demographic1$event[is.na(demographic1$survival)] = 0
demographic1$event[which(demographic1$followup < demographic1$survival)] = 0
demographic1 = demographic1[demographic1$followup>30,]

demographic = demographic1[!is.na(demographic1$date_first_phecode),]
demographic = demographic[!is.na(demographic$date_last_followup_phecode),]

demographic$followup[which(demographic$followup > demographic$survival)] = demographic$survival[which(demographic$followup > demographic$survival)]
demographic = demographic[ is.na(demographic$survival)| (demographic$survival>=0),]

demographic$followup = as.integer(demographic$followup)
demographic$patient_num = as.character(demographic$patient_num)

PE.mor = PE[match(demographic$patient_num,rownames(PE)),]
PE.mor = PE.mor/demographic$followup


##Run k-means for different k
set.seed(1)
kmeans.fit = kmeans(PE.mor, centers=2, iter.max = 20, nstart = 10)
demographic$cluster = kmeans.fit$cluster

save(demographic,file='demographic.Rdata')
save(UMPCsub,file='UPMCsub.Rdata')
save(PE.mor,file='PE.mor.Rdata')

load(file='UPMC_Final.Rdata')

size = NULL
follow = NULL
survival = NULL
mortality = NULL 
for(j in 1:2){
  size = c(size,sum(demographic$cluster==j))
  follow = c(follow,mean(demographic$followup[demographic$cluster==j]))
  survival = c(survival,mean(demographic$survival[demographic$cluster==j],na.rm=T))
  mortality = c(mortality,mean(demographic$event[demographic$cluster==j]))
}
size = c(size,length(demographic$cluster))
follow = c(follow,mean(demographic$followup))
survival = c(survival,mean(demographic$survival,na.rm=T))
mortality = c(mortality,mean(demographic$event))

Res = rbind(size,follow,survival,mortality)
colnames(Res) = c(paste0('Group ', 1:2),'overall')

load("dictfull.Rdata")
load("embedding.Rdata")
load("PreData.Rdata")

ADfeatures = sort(unique(UMPCsub$feature_id))
distAD = dict[match(ADfeatures,dict$code),]
embed = embed[match(distAD$code,rownames(embed)),]
