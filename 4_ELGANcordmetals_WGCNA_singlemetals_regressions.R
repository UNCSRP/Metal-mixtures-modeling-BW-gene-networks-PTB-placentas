
#################################################################################################
#################################################################################################
#### ELGAN cord metals gene expression: evaluating placental transcriptomic responses to cord 
#### metal levels 
#### Part 3: WGCNA analysis. A) Single metals regression analysis  
####
#### Code drafted by Lauren Eaves
#### Lasted updated: October 13th 2021
#################################################################################################
#################################################################################################

# Clean your working environment
rm(list=ls())

#################################################################################################
#################################################################################################
#### Installing and activating packages, and setting your working directory
#################################################################################################
#################################################################################################

# Activate the appropriate packages:
library(tidyverse)

# Set your working directory to the folderpath containing your input files:
setwd("#yourwd")

# Create an output folder
Output_Folder <- ("#youroutputfolder")

# Create a current date variable to name outputfiles
cur_date <- str_replace_all(Sys.Date(),"-","")


#################################################################################################
#################################################################################################
#### Loading, organizing, and initial visualizations of data
#################################################################################################
#################################################################################################

MEs <- read.csv(file="#input file generated from script 3_ELGANcordmetals_WGCNA_generatemodules, titled DATE_MEs_AggregateValues_bysubject.csv") #module values for each subject 
cohort <-read.csv(file="#input data file generated from script 2_ELGANcordmetals_multimetals_DESeq2 names DATE_SubjectInfo_SubjectsIncludedinModel") 

#merge the datasets 
MEs <- MEs %>% dplyr::rename(ID = X)
colnames(cohort)
cohort <- as.data.frame(cohort) %>% 
  dplyr::select(ID, sex, score, bmicat, smoke, Mn_ugg_log, Cu_ugg_log, Zn_ugg_log,As_ngg_log,Se_ugg_log,Sr_ugg_log,Cd_ngg_log,
        Sb_ngg_log,Ba_ngg_log,Hg_ngg_log,Pb_ngg_log,As_cat,Ba_cat,Cd_cat,Cu_cat,Hg_cat, Mn_cat,Pb_cat,Sb_cat,Se_cat,Sr_cat,Zn_cat, SV1, SV2, SV3)
full <- full_join(MEs, cohort, by="ID")

#################################################################################################
#################################################################################################
#### Running crude models with single metals 
#################################################################################################
#################################################################################################

#continuous metals
metalslog <- c("Mn_ugg_log", "Cu_ugg_log", "Zn_ugg_log","As_ngg_log","Se_ugg_log","Sr_ugg_log","Cd_ngg_log",
               "Sb_ngg_log","Ba_ngg_log","Hg_ngg_log","Pb_ngg_log")

modules <- colnames(MEs)[2:26]
print(modules)

results_cont <- data.frame()
for (i in 1:length(metalslog)) {
  metal <- metalslog[[i]]
  metal <- paste0(as.name(metal))
  print(metal)
  
  for (j in 1:length(modules)){
    module <- modules[[j]]
    print(module)
    
    output.lm.1 <- summary(lm(eval(parse(text = paste0(module,"~",metal))), data=full))$coefficients[2,1:4] %>%as.data.frame()
    colnames(output.lm.1)[1] <- paste0(metal,".",module)
    output.lm.1 <- output.lm.1 %>%  as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column(var="MetalModule")
    results_cont<-rbind(results_cont,output.lm.1)
  }
}
results_cont <- results_cont %>% 
  dplyr::rename("p.value" ="Pr(>|t|)") %>% 
  dplyr::rename("t.value" ="t value") %>% 
  mutate(sigchange.pvalue = ifelse(p.value <0.05,1,0)) %>% 
  mutate(adj.p = p.adjust(p.value, method="BH", n=275)) %>% 
  mutate(sigchange.adjpvalue = ifelse(adj.p <0.05,1,0))

#categorical metals
metalscat <- c("As_cat","Ba_cat","Cd_cat","Cu_cat","Hg_cat", "Mn_cat","Pb_cat","Sb_cat","Se_cat","Sr_cat","Zn_cat")

results_cat <- data.frame()

for (i in 1:length(metalscat)) {
  metal <- metalscat[[i]]
  metal <- paste0(as.name(metal))
  print(metal)
  
  for (j in 1:length(modules)){
    module <- modules[[j]]
    print(module)
    
    output.lm.1 <- summary(lm(eval(parse(text = paste0(module,"~",metal))), data=full))$coefficients[2,1:4] %>%as.data.frame()
    colnames(output.lm.1)[1] <- paste0(metal,".",module)
    output.lm.1 <- output.lm.1 %>%  as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column(var="MetalModule")
    results_cat<-rbind(results_cat,output.lm.1)
  }
}
results_cat <- results_cat %>% 
  dplyr::rename("p.value" ="Pr(>|t|)") %>% 
  dplyr::rename("t.value" ="t value") %>% 
  mutate(sigchange.pvalue = ifelse(p.value <0.05,1,0)) %>% 
  mutate(adj.p = p.adjust(p.value, method="BH", n=275)) %>% 
  mutate(sigchange.adjpvalue = ifelse(adj.p <0.05,1,0))

#################################################################################################
#################################################################################################
#### Running adjusted models with single metals 
#################################################################################################
#################################################################################################

#continous metals
results_cont_adj <- data.frame()
for (i in 1:length(metalslog)) {
  metal <- metalslog[[i]]
  metal <- paste0(as.name(metal))
  print(metal)
  
  for (j in 1:length(modules)){
    module <- modules[[j]]
    print(module)
    
    output.lm.1 <- summary(lm(eval(parse(text = paste0(module,"~",metal,"+ score + bmicat + smoke + sex + SV1 +SV2"))), data=full))$coefficients[2,1:4] %>%as.data.frame()
    colnames(output.lm.1)[1] <- paste0(metal,".",module)
    output.lm.1 <- output.lm.1 %>%  as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column(var="MetalModule")
    results_cont_adj<-rbind(results_cont_adj,output.lm.1)
  }
}
results_cont_adj <- results_cont_adj %>% 
  dplyr::rename("p.value" ="Pr(>|t|)") %>% 
  dplyr::rename("t.value" ="t value") %>% 
  mutate(sigchange.pvalue = ifelse(p.value <0.05,1,0)) %>% 
  mutate(adj.p = p.adjust(p.value, method="BH", n=275)) %>% 
  mutate(sigchange.adjpvalue = ifelse(adj.p <0.05,1,0))


#categorical metals
results_cat_adj <- data.frame()
for (i in 1:length(metalscat)) {
  metal <- metalscat[[i]]
  metal <- paste0(as.name(metal))
  print(metal)
  
  for (j in 1:length(modules)){
    module <- modules[[j]]
    print(module)
    
    output.lm.1 <- summary(lm(eval(parse(text = paste0(module,"~",metal,"+ score + bmicat + smoke + sex+ SV1 +SV2"))), data=full))$coefficients[2,1:4] %>%as.data.frame()
    colnames(output.lm.1)[1] <- paste0(metal,".",module)
    output.lm.1 <- output.lm.1 %>%  as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column(var="MetalModule")
    results_cat_adj<-rbind(results_cat_adj,output.lm.1)
  }
}
results_cat_adj <- results_cat_adj %>% 
  dplyr::rename("p.value" ="Pr(>|t|)") %>% 
  dplyr::rename("t.value" ="t value") %>% 
  mutate(sigchange.pvalue = ifelse(p.value <0.05,1,0)) %>% 
  mutate(adj.p = p.adjust(p.value, method="BH", n=275)) %>% 
  mutate(sigchange.adjpvalue = ifelse(adj.p <0.05,1,0))



#################################################################################################
#################################################################################################
#### Running crude models with single metals, using quantized exposures  
#################################################################################################
#################################################################################################

library(qgcomp)

#use qgcomp to generate quantized exposures 
#vector of exposures
Xnm <- c("Mn_ugg_log", "Cu_ugg_log", "Zn_ugg_log","As_ngg_log","Se_ugg_log","Sr_ugg_log","Cd_ngg_log",
         "Sb_ngg_log","Ba_ngg_log","Hg_ngg_log","Pb_ngg_log")
#vector of covariates
covars = c("sex", "score", "bmicat", "smoke")

qc.fit <- qgcomp.noboot(ME1 ~ ., expnms=Xnm,
                        dat=full[,c(Xnm,'ME1')],family=gaussian(),q=4)
  
head(qc.fit$qx)
quantizedexposures <- as.data.frame(qc.fit$qx)
full <- cbind(full,quantizedexposures)

metalsq <- colnames(quantizedexposures)

#unadjusted quantized exposures 
results_q <- data.frame()
for (i in 1:length(metalsq)) {
  metal <- metalsq[[i]]
  metal <- paste0(as.name(metal))
  print(metal)
  
  for (j in 1:length(modules)){
    module <- modules[[j]]
    print(module)
    
    output.lm.1 <- summary(lm(eval(parse(text = paste0(module,"~",metal))), data=full))$coefficients[2,1:4] %>%as.data.frame()
    colnames(output.lm.1)[1] <- paste0(metal,".",module)
    output.lm.1 <- output.lm.1 %>%  as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column(var="MetalModule")
    results_q<-rbind(results_q,output.lm.1)
  }
}
results_q <- results_q %>% 
  dplyr::rename("p.value" ="Pr(>|t|)") %>% 
  dplyr::rename("t.value" ="t value") %>% 
  mutate(sigchange.pvalue = ifelse(p.value <0.05,1,0)) %>% 
  mutate(adj.p = p.adjust(p.value, method="BH", n=275)) %>% 
  mutate(sigchange.adjpvalue = ifelse(adj.p <0.05,1,0))


#adjusted
results_q_adj <- data.frame()
for (i in 1:length(metalsq)) {
  metal <- metalsq[[i]]
  metal <- paste0(as.name(metal))
  print(metal)
  
  for (j in 1:length(modules)){
    module <- modules[[j]]
    print(module)
    
    output.lm.1 <- summary(lm(eval(parse(text = paste0(module,"~",metal,"+ score + bmicat + smoke + sex+ SV1 +SV2"))), data=full))$coefficients[2,1:4] %>%as.data.frame()
    colnames(output.lm.1)[1] <- paste0(metal,".",module)
    output.lm.1 <- output.lm.1 %>%  as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column(var="MetalModule")
    results_q_adj<-rbind(results_q_adj,output.lm.1)
  }
}
results_q_adj <- results_q_adj %>% 
  dplyr::rename("p.value" ="Pr(>|t|)") %>% 
  dplyr::rename("t.value" ="t value") %>% 
  mutate(sigchange.pvalue = ifelse(p.value <0.05,1,0)) %>% 
  mutate(adj.p = p.adjust(p.value, method="BH", n=2775)) %>% 
  mutate(sigchange.adjpvalue = ifelse(adj.p <0.05,1,0))

#################################################################################################
#################################################################################################
#### Export results as .csv's
#################################################################################################
#################################################################################################

write.csv(results_cont,paste0(Output_Folder,"/", cur_date,"_linearregression_singlemetals_continous_crude.csv"))
write.csv(results_cont_adj,paste0(Output_Folder,"/", cur_date,"_linearregression_singlemetals_continous_adjusted.csv"))
write.csv(results_cat,paste0(Output_Folder,"/", cur_date,"_linearregression_singlemetals_categorical_crude.csv"))
write.csv(results_cat_adj,paste0(Output_Folder,"/", cur_date,"_linearregression_singlemetals_categorical_adjusted.csv"))
write.csv(results_q,paste0(Output_Folder,"/", cur_date,"_linearregression_singlemetals_quantized_crude.csv"))
write.csv(results_q_adj,paste0(Output_Folder,"/", cur_date,"_linearregression_singlemetals_quantized_adjusted.csv"))

