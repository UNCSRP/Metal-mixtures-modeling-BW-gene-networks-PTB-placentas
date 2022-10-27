
#################################################################################################
#################################################################################################
#### ELGAN cord metals gene expression: evaluating placental transcriptomic responses to cord 
#### metal levels 
#### Part 3: WGCNA analysis. B) Q-gcomp multi-metal analysis 
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
library("qgcomp")
library(ggpubr)

# Set your working directory to the folderpath containing your input files:
setwd("#yourwd")

# Create an output folder
Output_Folder <- ("#youroutputfolder")

# Create a current date variable to name outputfiles
cur_date <- paste0(str_replace_all(Sys.Date(),"-",""),"allmetals")


#################################################################################################
#################################################################################################
#### Loading, organizing, and initial visualizations of data
#################################################################################################
#################################################################################################

MEs <- read.csv(file="#input file generated from script 3_ELGANcordmetals_WGCNA_generatemodules, titled DATE_MEs_AggregateValues_bysubject.csv") #module values for each subject 
cohort <-read.csv(file="#input data file generated from script 2_ELGANcordmetals_multimetals_DESeq2 names DATE_SubjectInfo_SubjectsIncludedinModel") 

#merge the datasets 
MEs <- MEs %>% dplyr::rename(ID = X)
cohort <- cohort %>% 
  dplyr::select(ID, sex, score, bmicat, smoke, Mn_ugg_log, Cu_ugg_log, Zn_ugg_log,As_ngg_log,Se_ugg_log,Sr_ugg_log,Cd_ngg_log,
         Sb_ngg_log,Ba_ngg_log,Hg_ngg_log,Pb_ngg_log,As_cat,Ba_cat,Cd_cat,Cu_cat,Hg_cat, Mn_cat,Pb_cat,Sb_cat,Se_cat,Sr_cat,Zn_cat, 
         PC1, PC2, SV1, SV2, SV3, toxicindex,essentialindex)
full <- full_join(MEs, cohort, by="ID")

#################################################################################################
#################################################################################################
#### Mixtures modelling utilizing quantile-based g-computation 
#################################################################################################
#################################################################################################

#################################################################################################
#### Websites describing qgcomp:
#### https://cran.r-project.org/web/packages/qgcomp/vignettes/qgcomp-vignette.html
#### https://cran.r-project.org/web/packages/qgcomp/index.html
#### https://cran.r-project.org/web/packages/qgcomp/qgcomp.pdf
#### https://arxiv.org/pdf/1902.04200.pdf
#### https://www.rdocumentation.org/packages/qgcomp/versions/1.0.0
#################################################################################################

#vector of exposures
Xnm <- c("Mn_ugg_log", "Cu_ugg_log", "Zn_ugg_log","As_ngg_log","Se_ugg_log","Sr_ugg_log","Cd_ngg_log",
         "Sb_ngg_log","Ba_ngg_log","Hg_ngg_log","Pb_ngg_log")
toxicXnm <- c("Mn_ugg_log","As_ngg_log","Sr_ugg_log","Cd_ngg_log",
              "Sb_ngg_log","Ba_ngg_log","Hg_ngg_log","Pb_ngg_log")
essentialXnm <- c("Cu_ugg_log", "Zn_ugg_log","Se_ugg_log")

#vector of covariates
covars = c("sex", "score", "bmicat", "smoke", "SV1","SV2")
#vector of modules 
modules <- colnames(MEs)[2:26]
print(modules)

#loop through all modules, performing crude qgcomp, allmetals
for (j in 1:length(modules)){
  module <- modules[[j]]
  print(module)
    
    crude <- qgcomp.noboot(eval(parse(text = paste0(module,"~."))), expnms=Xnm,
                            dat=eval(parse(text = paste0("full[,c(Xnm,'",module,"')]"))),family=gaussian(),q=4)
    crude
    plot(crude)
    assign(paste0(module,".crude"),crude) 
}


#loop through all modules, performing crude qgcomp, toxics
for (j in 1:length(modules)){
  module <- modules[[j]]
  print(module)
  
  crude <- qgcomp.noboot(eval(parse(text = paste0(module,"~."))), expnms=toxicXnm,
                         dat=eval(parse(text = paste0("full[,c(toxicXnm,'",module,"')]"))),family=gaussian(),q=4)
  crude
  plot(crude)
  assign(paste0(module,".crude.toxic"),crude) 
}

#loop through all modules, performing crude qgcomp, essentials
for (j in 1:length(modules)){
  module <- modules[[j]]
  print(module)
  
  crude <- qgcomp.noboot(eval(parse(text = paste0(module,"~."))), expnms=essentialXnm,
                         dat=eval(parse(text = paste0("full[,c(essentialXnm,'",module,"')]"))),family=gaussian(),q=4)
  crude
  plot(crude)
  assign(paste0(module,".crude.essential"),crude) 
}


#loop through all modules, performing adjusted qgcomp, allmetals 
for (j in 1:length(modules)){
  module <- modules[[j]]
  print(module)
  
  adjusted <- qgcomp.noboot(eval(parse(text = paste0(module,"~."))), expnms=Xnm,
                         dat=eval(parse(text = paste0("full[,c(Xnm,covars,'",module,"')]"))),family=gaussian(),q=4)
  adjusted
  plot(adjusted)
  assign(paste0(module,".adjusted"),adjusted) 
}

#loop through all modules, performing adjusted qgcomp, toxics
for (j in 1:length(modules)){
  module <- modules[[j]]
  print(module)
  
  adjusted <- qgcomp.noboot(eval(parse(text = paste0(module,"~."))), expnms=toxicXnm,
                            dat=eval(parse(text = paste0("full[,c(toxicXnm,covars,'",module,"')]"))),family=gaussian(),q=4)
  adjusted
  plot(adjusted)
  assign(paste0(module,".adjusted.toxic"),adjusted) 
}

#loop through all modules, performing adjusted qgcomp, essential
for (j in 1:length(modules)){
  module <- modules[[j]]
  print(module)
  
  adjusted <- qgcomp.noboot(eval(parse(text = paste0(module,"~."))), expnms=essentialXnm,
                            dat=eval(parse(text = paste0("full[,c(essentialXnm,covars,'",module,"')]"))),family=gaussian(),q=4)
  adjusted
  plot(adjusted)
  assign(paste0(module,".adjusted.essential"),adjusted) 
}

#################################################################################################
#################################################################################################
#### Exporting results   
#################################################################################################
#################################################################################################
crudemodels <- paste0(modules,".crude")
crudemodelstoxic <- paste0(modules,".crude.toxic")
crudemodelsessential <- paste0(modules,".crude.essential")
adjustedmodels <-paste0(modules,".adjusted")
adjustedmodelstoxic <- paste0(modules,".adjusted.toxic")
adjustedmodelsessential <- paste0(modules,".adjusted.essential")
allmodels <- c(crudemodels, crudemodelstoxic, crudemodelsessential, 
               adjustedmodels, adjustedmodelstoxic, adjustedmodelsessential)

# SLOPE PARAMETERS
clean_print <- function(x){
  output  = data.frame(
    x$coef,
    sqrt(x$var.coef),
    x$ci.coef,
    x$pval
  )
  names(output) = c("Estimate", "Std. Error", "Lower CI", "Upper CI", "p value")
  return(output)
}

Results_SlopeParams <- data.frame() #empty vector to append dfs to
for (i in allmodels){
  print(i)
  df <- eval(parse(text = paste0("clean_print(",i,")"))) %>%
    rownames_to_column("Parameter") %>%
    mutate("Model" = i) 
  Results_SlopeParams <- rbind(Results_SlopeParams,df)
}

# METAL COEFFICIENTS
Results_MetalCoeffs <- data.frame()
for (i in allmodels){
  print(i)
  df <- eval(parse(text = paste0("as.data.frame(summary(",i,"$fit)$coefficients[,])"))) %>% 
    mutate("Model" = i)
  #Create variable containing metal identities
  df$Metal <- eval(parse(text = paste0("row.names(summary(",i,"$fit)$coefficients)")))
  df <- df %>% dplyr::select(c("Metal",1:5))
  Results_MetalCoeffs<- rbind(Results_MetalCoeffs,df)
}


# WEIGHTS
# Organizing the weights into one dataframe:
Results_MetalWeights <- data.frame()
for (i in allmodels){
  Results_PWeights <- eval(parse(text = paste0("as.data.frame(",i,"$pos.weights)"))) %>%
    rownames_to_column("Metal") %>%
    dplyr::rename("Weight" = 2) %>%
    mutate("Weight Direction" = "Positive")
  Results_NWeights <- eval(parse(text = paste0("as.data.frame(",i,"$neg.weights)"))) %>%
    rownames_to_column("Metal") %>%
    dplyr::rename("Weight" = 2) %>%
    mutate("Weight Direction" = "Negative")
  Results_Weights <- rbind(Results_PWeights, Results_NWeights) %>%
    mutate("Model" = i) %>% as.data.frame()
  #Merge
  Results_MetalWeights <- rbind(Results_MetalWeights, Results_Weights)
}

write.csv(Results_SlopeParams, paste0(Output_Folder,"/", cur_date, "_qgcomp_Results_SlopeParams.csv"), row.names=TRUE)
write.csv(Results_MetalCoeffs, paste0(Output_Folder,"/", cur_date, "_qgcomp_Results_MetalCoeffs.csv"), row.names=TRUE)
write.csv(Results_MetalWeights, paste0(Output_Folder,"/", cur_date, "_qgcomp_Results_MetalWeights.csv"), row.names=TRUE)

plot(ME2.crude)


#################################################################################################
#################################################################################################
#### Mixtures modeling with PCs   
#################################################################################################
#################################################################################################

PCs <- c("PC1","PC2") 
print(modules)

#crude models
results_pc <- data.frame()
for (i in 1:length(PCs)) {
  pc <- PCs[[i]]
  pc <- paste0(as.name(pc))
  print(pc)
  
  for (j in 1:length(modules)){
    module <- modules[[j]]
    print(module)
    
    output.lm.1 <- summary(lm(eval(parse(text = paste0(module,"~",pc))), data=full))$coefficients[2,1:4] %>%as.data.frame()
    colnames(output.lm.1)[1] <- paste0(pc,".",module)
    output.lm.1 <- output.lm.1 %>%  as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column(var="PCModule")
    results_pc<-rbind(results_pc,output.lm.1)
  }
}
results_pc <- results_pc %>% 
  dplyr::rename("p.value" ="Pr(>|t|)") %>% 
  dplyr::rename("t.value" ="t value") %>% 
  mutate(sigchange.pvalue = ifelse(p.value <0.05,1,0)) %>% 
  mutate(adj.p = p.adjust(p.value, method="BH", n=50)) %>% 
  mutate(sigchange.adjpvalue = ifelse(adj.p <0.05,1,0))


#adjusted models
results_pc_adj <- data.frame()
for (i in 1:length(PCs)) {
  pc <- PCs[[i]]
  pc <- paste0(as.name(pc))
  print(pc)
  
  for (j in 1:length(modules)){
    module <- modules[[j]]
    print(module)
    
    output.lm.1 <- summary(lm(eval(parse(text = paste0(module,"~",pc,"+ score + bmicat + smoke + sex + SV1 + SV2"))), data=full))$coefficients[2,1:4] %>%as.data.frame()
    colnames(output.lm.1)[1] <- paste0(pc,".",module)
    output.lm.1 <- output.lm.1 %>%  as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column(var="PCModule")
    results_pc_adj<-rbind(results_pc_adj,output.lm.1)
  }
}
results_pc_adj <- results_pc_adj %>% 
  dplyr::rename("p.value" ="Pr(>|t|)") %>% 
  dplyr::rename("t.value" ="t value") %>% 
  mutate(sigchange.pvalue = ifelse(p.value <0.05,1,0)) %>% 
  mutate(adj.p = p.adjust(p.value, method="BH", n=50)) %>% 
  mutate(sigchange.adjpvalue = ifelse(adj.p <0.05,1,0))

write.csv(results_pc, paste0(Output_Folder,"/", cur_date, "_PCA_WGCNA_Results_crude.csv"), row.names=TRUE)
write.csv(results_pc_adj, paste0(Output_Folder,"/", cur_date, "_PCA_WGCNA_Results_adjusted.csv"), row.names=TRUE)

#################################################################################################
#################################################################################################
#### Mixtures modeling with index vars
#################################################################################################
#################################################################################################

indexvars <- c("toxicindex","essentialindex") 
print(modules)

#crude models
results_index<- data.frame()
for (i in 1:length(indexvars)) {
  index <- indexvars[[i]]
  index <- paste0(as.name(index))
  print(index)
  
  for (j in 1:length(modules)){
    module <- modules[[j]]
    print(module)
    
    output.lm.1 <- summary(lm(eval(parse(text = paste0(module,"~",index))), data=full))$coefficients[2,1:4] %>%as.data.frame()
    colnames(output.lm.1)[1] <- paste0(index,".",module)
    output.lm.1 <- output.lm.1 %>%  as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column(var="IndexModule")
    results_index<-rbind(results_index,output.lm.1)
  }
}
results_index <- results_index %>% 
  dplyr::rename("p.value" ="Pr(>|t|)") %>% 
  dplyr::rename("t.value" ="t value") %>% 
  mutate(sigchange.pvalue = ifelse(p.value <0.05,1,0)) %>% 
  mutate(adj.p = p.adjust(p.value, method="BH", n=50)) %>% 
  mutate(sigchange.adjpvalue = ifelse(adj.p <0.05,1,0))

#adjusted models
results_index_adj <- data.frame()
for (i in 1:length(indexvars)) {
  index <- indexvars[[i]]
  index <- paste0(as.name(index))
  print(index)
  
  for (j in 1:length(modules)){
    module <- modules[[j]]
    print(module)
    
    output.lm.1 <- summary(lm(eval(parse(text = paste0(module,"~",index,"+ score + bmicat + smoke + sex +SV1 +SV2"))), data=full))$coefficients[2,1:4] %>%as.data.frame()
    colnames(output.lm.1)[1] <- paste0(index,".",module)
    output.lm.1 <- output.lm.1 %>%  as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column(var="IndexModule")
    results_index_adj<-rbind(results_index_adj,output.lm.1)
  }
}
results_index_adj <- results_index_adj %>% 
  dplyr::rename("p.value" ="Pr(>|t|)") %>% 
  dplyr::rename("t.value" ="t value") %>% 
  mutate(sigchange.pvalue = ifelse(p.value <0.05,1,0)) %>% 
  mutate(adj.p = p.adjust(p.value, method="BH", n=50)) %>% 
  mutate(sigchange.adjpvalue = ifelse(adj.p <0.05,1,0))

write.csv(results_index, paste0(Output_Folder,"/", cur_date, "_Index_WGCNA_Results_crude.csv"), row.names=TRUE)
write.csv(results_index_adj, paste0(Output_Folder,"/", cur_date, "_Index_WGCNA_Results_adjusted.csv"), row.names=TRUE)

