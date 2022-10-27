#################################################################################################
#################################################################################################
#### ELGAN cord metals gene expression: evaluating placental transcriptomic responses to cord 
#### metal levels 
#### Part 4:Evaluating placental transcriptomic responses to cord metal levels 
#### Connection to key phenotypes: bw (birthweight), placwt (placental weight), gadays (gestational age), and z (fetal growth)
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
library(ggpubr)
library(ggplot2)
library(Hmisc)
library(corrplot)


# Set your working directory to the folderpath containing your input files:
setwd("#yourwd")

# Create an output folder
Output_Folder <- ("#youroutputfolder")

# Create a current date variable to name outputfiles
cur_date <- paste0(str_replace_all(Sys.Date(),"-",""),"allmetals")



#################################################################################################
#################################################################################################
#### Load data and merge
#################################################################################################
#################################################################################################

mes<- read.csv(file="#input file generated from script 3_ELGANcordmetals_WGCNA_generatemodules, titled DATE_MEs_AggregateValues_bysubject.csv") #module values for each subject 
ids <- read.csv(file="idsfile")
subjectinfo <-read.csv(file="#ELGAN_neonatalphenotypes")

subjectinfo <- subjectinfo %>% 
  select(Elgan_ID=id, bw, gadays, placwt) %>% 
  distinct(Elgan_ID, .keep_all=TRUE)

ids <- ids %>% 
  select(Elgan_ID,Fry_ID=mRNA.seq_datafile_ID)

mes <- mes %>% 
  select(Fry_ID=X,ME2,ME6,ME13,ME14)  
  
subjectinfo <- left_join(subjectinfo,ids, by="Elgan_ID") %>% 
  filter(!is.na(Fry_ID))

data <- left_join(subjectinfo,mes,by="Fry_ID") %>% 
  filter(!is.na(ME2))

#################################################################################################
#################################################################################################
#### Look at distributions of variables
#################################################################################################
#################################################################################################

hist(data$ME2)
hist(data$ME6)
hist(data$ME13)
hist(data$ME14)

summary(data$ME2)
summary(data$ME6)
summary(data$ME13)
summary(data$ME14)

summary(data$bw)
summary(data$placwt)
summary(data$z)
summary(data$gadays)

#################################################################################################
#################################################################################################
#### Correlations between continuous and MEs
#################################################################################################
#################################################################################################

#spearman correlations
correlationset <- data %>% select(bw,placwt,z,gadays,ME2,ME6,ME13,ME14) %>% 
  as.matrix()

spearman <-rcorr(correlationset, type=c("spearman"))
#access correlation matrix 
spearman_corrs <- as.matrix(spearman[[1]])[1:4,5:8]
#access p values 
spearman_p <-as.matrix(spearman[[3]])[1:4,5:8]

spearman_corrs <- as.data.frame(t(spearman_corrs)) %>% 
  dplyr::rename("Birth weight" = "bw", 
                "Placenta weight"="placwt",
                "Fetal growth"="z",
                "Gestational Age"="gadays") 
spearman_corrs <- as.data.frame(t(spearman_corrs)) %>% as.matrix()

spearman_p <- as.data.frame(t(spearman_p)) %>% 
  dplyr::rename("Birth weight" = "bw", 
                "Placenta weight"="placwt",
                "Fetal growth"="z",
                "Gestational Age"="gadays")  
spearman_p <- as.data.frame(t(spearman_p)) %>% as.matrix()

png(file = (paste0(Output_Folder,"/", cur_date, "_corrplot_MEs_contvars.png")), width = 4, height = 4, units = "in", pointsize = 8, res = 300)
corrplot(spearman_corrs, method="color", tl.col="black",
         tl.cex = 1, tl.srt=45, insig = 'label_sig', sig.level =0.05, p.mat=spearman_p)
dev.off()

#reorganize data into one dataframe
coeff <- as.data.frame(spearman_corrs) %>% 
  rownames_to_column(var="Phenotype") %>%
  pivot_longer(2:5,names_to="ME") %>% 
  dplyr::rename("Spearman_correlation"="value")

pvalues <- as.data.frame(spearman_p) %>% 
  rownames_to_column(var="Phenotype") %>%
  pivot_longer(2:5,names_to="ME") %>% 
  dplyr::rename("p_value"="value")

joined <- inner_join(coeff, pvalues, by=c("Phenotype", "ME")) 
write.csv(joined, paste0(Output_Folder,"/", cur_date, "_correlations_MEs_contvars.csv"), row.names= TRUE)

#################################################################################################
#################################################################################################
#### Correlations between continuous and MEs: scatterplots
#################################################################################################
#################################################################################################


#Scatterplots with birth weight
plot1 <-
  ggplot(data, aes(x=bw, y=ME2)) +
  geom_point(size=2, shape=23)+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ 
  labs(x="Birth weight (g)", y = "ME2") + 
  ylim(-0.032,0.1)+
  annotate("text", x=900, y=0.1, label = "Coefficient=-0.231, p-value <0.01", size=3)+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill="white", colour="black",
                                                                               size = 0.5, linetype = "solid"))+ 
  coord_cartesian(clip = "off")
plot1

plot2 <-
  ggplot(data, aes(x=bw, y=ME6)) +
  geom_point(size=2, shape=23)+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ 
  labs(x="Birth weight (g)", y = "ME6") + 
  ylim(-0.06,0.4)+
  annotate("text", x=900, y=0.4, label = "Coefficient=-0.255, p-value <0.01", size=3)+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill="white", colour="black",
                                                                               size = 0.5, linetype = "solid"))+ 
  coord_cartesian(clip = "off")
plot2


plot3 <-
  ggplot(data, aes(x=bw, y=ME13)) +
  geom_point(size=2, shape=23)+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ 
  labs(x="Birth weight (g)", y = "ME13") + 
  ylim(-0.03,0.2)+
  annotate("text", x=900, y=0.2, label = "Coefficient=-0.269, p-value <0.01", size=3)+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill="white", colour="black",
                                                                               size = 0.5, linetype = "solid"))+ 
  coord_cartesian(clip = "off")
plot3


plot4 <-
  ggplot(data, aes(x=bw, y=ME14)) +
  geom_point(size=2, shape=23)+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ 
  labs(x="Birth weight(g)", y = "ME14") + 
  ylim(-0.055,0.4)+
  annotate("text", x=900, y=0.4, label = "Coefficient=-0.231, p-value <0.01", size=3)+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill="white", colour="black",
                                                                               size = 0.5, linetype = "solid"))+ 
  coord_cartesian(clip = "off")
plot4


combinedbw <- ggarrange(plot1,plot2,plot3,plot4,
                      nrow=2,ncol=2)
plot(combinedbw)

png(file = (paste0(Output_Folder,"/", cur_date, "_scatterplots_MEs_bw.png")), width = 8, height = 6, units = "in", pointsize = 12, res = 600)
plot(combinedbw)
dev.off()





#Scatterplots with placenta weight
plot5 <-
  ggplot(data, aes(x=placwt, y=ME2)) +
  geom_point(size=2, shape=23)+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ 
  labs(x="Placenta weight (g)", y = "ME2") + 
  ylim(-0.032,0.1)+
  annotate("text", x=650, y=0.1, label = "Coefficient=-0.145, p-value=0.03", size=3)+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill="white", colour="black",
                                                                               size = 0.5, linetype = "solid"))+ 
  coord_cartesian(clip = "off")
plot5

plot6 <-
  ggplot(data, aes(x=placwt, y=ME6)) +
  geom_point(size=2, shape=23)+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ 
  labs(x="Placenta weight (g)", y = "ME6") + 
  ylim(-0.06,0.4)+
  annotate("text", x=650, y=0.4, label = "Coefficient=-0.152, p-value=0.03", size=3)+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill="white", colour="black",
                                                                               size = 0.5, linetype = "solid"))+ 
  coord_cartesian(clip = "off")
plot6


plot7 <-
  ggplot(data, aes(x=placwt, y=ME13)) +
  geom_point(size=2, shape=23)+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ 
  labs(x="Placenta weight (g)", y = "ME13") + 
  ylim(-0.03,0.2)+
  annotate("text", x=650, y=0.2, label = "Coefficient=-0.124, p-value=0.07", size=3)+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill="white", colour="black",
                                                                               size = 0.5, linetype = "solid"))+ 
  coord_cartesian(clip = "off")
plot7


plot8 <-
  ggplot(data, aes(x=placwt, y=ME14)) +
  geom_point(size=2, shape=23)+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ 
  labs(x="Placenta weight (g)", y = "ME14") + 
  ylim(-0.055,0.4)+
  annotate("text", x=650, y=0.4, label = "Coefficient=-0.167, p-value=0.01", size=3)+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill="white", colour="black",
                                                                               size = 0.5, linetype = "solid"))+ 
  coord_cartesian(clip = "off")
plot8



#Scatterplots with gestational age
plot9 <-
  ggplot(data, aes(x=gadays, y=ME2)) +
  geom_point(size=2, shape=23)+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ 
  labs(x="Gestational age (days)", y = "ME2") + 
  ylim(-0.032,0.1)+
  annotate("text", x=180, y=0.1, label = "Coefficient=-0.127, p-value=0.06", size=3)+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill="white", colour="black",
                                                                               size = 0.5, linetype = "solid"))+ 
  coord_cartesian(clip = "off")
plot9

plot10 <-
  ggplot(data, aes(x=gadays, y=ME6)) +
  geom_point(size=2, shape=23)+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ 
  labs(x="Gestational age (days)", y = "ME6") + 
  ylim(-0.06,0.4)+
  annotate("text", x=180, y=0.4, label = "Coefficient=-0.147, p-value=0.03", size=3)+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill="white", colour="black",
                                                                               size = 0.5, linetype = "solid"))+ 
  coord_cartesian(clip = "off")
plot10


plot11 <-
  ggplot(data, aes(x=gadays, y=ME13)) +
  geom_point(size=2, shape=23)+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ 
  labs(x="Gestational age (days)", y = "ME13") + 
  ylim(-0.03,0.2)+
  annotate("text", x=180, y=0.2, label = "Coefficient=-0.107, p-value=0.11", size=3)+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill="white", colour="black",
                                                                               size = 0.5, linetype = "solid"))+ 
  coord_cartesian(clip = "off")
plot11


plot12 <-
  ggplot(data, aes(x=gadays, y=ME14)) +
  geom_point(size=2, shape=23)+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ 
  labs(x="Gestational age (days)", y = "ME14") + 
  ylim(-0.055,0.4)+
  annotate("text", x=180, y=0.4, label = "Coefficient=-0.182, p-value=0.01", size=3)+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill="white", colour="black",
                                                                               size = 0.5, linetype = "solid"))+ 
  coord_cartesian(clip = "off")
plot12

#Scatterplots with fetal growth
plot13 <-
  ggplot(data, aes(x=z, y=ME2)) +
  geom_point(size=2, shape=23)+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ 
  labs(x="Fetal growth (BW by GA, z score)", y = "ME2") + 
  ylim(-0.032,0.1)+
  annotate("text", x=0, y=0.1, label = "Coefficient=-0.191, p-value<0.01", size=3)+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill="white", colour="black",
                                                                               size = 0.5, linetype = "solid"))+ 
  coord_cartesian(clip = "off")
plot13

plot14 <-
  ggplot(data, aes(x=z, y=ME6)) +
  geom_point(size=2, shape=23)+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ 
  labs(x="Fetal growth (BW by GA, z score)", y = "ME6") + 
  ylim(-0.06,0.4)+
  annotate("text", x=0, y=0.4, label = "Coefficient=-0.209, p-value<0.01", size=3)+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill="white", colour="black",
                                                                               size = 0.5, linetype = "solid"))+ 
  coord_cartesian(clip = "off")
plot14


plot15 <-
  ggplot(data, aes(x=z, y=ME13)) +
  geom_point(size=2, shape=23)+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ 
  labs(x="Fetal growth (BW by GA, z score)", y = "ME13") + 
  ylim(-0.03,0.2)+
  annotate("text", x=0, y=0.2, label = "Coefficient=-0.277, p-value<0.01", size=3)+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill="white", colour="black",
                                                                               size = 0.5, linetype = "solid"))+ 
  coord_cartesian(clip = "off")
plot15


plot16 <-
  ggplot(data, aes(x=z, y=ME14)) +
  geom_point(size=2, shape=23)+ 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+ 
  labs(x="Fetal growth (BW by GA, z score)", y = "ME14") + 
  ylim(-0.055,0.4)+
  annotate("text", x=0, y=0.4, label = "Coefficient=-0.212, p-value<0.01", size=3)+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill="white", colour="black",
                                                                               size = 0.5, linetype = "solid"))+ 
  coord_cartesian(clip = "off")
plot16


combined <- ggarrange(plot1,plot2,plot3,plot4,
                      plot5,plot6,plot7,plot8,
                      plot9,plot10,plot11,plot12,
                      plot13,plot14,plot15,plot16,
                      nrow=4,ncol=4)
plot(combined)

png(file = (paste0(Output_Folder,"/", cur_date, "_scatterplots_MEs_contvars.png")), width = 16, height = 12, units = "in", pointsize = 12, res = 300)
plot(combined)
dev.off()


#################################################################################################
#################################################################################################
#### Correlations between continuous and PC-genes 
#################################################################################################
#################################################################################################

normcounts <- read.csv(file="#input normalized counts generated from script 1_ELGANcordmetals_singlemetals_DESeq2, titled DATE_NormCounts_IncludedSubjects")
pc1genes <- read.csv(file = "#input results filed generated from script2_ELGANcordmetals_multimetals_DESeq2, titled DATE_PC1_AllStatResults.csv") %>% 
  filter(padj <0.1)
pc1genes <- pc1genes$X
pc2genes <- read.csv(file = "#input results filed generated from script2_ELGANcordmetals_multimetals_DESeq2, titled DATE_PC2_AllStatResults.csv")%>% 
  filter(padj <0.1)
pc2genes <- pc2genes$X

normcounts <- normcounts %>% 
  column_to_rownames(var="X")
normcounts <- t(normcounts)

pc1genescounts <- as.data.frame(normcounts) %>% select(all_of(c(pc1genes))) %>% 
  rownames_to_column(var="Fry_ID")
pc2genescounts <- as.data.frame(normcounts) %>% select(all_of(c(pc2genes)))%>% 
  rownames_to_column(var="Fry_ID")
pc1genescounts <- left_join(pc1genescounts,ids,by="Fry_ID")
pc2genescounts <- left_join(pc2genescounts,ids,by="Fry_ID")


correlationset_pc1 <- data %>% select(bw,placwt,z,gadays,Elgan_ID)  
correlationset_pc1 <- left_join(correlationset_pc1,pc1genescounts,by="Elgan_ID") %>% 
  select(-"Elgan_ID",-"Fry_ID") %>% as.matrix()

correlationset_pc2 <- data %>% select(bw,placwt,z,gadays,Elgan_ID)  
correlationset_pc2 <- left_join(correlationset_pc2,pc2genescounts,by="Elgan_ID")%>% 
  select(-"Elgan_ID",-"Fry_ID") %>% as.matrix()

#PC1 gene correlations 
spearman <-rcorr(correlationset_pc1, type=c("spearman"))
#access correlation matrix 
spearman_corrs <- as.matrix(spearman[[1]])[1:4,5:44]
#access p values 
spearman_p <-as.matrix(spearman[[3]])[1:4,5:44]

spearman_corrs <- as.data.frame(t(spearman_corrs)) %>% 
  dplyr::rename("Birth weight" = "bw", 
                "Placenta weight"="placwt",
                "Fetal growth"="z",
                "Gestational Age"="gadays") 
spearman_corrs <- as.data.frame(t(spearman_corrs)) %>% as.matrix()

spearman_p <- as.data.frame(t(spearman_p)) %>% 
  dplyr::rename("Birth weight" = "bw", 
                "Placenta weight"="placwt",
                "Fetal growth"="z",
                "Gestational Age"="gadays")  
spearman_p <- as.data.frame(t(spearman_p)) %>% as.matrix()

png(file = (paste0(Output_Folder,"/", cur_date, "_corrplot_PC1genes_contvars.png")), width = 15, height = 4, units = "in", pointsize = 8, res = 300)
corrplot(spearman_corrs, method="color", tl.col="black",
         tl.cex = 1, tl.srt=45, insig = 'label_sig', sig.level =0.05, p.mat=spearman_p)
dev.off()

#reorganize data into one dataframe
coeff <- as.data.frame(spearman_corrs) %>% 
  rownames_to_column(var="Phenotype") %>%
  pivot_longer(2:41,names_to="ME") %>% 
  dplyr::rename("Spearman_correlation"="value")

pvalues <- as.data.frame(spearman_p) %>% 
  rownames_to_column(var="Phenotype") %>%
  pivot_longer(2:41,names_to="ME") %>% 
  dplyr::rename("p_value"="value")

joined <- inner_join(coeff, pvalues, by=c("Phenotype", "ME")) 
write.csv(joined, paste0(Output_Folder,"/", cur_date, "_correlations_PC1genes_contvars.csv"), row.names= TRUE)



#PC2 gene correlations 
spearman <-rcorr(correlationset_pc2, type=c("spearman"))
#access correlation matrix 
spearman_corrs <- as.matrix(spearman[[1]])[1:4,5:661]
#access p values 
spearman_p <-as.matrix(spearman[[3]])[1:4,5:661]

spearman_corrs <- as.data.frame(t(spearman_corrs)) %>% 
  dplyr::rename("Birth weight" = "bw", 
                "Placenta weight"="placwt",
                "Fetal growth"="z",
                "Gestational Age"="gadays") 
spearman_corrs <- as.data.frame(t(spearman_corrs)) %>% as.matrix()

spearman_p <- as.data.frame(t(spearman_p)) %>% 
  dplyr::rename("Birth weight" = "bw", 
                "Placenta weight"="placwt",
                "Fetal growth"="z",
                "Gestational Age"="gadays")  
spearman_p <- as.data.frame(t(spearman_p)) %>% as.matrix()

png(file = (paste0(Output_Folder,"/", cur_date, "_corrplot_PC2genes_contvars.png")), width = 30, height = 4, units = "in", pointsize = 8, res = 300)
corrplot(spearman_corrs, method="color", tl.col="black",
         tl.cex = 1, tl.srt=45, insig = 'label_sig', sig.level =0.05, p.mat=spearman_p)
dev.off()

#reorganize data into one dataframe
coeff <- as.data.frame(spearman_corrs) %>% 
  rownames_to_column(var="Phenotype") %>%
  pivot_longer(2:658,names_to="ME") %>% 
  dplyr::rename("Spearman_correlation"="value")

pvalues <- as.data.frame(spearman_p) %>% 
  rownames_to_column(var="Phenotype") %>%
  pivot_longer(2:658,names_to="ME") %>% 
  dplyr::rename("p_value"="value")

joined <- inner_join(coeff, pvalues, by=c("Phenotype", "ME")) 
write.csv(joined, paste0(Output_Folder,"/", cur_date, "_correlations_PC2genes_contvars.csv"), row.names= TRUE)

