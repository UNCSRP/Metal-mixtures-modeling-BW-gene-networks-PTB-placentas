#################################################################################################
#################################################################################################
#### ELGAN cord metals gene expression: evaluating placental transcriptomic responses to cord 
#### metal levels 
#### Part 3: WGCNA analysis. A) Generating modules  
####
#### Code drafted by Lauren Eaves with refernce to code orginially written by Julia Rager 
#### Lasted updated: July 14th 2021, October 13th 2021
#################################################################################################
#################################################################################################
# Clean your working environment
rm(list=ls())

#################################################################################################
#################################################################################################
#### Installing and activating packages, and setting your working directory
#################################################################################################
#################################################################################################

#The WGCNA package requires the following packages to be installed: stats, fields, impute, grDevices, dynamicTreeCut (1.20 or higher), qvalue, utils, and flashClust.
#install.packages(c("fields", "impute", "dynamicTreeCut", "qvalue", "flashClust", "Hmisc") )

# If this is your first time running WGCNA, install it using:
#BiocManager::install(c("AnnotationDbi", "impute", "GO.db" , "preprocessCore"))
#BiocManager::install("org.Hs.eg.db") #this is the human specific annotation package for running the functional enrichment analysis 
#install.packages("WGCNA")

# Activate the appropriate packages:
library(tidyverse)
library(WGCNA)


# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Set your working directory to the folderpath containing your input files:
setwd("#yourwd")

# Create an output folder
Output_Folder <- ("youroutputfolder")

# Create a current date variable to name outputfiles
cur_date <- str_replace_all(Sys.Date(),"-","")


#################################################################################################
#################################################################################################
#### Loading, organizing, and initial visualizations of data
####
#################################################################################################
#################################################################################################

# Read in the transcriptomic/ mRNA expression data 
#Note that these data have already been QA/QC, filtered to subjects for which there is metal data, normalized in script: 1_genome-wide_DeSeq2/20210714_ELGAN_cordmetals_singlemetal_crudemodels_fullcohort_DeSeq2.R
countdata <- read.csv(file = '#countdata', header = TRUE)
dim(countdata)
#224 samples (plus Gene ID column), 11402 genes
# visualize these data quickly by viewing top left corner:
countdata[1:3,1:6]
colnames(countdata)[1] <- "Gene"
countdata <- countdata %>% column_to_rownames(var="Gene")

# Read in metadata (subject info/covariates and metals data)
subjectinfo <- read.csv("#input data file generated from script 2_ELGANcordmetals_multimetals_DESeq2 names DATE_SubjectInfo_SubjectsIncludedinModel", check.names=FALSE)
dim(subjectinfo)
#224 samples, 87 variables 
# Visualize these data quickly by viewing top left corner:
subjectinfo[1:3,1:10]
colnames(subjectinfo)

#Create a dataset just of metals and ID
metals <- subjectinfo %>% 
  dplyr::select(ID,contains("_ngg"),contains("_ugg"),-contains("log")) %>% 
  column_to_rownames(var="ID") %>% 
  mutate_all(as.numeric)
logmetals <- subjectinfo %>% 
  dplyr::select(ID,contains("log")) %>% 
  column_to_rownames(var="ID") %>% 
  mutate_all(as.numeric)
# Double checking that the same IDs are in the two datasets 
setequal(as.character(subjectinfo$ID), as.character(colnames(countdata)))

# Create a transposed dataframe of the mRNA data for future steps
t_countdata=as.data.frame(t(countdata))

#################################################################################################
#################################################################################################
## Further data cleaning (most of the time, this is not necessary for our data, since we do so much pre-processing)
## Modified based on WGCNA tutorial I.1. Data input and cleaning
#################################################################################################
#################################################################################################

# Check for genes and samples with too many missing values
# If the last statement returns TRUE, all genes have passed the cuts.
# If not, we can remove the offending genes and samples from the data
gsg = goodSamplesGenes(t_countdata, verbose = 3);
gsg$allOK
#this returned a TRUE, meaning all good genes (so didn't need the next code, lines 88-97)

#if (!gsg$allOK)
#{
#  # Optionally, print the gene and sample names that were removed:
#  if (sum(!gsg$goodGenes)>0)
#    printFlush(paste("Removing genes:", paste(names(t_countdata)[!gsg$goodGenes], collapse = ", ")));
#  if (sum(!gsg$goodSamples)>0)
#    printFlush(paste("Removing samples:", paste(rownames(t_countdata)[!gsg$goodSamples], collapse = ", ")));
#  # Remove the offending genes and samples from the data:
#  t_countdata = t_countdata[gsg$goodSamples, gsg$goodGenes]
#}


#################################################################################################
#################################################################################################
## Initial visualization of data
## Specifically making a dendogram of how trait (metals) data relate to gene data
#################################################################################################
#################################################################################################


# Calculating sample clusters
sampleTree = hclust(dist(t_countdata), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
metalColors = numbers2colors(metals, signed = FALSE);
metalColorslog = numbers2colors(logmetals, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
tiff(file = (paste0(Output_Folder,"/", cur_date, "_Dendogram_heatmap.tiff")), width = 18, height = 8, units = "in", pointsize = 12, res = 600)
plotDendroAndColors(sampleTree,
                    metalColors,
                    groupLabels = names(metals),
                    addGuide = FALSE,
                    guideAll = FALSE, 
                    guideCount = 50,
                    guideHang = 0.5, 
                    cex.colorLabels = 0.6,
                    cex.dendroLabels = 0.6, 
                    cex.rowText = 0.6,
                    main = "Sample dendrogram and metals heatmap")
dev.off()

tiff(file = (paste0(Output_Folder,"/", cur_date, "_Dendogram_heatmap_logmetals.tiff")), width = 18, height = 8, units = "in", pointsize = 12, res = 600)
plotDendroAndColors(sampleTree,
                    metalColorslog,
                    groupLabels = names(logmetals),
                    addGuide = FALSE,
                    guideAll = FALSE, 
                    guideCount = 50,
                    guideHang = 0.5, 
                    cex.colorLabels = 0.6,
                    cex.dendroLabels = 0.6, 
                    cex.rowText = 0.6,
                    main = "Sample dendrogram and log-transformed metals heatmap")
dev.off()


#################################################################################################
#################################################################################################
#### First step in WGCNA: Selecting a soft threshold
#### Code was developed (with modifications) based on WGCNA tutorial I.2c.1 
#### "Dealing with large data sets: block-wise network construction and module detection, Selecting soft threshold"
#################################################################################################
#################################################################################################

# Choosing the soft-thresholding power: analysis of network topology
# The function pickSoftThreshold is used to performs the analysis of network topology, and
# aids the user in choosing a proper soft-thresholding power

# Choose a set of soft-thresholding powers to test
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(t_countdata, powerVector = powers, verbose = 5)


# Plot the results
# First prep the plot window:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Then plot the scale-free topology fit index as a function of the soft-thresholding power:
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# Here we choose the power 14, which is the lowest power for which the scale-free topology fit index reaches 0.90.

#################################################################################################
#################################################################################################
#### Second step in WGCNA: Network construction and module detection
#### Code was developed (with modifications) based on WGCNA tutorial I.2a. 
#### "Automatic network construction and module detection"
#################################################################################################
#################################################################################################


# Calculating module eigengenes 
net = blockwiseModules(t_countdata,
                       power = 14,
                       TOMType = "unsigned",
                       minModuleSize = 30,
                       mergeCutHeight = 0.1,
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "GenesTOM",
                       verbose = 3)
# TOMType should be "signed" when you want to include genes/chemicals that are co-modulated in the same direction in each module


# Plot the results by first prepping the plot window
sizeGrWindow(12, 9)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath
tiff(file = (paste0(Output_Folder,"/", cur_date, "_Dendogram_moduleassignments.tiff")), width = 18, height = 8, units = "in", pointsize = 12, res = 600)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
table(moduleLabels)
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

#################################################################################################
#################################################################################################
# Organize and save these module eigengene assignments, to analyze/generate future visualizations:
#################################################################################################
#################################################################################################

write.csv(MEs,paste0(Output_Folder,"/", cur_date,"_MEs_AggregateValues_bysubject.csv"))

moduleLabels_df <- as.data.frame(moduleLabels)
moduleColors_df <- as.data.frame(moduleColors)
GeneModuleAssignments <- cbind(moduleLabels_df, moduleColors_df)
  
write.csv(GeneModuleAssignments,paste0(Output_Folder,"/", cur_date,"_Gene_to_Module_Assignments.csv"))

#################################################################################################
#################################################################################################
#### Third step in WGCNA: Relating modules to metal information and ranking genes of most importance in each module
#### Code was developed (with modifications) based on WGCNA tutorial I.3.
#### "Relating modules to external information and identifying important genes"
#################################################################################################
#################################################################################################

###########
# Quantifying module-metal associations
###########

# Define the numbers of genes and samples
nGenes = ncol(t_countdata);
nSamples = nrow(t_countdata);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(t_countdata, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, logmetals, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# Plot the results by first prepping the plot window
sizeGrWindow(10,6)

# We will display correlation R values and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));


# Display the correlation values within a heatmap plot -- FIX SO THAT IT PRINTS FULL LABELS 
tiff(file = (paste0(Output_Folder,"/", cur_date, "_Module_Metal_Heatmap.tiff")), width = 12, height = 8, units = "in", pointsize = 12, res = 600)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(logmetals),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-Metal relationships"))
dev.off()                           

heatmap <- labeledHeatmap(Matrix = moduleTraitCor,
                          xLabels = names(logmetals),
                          yLabels = names(MEs),
                          ySymbols = names(MEs),
                          colorLabels = FALSE,
                          colors = blueWhiteRed(50),
                          textMatrix = textMatrix,
                          setStdMargins = FALSE,
                          cex.text = 0.7,
                          zlim = c(-1,1),
                          main = paste("Module-Metal relationships"))

heatmap

#################################################################################################
#################################################################################################
# Quantifying gene-module-trait associations through two values:
# (1) Gene Significance (GS)  
# (2) Module Membership (MM)
# Code developed from WGCNA tutorial step 3.b. Gene relationship to trait and important modules: 
# "Gene Signicance and Module Membership"
# And 3.c. Intramodular analysis: 
# "Identifying genes with high gene significance (GS) and module membership (MM)"
#################################################################################################
#################################################################################################

# Define variable weight containing the weight column of subject information / trait data
Pb_ngg_log = as.data.frame(logmetals$Pb_ngg_log);
names(Pb_ngg_log) = "Pb_ngg_log"
Hg_ngg_log = as.data.frame(logmetals$Hg_ngg_log);
names(Hg_ngg_log) = "Hg_ngg_log"
Cd_ngg_log = as.data.frame(logmetals$Cd_ngg_log);
names(Cd_ngg_log) = "Cd_ngg_log"
As_ngg_log = as.data.frame(logmetals$As_ngg_log);
names(As_ngg_log) = "As_ngg_log"
Zn_ugg_log = as.data.frame(logmetals$Zn_ugg_log);
names(Zn_ugg_log) = "Zn_ugg_log"
Cu_ugg_log = as.data.frame(logmetals$Cu_ugg_log);
names(Cu_ugg_log) = "Cu_ugg_log"
Mn_ugg_log = as.data.frame(logmetals$Mn_ugg_log);
names(Mn_ugg_log) = "Mn_ugg_log"

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(t_countdata, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance_Pb = as.data.frame(cor(t_countdata, Pb_ngg_log, use = "p"));
GSPvalue_Pb = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_Pb), nSamples));
names(geneTraitSignificance_Pb) = paste("GS.", names(Pb_ngg_log), sep="");
names(GSPvalue_Pb) = paste("p.GS.", names(Pb_ngg_log), sep="");

geneTraitSignificance_Hg = as.data.frame(cor(t_countdata, Hg_ngg_log, use = "p"));
GSPvalue_Hg = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_Hg), nSamples));
names(geneTraitSignificance_Hg) = paste("GS.", names(Hg_ngg_log), sep="");
names(GSPvalue_Hg) = paste("p.GS.", names(Hg_ngg_log), sep="");

geneTraitSignificance_Cd = as.data.frame(cor(t_countdata, Cd_ngg_log, use = "p"));
GSPvalue_Cd = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_Cd), nSamples));
names(geneTraitSignificance_Cd) = paste("GS.", names(Cd_ngg_log), sep="");
names(GSPvalue_Cd) = paste("p.GS.", names(Cd_ngg_log), sep="");

geneTraitSignificance_As = as.data.frame(cor(t_countdata, As_ngg_log, use = "p"));
GSPvalue_As = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_As), nSamples));
names(geneTraitSignificance_As) = paste("GS.", names(As_ngg_log), seo="");
names(GSPvalue_As) = paste("p.GS.", names(As_ngg_log), sep="");

geneTraitSignificance_Zn = as.data.frame(cor(t_countdata, Zn_ugg_log, use = "p"));
GSPvalue_Zn = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_Zn), nSamples));
names(geneTraitSignificance_Zn) = paste("GS.", names(Zn_ugg_log), seo="");
names(GSPvalue_Zn) = paste("p.GS.", names(Zn_ugg_log), sep="");

geneTraitSignificance_Cu = as.data.frame(cor(t_countdata, Cu_ugg_log, use = "p"));
GSPvalue_Cu = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_Cu), nSamples));
names(geneTraitSignificance_Cu) = paste("GS.", names(Cu_ugg_log), seo="");
names(GSPvalue_Cu) = paste("p.GS.", names(Cu_ugg_log), sep="");

geneTraitSignificance_Mn = as.data.frame(cor(t_countdata, Mn_ugg_log, use = "p"));
GSPvalue_Mn = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_Mn), nSamples));
names(geneTraitSignificance_Mn) = paste("GS.", names(Mn_ugg_log), seo="");
names(GSPvalue_Mn) = paste("p.GS.", names(Mn_ugg_log), sep="");

# Identifying genes with high gene significance (GS) and module membership (MM)
# Here, focusing on MEblue + Pb 

#module = "blue"
#column = match(module, modNames);
#moduleGenes = moduleColors==module;
#sizeGrWindow(7, 7);
#par(mfrow = c(1,1));
#verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                   abs(geneTraitSignificance[moduleGenes, 1]),
#                   xlab = paste("Module Membership in", module, "module"),
#                   ylab = "Gene significance for Pb_ngg_log",
#                   main = paste("Module membership vs. gene significance\n"),
#                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



###########
# Organizing and summarizing output
###########

# Create the starting data frame
IDcol = names(t_countdata)
geneInfo0 = data.frame(ID = IDcol,
                       moduleColor = moduleColors,
                       geneTraitSignificance_Pb,
                       GSPvalue_Pb,
                       geneTraitSignificance_Hg,
                       GSPvalue_Hg,
                       geneTraitSignificance_Cd,
                       GSPvalue_Cd,
                       geneTraitSignificance_As,
                       GSPvalue_As,
                       geneTraitSignificance_Zn,
                       GSPvalue_Zn,
                       geneTraitSignificance_Cu,
                       GSPvalue_Cu,
                       geneTraitSignificance_Mn,
                       GSPvalue_Mn
                       )

# Order modules by their significance for Pb_ngg_log
modOrder = order(-abs(cor(MEs, Pb_ngg_log, use = "p")));

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Pb_ngg_log));
geneInfo = geneInfo0[geneOrder, ]

# Export gene significance (GS) and modular membership (MM) for every gene
# Note that this also includes all genes in every module
# I QCed that the same results occur when focusing on other MEs in the above code - the same data are organized/exported; just reordered

# Before exporting, let's merge with module membership numbers
ModuleMembership <- merge(GeneModuleAssignments, geneInfo, by=0, all.x=TRUE, all.y=TRUE)
# removing some duplicate columns
ModuleMembership$ID <- NULL
ModuleMembership$moduleColor <- NULL
ModuleMembership <- ModuleMembership %>% dplyr::rename("Gene"="Row.names")

write.csv(ModuleMembership, file = paste0(Output_Folder,"/", cur_date,"_GeneSignificance_ModuleMembership_permetal.csv"), row.names=FALSE)

#################################################################################################
#################################################################################################
# Interfacing network analysis with functional annotation and gene ontology 
# Code developed from AnRichment tutorial:https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/Tutorials/anRichment-Tutorial1.pdf
#################################################################################################
#################################################################################################

#practicing based on tutorial: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/Tutorials/anRichment-Tutorial1.pdf
#source("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/installAnRichment.R");
#installAnRichment();
library("anRichment");

#need to generate a vector of the genes and a vector of their corresponding module colors 
symbol.0 = geneInfo0$ID;
moduleColor = geneInfo$moduleColor;
table(moduleColor) #note that grey is the "non-module" module, where genes that were not assigned to any module are stored 

# Some gene symbols have the form "XYZ /// ABC". Keep only the first symbol of all such multi-symbols.
split = strsplit(symbol.0, split = " /// ", fixed = TRUE);
symbol = sapply(split, function(x) x[1]);
# Convert symbols to Entrez IDs
entrez = convert2entrez(organism = "human", symbol = symbol);
# How many conversions were successful?
table(is.finite(entrez))

#construct GO library and run enrichment analysis
GOcollection = buildGOcollection(organism = "human")
GOenrichment = enrichmentAnalysis(
  classLabels = moduleColor, identifiers = entrez,
  refCollection = GOcollection,
  useBackground = "given",
  threshold = 1e-4,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "grey") #ignore this module, ie. a non-module 

table.display = GOenrichment$enrichmentTable
table.display$overlapGenes = shortenStrings(table.display$overlapGenes, maxLength = 70,
                                            split = "|")
head(table.display)

write.csv(GOenrichment$enrichmentTable, paste0(Output_Folder,"/", cur_date,"_GOenrichment_enrichmentTable.csv"),
          row.names = FALSE)
