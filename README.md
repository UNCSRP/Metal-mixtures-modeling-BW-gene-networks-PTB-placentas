# Metal-mixtures-modeling-BW-gene-networks-PTB-placentas
Code and data for manuscript: Metals mixtures modeling identifies birth weight-associated gene networks in the placentas of children born extremely preterm, which was published in Chemosphere in 2024 (PMID: 36493891). doi: 10.1016/j.chemosphere.2022.137469

# Abstract 
Prenatal exposure to toxic metals is linked to numerous adverse birth and later-in-life outcomes. These outcomes are tied to disrupted biological processes in fetal-derived tissues such as the placenta and umbilical cord yet these pathways are understudied in these target tissues. In the present study, we set out to examine the relationship between metal concentrations in umbilical cord and altered gene expression networks in placental tissue. These novel relationships were investigated in a subset of the Extremely Low Gestational Age Newborn (ELGAN) cohort (n=226). Prenatal exposure to 11 metals/metalloids was measured using inductively coupled plasma tandem-mass spectrometry (ICP-MS/MS) in cord tissue, ensuring passage through the placental barrier. RNA-sequencing was used to quantify >37,000 mRNA transcripts. Differentially expressed genes (DEGs) were identified with respect to each metal. Weighted gene co-expression analysis identified gene networks modulated by metals. Two innovative mixtures modeling techniques - principal components analysis as well as quantile-based g-computation - were employed to identify genes/gene networks associated with multi-metal exposure. Individually, lead was associated with the strongest genomic response of 191 DEGs. 657 DEGs were related to joint lead and cadmium exposure, including DNA Methyl Transferase 1 (DNMT1). These genes were enriched for the Eukaryotic Initiation Factor 2 (EIF2) pathway. Four gene networks, each containing genes within a NF-kB-mediated network, were significantly increased in average expression level when increasing all metal concentrations by one quartile. All four of these metal mixture-associated gene networks were negatively correlated with important predictors of neonatal health including birth weight, placenta weight and fetal growth. Bringing together novel methodologies from epidemiological mixtures analyses and toxicogenomics, applied to a unique cohort of extremely preterm children, the present study highlighted critical genes and pathways in the placenta dysregulated by prenatal metal mixtures. These represent potential underlying mechanisms of effect in the developmental origins of metal-induced disease.  

# Files 
Data: 
1. ELGAN_cohortdata_cordmetals_surrogatevars_imputationcomplete: ELGAN cord metal concentrations, clinical and demographic data with missing values imputed using random forest modeling, surrogate variables to capture unwanted variance generated using svaseq 
2. ELGAN_MasterIDs_file: ID file to connect mRNA count data and clinical/demographic data 
3. ELGAN_neonatalphenotypes: extra ELGAN clinical data including four neonatal phenotypes 

Code: 
1. 1_ELGANcordmetals_singlemetals_DESeq2: identifying DEGs in relation to individual metals 
2. 2_ELGANcordmetals_multimetals_DESeq2: identifying DEGs in relation to metal mixtures (principal components) 
3. 3_ELGANcordmetals_WGCNA_generatemodules: generating gene networks using WGCNA
4. 4_ELGANcordmetals_WGCNA_singlemetals_regressions: identifying gene networks associated with individual metals
5. 5_ELGANcordmetals_WGCNA_multimetals_qgcomp: identifying gene networks associated with metal mixtures (quantile based g-computation)
6. 6_ELGANcordmetals_connectiontopheno: evaluating correlations between metal mixture responsive genes and gene networks and neonatal phenotypes 
