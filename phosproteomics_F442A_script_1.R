### Data analysis - 1, phosphoproteomics, F442A cells
# Script used for initial data analysis and normalization
# Niels Banhos Danneskiold-Samsoee, Sep 16 2022

.libPaths()
# clear workspace
rm(list=ls())

# Load necessary packages
library(Matrix)
library(readxl)
library(DEP)
library(dplyr)
library(SummarizedExperiment)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(purrr)
library(ggplot2)
library(rentrez)
library(ggvenn)
library(RColorBrewer)
library(car)

# set working directory
setwd("~/Isthmin phoshoproteomics")

# Loading raw data
raw_dta <- read_excel("20-2114_Phospho_Analyzed_Raw data from Northwest.xlsx", sheet = "All-Phospho")
#View(raw_dta)
colnames(raw_dta)

# Reformating column names
colnames(raw_dta) <- gsub(".","_", colnames(raw_dta), fixed=TRUE)
colnames(raw_dta) <- gsub(" ","_", colnames(raw_dta), fixed=TRUE)
colnames(raw_dta) <- gsub("-","_", colnames(raw_dta), fixed=TRUE)
colnames(raw_dta) <- gsub("(r\\d)_+1\\d","\\1_raw", colnames(raw_dta))
colnames(raw_dta) <- gsub("(r\\d)_+2\\d","\\1_norm", colnames(raw_dta))
colnames(raw_dta) <- gsub("\\(|\\)","", colnames(raw_dta))
grep("(r\\d)_+1",colnames(raw_dta))

# Count number of unique proteins
length(unique(raw_dta$Proteins))
length(unique(raw_dta$Protein)) 
length(unique(raw_dta$Protein_names))
# Summary: there are multiple variations (modifications) per protein

# Count number of unique genes?
length(unique(raw_dta$Gene_names)) # 1) 2447

## Preprocessing data
# Proteins must get unique names. Additionally, 
# some proteins do not have an annotated gene name and for those we 
# assign a unique ID based on protein name.  

# Investigate if duplicated gene names exist?
raw_dta$Gene_names %>% duplicated() %>% any() # TRUE

# Make a table of duplicated gene names
raw_dta %>% group_by(Gene_names) %>% dplyr::summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)

# Counting number of proteins without an annotated gene name
sum(is.na(raw_dta$Gene_names)) # 1) 61 proteins with missing gene names
missingGenes <- raw_dta[which(is.na(raw_dta$Gene_names)),]
dim(missingGenes) # 61 29

# Test for duplicated gene names?
raw_dta$Gene_names %>% duplicated() %>% any() # TRUE

# Finding symbols for proteins with missing gene and retrieving them
misymbol <- AnnotationDbi::select(org.Mm.eg.db, 
  keys=unique(missingGenes$Protein),columns="SYMBOL", keytype="UNIPROT")
# Collapse genes for different uniprot entries
misymbol <- distinct(misymbol %>% 
  group_by(UNIPROT) %>%
  mutate(genes = paste0(SYMBOL, collapse = ";")),UNIPROT, .keep_all = TRUE)
misymbol <- misymbol[!is.na(misymbol$SYMBOL),]
data_unique_no_mis <- raw_dta
# Set gene names where missing 
data_unique_no_mis$Gene_names[which(data_unique_no_mis$Protein %in% 
    misymbol$UNIPROT)] <- purrr::discard(misymbol$SYMBOL[
      match(data_unique_no_mis$Protein,misymbol$UNIPROT)],is.na) 
# Find uniprot ID for genes not found in org.Mm.eg.db
misymbol <-unique(data_unique_no_mis$Protein
                  [is.na(data_unique_no_mis$Gene_names)])
esearch <- NULL
for (i in 1:length(misymbol)){
  esearch <- rbind(esearch, entrez_search(db="protein", term=misymbol[i])$ids)
}
# Add missing gene for uniprot ids without gene symbol
misymbol <- data.frame(UNIPROT=misymbol, SYMBOL=AnnotationDbi::
                         select(org.Mm.eg.db, keys=esearch[,1], columns=c('SYMBOL'), 
                                keytype="ENTREZID")$SYMBOL)
misymbol <- data.frame(UNIPROT=misymbol, SYMBOL=AnnotationDbi::
            select(org.Mm.eg.db, keys=esearch[,1], columns=c('SYMBOL'), 
              keytype="ENTREZID")$SYMBOL)
data_unique_no_mis$Gene_names[which(data_unique_no_mis$Protein %in% 
    misymbol$UNIPROT)] <- purrr::discard(misymbol$SYMBOL[
      match(data_unique_no_mis$Protein,misymbol$UNIPROT)],is.na) 

# Making unique names using the annotation in the "Gene_names" column as primary 
# names 
data_unique <- make_unique(data_unique_no_mis,
                           "Gene_names", "Protein", delim = ";")
dim(data_unique) # 7754   23
head(data_unique[,c("Gene_names","Protein","name")])

# Testing for any duplicated names?
data_unique$name %>% duplicated() %>% any() # FALSE

colnames(data_unique)
# Removing rows/proteins with all zeros across samples
n <- data_unique[rowSums(data_unique[
  ,grep("raw",colnames(data_unique))] == 0) < 6, ]
dim(n) # 1) 7754   31

# Generating a SummarizedExperiment object using the experimental design
experimental_design <- read.csv("ed.txt", header=T,sep = "\t")

# Generating a SummarizedExperiment object using the experimental design
intensity_columns <- grep("raw", colnames(n)) # get intensity column numbers
data_se <- make_se(n, intensity_columns, experimental_design)

# Counting number of identified phosphopeptides per sample
n_prot <- plot_numbers(data_se, plot=FALSE)
n_prot <- n_prot[c(1,2,5,6,3,4),]
n_prot$condition <- factor(n_prot$condition,levels = c("BSA", "ISM1", "INS"))
ggplot(n_prot, aes(x = condition, y = proteins_in_sample,fill=condition)) + 
  geom_bar(stat = "summary", fun=mean) + 
  labs(y="Number of detected phosphosites") +
  scale_fill_manual(values = brewer.pal(9, "Set1")[c(1,3,2)]) +
  geom_dotplot(binaxis='y', stackdir='center',binwidth = 200) + theme_bw() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=10))
ggsave(paste0("plots/",paste0("number_phosphosites",".jpeg")),
       units="cm", width=4, height=8, dpi=600)
ggsave(file=paste0("plots/",paste0("number_phosphosites",".svg")),
       width=4, height=8)

## Make venn diagrams of identified phosphopeptides between treatments
# Make diagram with no missing value for a treatment
peptide_found <- data.frame(BSA=
                              apply(data_se@assays@data@listData[[1]][,1:2], 
                                     1, function(x) sum(is.na(x))==0),
                            INS=
                              apply(data_se@assays@data@listData[[1]][,3:4], 
                                    1, function(x) sum(is.na(x))==0),
                            ISM1=
                              apply(data_se@assays@data@listData[[1]][,5:6], 
                                    1, function(x) sum(is.na(x))==0))
ggvenn(peptide_found) 
ggsave(file=paste0("plots/",paste0("number_shared_proteins_venn_no_na",".svg")),
       width=4, height=4)
# Make diagram with one allowed missing value for a treatment
peptide_found <- data.frame(BSA=
                              apply(data_se@assays@data@listData[[1]][,1:2], 
                                    1, function(x) sum(is.na(x))<2),
                            INS=
                              apply(data_se@assays@data@listData[[1]][,3:4], 
                                    1, function(x) sum(is.na(x))<2),
                            ISM1=
                              apply(data_se@assays@data@listData[[1]][,5:6], 
                                    1, function(x) sum(is.na(x))<2))
ggvenn(peptide_found) 
ggsave(file=paste0("plots/",paste0("number_shared_proteins_venn_1_na_ok",".svg")),
       width=4, height=4)

# Checking that variation is similar across proteins
meanSdPlot(data_se) + theme_bw()

# Normalizing the data. N.B. Normalization of data introduces missing values 
# instead of zeros
data_norm <- normalize_vsn(data_se)  
meanSdPlot(data_norm) + theme_bw() 

# Checking data distribution after normalization with VSN
qqPlot(assay(data_norm)[1,])
qqPlot(assay(data_norm)[,])
qqPlot(assay(data_norm)[1:1000,])
qqPlot(assay(data_norm)[,2],id=FALSE,ylab="normalized value", main="normalized data") +
  theme(panel.grid = element_blank())
write.table(assay(data_norm)[,2],file=paste0("review/tables/",
                                paste0("qqplot",".tsv")),sep="\t")
shapiro.test(sample_n(as.data.frame(assay(data_norm)[,2][!is.na(assay(data_norm)[,2])]), 500)[,1])
# Data are trending towards a normal distribution except for extremes of qqplot 
# However, using Shapiro-wilks test, data are not normally distributed.

# Checking data distribution after normalization with log2
qqPlot(log2(assay(data_norm)[1,]))
qqPlot(log2(assay(data_norm)[,]))
qqPlot(log2(assay(data_norm))[1:1000,])
qqPlot(log2(assay(data_norm))[,2],id=FALSE,ylab="normalized value", main="normalized and log2 transformed data") + 
  theme(panel.grid = element_blank())
write.table(log2(assay(data_norm))[,2],file=paste0("review/tables/",
                                             paste0("qqplot_log",".tsv")),sep="\t")
# qqplot is improved using log2 for normalization. Data are quazi-linear. 
shapiro.test(sample_n(as.data.frame(log2(assay(data_se))[,2][!is.na(assay(data_norm)[,2])]), 500)[,1])
# After log transformation data are not consistingly deviating form normality.

# 
data_norm_log <- data_norm  
assay(data_norm_log) <- log2(assay(data_norm_log))
meanSdPlot(data_norm_log) + theme_bw() 

# Visualizing normalization by boxplots for all samples before and after 
# normalization
plot_normalization(data_se, data_norm)
p <- plot_normalization(data_se, data_norm_log) + facet_wrap(~var, scales = "free")
p
saveRDS(data_se,file="review/tables/data_se.rds")
saveRDS(data_se,file="review/tables/data_norm_log.rds")
# Normalization and log2 transformation leads to equal variance across samples

# Plotting a heatmap of proteins with missing values
plot_missval(data_se)

# Plotting intensity distributions and cumulative fraction of proteins 
plot_detect(data_se) 
# Since data are missing not at random, missing data are imputed 
# using a left-censored imputation method

saveRDS(data_norm,file="data_norm.rds")
saveRDS(data_norm_log,file="data_norm_log.rds")


