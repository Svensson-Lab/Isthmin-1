### Data analysis - 2, phosphoproteomics, F442A cells
# Script used for analyzing phophoproteome composition and calculate 
# significant differences
# Niels Banhos Danneskiold-Samsoee, May 16 2022

.libPaths(c("C:/R/R-4.1.1/lib"))
.libPaths()
# clear workspace
rm(list=ls())

library(Matrix)
library(readxl)
library(DEP)
library(dplyr)
library(SummarizedExperiment)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(purrr)
library(ggplot2)
library(RColorBrewer)

# set working directory
setwd("C:/Postdoc Stanford/Projects/Isthmin phoshoproteomics")

# Loading raw data and experimental design
data_norm <- readRDS("data_norm.rds")
experimental_design <- read.csv("ed.txt", header=T,sep = "\t")
pp_norm <- data_norm

# Calculating PCA scores using all proteins
d <- assay(pp_norm) %>% data.frame()
d[is.na(d)] <- 0
phosprot_pca <- prcomp(as.matrix(t(d)),center = TRUE, scale. = TRUE)

#drawing PCA scores using all proteins as input
pc1<-phosprot_pca$x[,1]
pc2<-phosprot_pca$x[,2]
scores<-data.frame(pc1,pc2)
scores$treatment <- experimental_design$condition
scores$replicate <- experimental_design$replicate
ggplot(data=scores, aes()) +
  geom_point(aes(pc1, pc2,color = factor(treatment),fill = factor(treatment)),
        size=4, alpha=0.75) +
  scale_shape_manual(values=c(22,21,24), name="replicate") +
  scale_fill_manual(values=brewer.pal(9, "Set1"), name="treatment") +
  scale_color_manual(values=brewer.pal(9, "Set1"), name="treatment") +
  xlab(paste("PC1  (",toString(round(
    summary(phosprot_pca)$importance[2]*100,2)),"%)", collapse = ",")) +
  ylab(paste("PC2  (",toString(round(
    summary(phosprot_pca)$importance[5]*100,2)),"%)", collapse = ",")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size=16,family="sans"),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size=16,family="sans"),
        legend.text = element_text(size=16,family="sans"))
ggsave(paste0("plots/",paste0("PCA",".jpeg")),
       units="cm", width=8, height=6, dpi=600)
ggsave(file=paste0("plots/",paste0("PCA",".svg")),
       width=8, height=6)

# Investigating effect of imputation and imputation algorithm on 
# significant phosphopeptides
data_imp <- impute(pp_norm, fun = "MinProb", q = 0.01)
#data_imp <- impute(pp_norm, fun = "QRILC")
#data_imp <- impute(pp_norm, fun = "man")
data_se@assays@data@listData[[1]][grep("Insr",
            rownames(data_se@assays@data@listData[[1]])),]
data_imp@assays@data@listData[[1]][grep("Insr",
            rownames(data_imp@assays@data@listData[[1]])),] 
# Summary: all imputation algorithms (QRILC, MinProb or man) 
# sets missing values to low intensity despite the lack of
# signal for both replicates. It is likely not reflecting the true value
# as the insulin receptor should not be phosphorylated by BSA.

# Plotting intensity distributions and cumulative fraction of 
# proteins with and without missing values
plot_detect(pp_norm)
# Summary: proteins with missing values have lower intensities

# Plotting intensity distributions before and after imputation
plot_imputation(pp_norm, data_imp)

# Testing whether insulin significantly phosporylates the insulin receptor 
# using imputation
data_imp <- impute(pp_norm, fun = "MinProb", q = 0.01) 
data_diff_BSA_vs_INS <- test_diff(data_imp, type = "manual",
                        test = "BSA_vs_INS")
# Denoting significant proteins
dep <- add_rejections(data_diff_BSA_vs_INS, alpha = 0.05, lfc = 0)

table(dep@elementMetadata@listData$significant) # ~244 significant polypeptides
as.data.frame(dep@elementMetadata@listData)[grep("Insr",
        dep@elementMetadata@listData[[1]]),] 
# As sample size is low, non-deterministic imputation results in
# insulin receptor not always significantly phosphorylated by insulin

# Using permutations to mitigate effects of non-deterministic variation in 
# results due to imputation. Testing number of permutations necessary to avoid 
# differences in significant genes.
p.vals <- data.frame()
padj.vals <- data.frame()
differences <- data.frame()
sig_permute_minprob <- NULL
pp_name <- data.frame(pp_name=data_norm@elementMetadata@listData$name)
permutations <- 200
# Loop that permutes imputation and significance testing. Compares if more
# permutations lead to further smoothing of results.
for (i in 1:permutations){
  data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
  data_diff_BSA_vs_ISM1 <- suppressMessages(test_diff(
    data_imp, type = "manual", test = "ISM1_vs_BSA"))
  data_diff_BSA_vs_ISM1 <- add_rejections(data_diff_BSA_vs_ISM1, 
                                          alpha = 0.05, lfc = 0)
  padj.val <- data_diff_BSA_vs_ISM1@elementMetadata@listData$ISM1_vs_BSA_p.adj
  if (i==1){
    padj.vals <- cbind(pp_name,padj.val)
  }
  if (i==2){
    p.vals.old <- data.frame(pp_name, padj.vals=padj.vals[,2])
    p.vals.old <- pp_name$pp_name[p.vals.old$padj.vals <0.05]
    padj.vals <- cbind(padj.vals,padj.val)
    p.vals.new <- data.frame(pp_name,
                             padj.vals=apply(padj.vals[,2:(3)], 1, 
                                             FUN = median,na.rm=TRUE))
    p.vals.new <- pp_name$pp_name[p.vals.new$padj.vals <0.05]
    differences <- rbind(differences,length(setdiff(p.vals.old,p.vals.new)))
  }
  if (i>2){
    p.vals.old <- data.frame(pp_name,
                             padj.vals=apply(padj.vals[,2:i], 1, 
                                             FUN = median,na.rm=TRUE))
    p.vals.old <- pp_name$pp_name[p.vals.old$padj.vals <0.05]
    padj.vals <- cbind(padj.vals,padj.val)
    # Finding out whether adding one test changes the conclusion on a peptide 
    p.vals.new <- data.frame(pp_name,
                             padj.vals=apply(padj.vals[,2:(i+1)], 1, 
                                             FUN = median,na.rm=TRUE))
    p.vals.new <- pp_name$pp_name[p.vals.new$padj.vals <0.05]
    differences <- rbind(differences,length(setdiff(p.vals.old,p.vals.new)))
  }  
  if (i %% 10 ==0) print(i)
  if (i==permutations){
    
  }
}
head(differences)
colnames(differences) <- c("diff")
differences$sum <- cumsum(differences$diff)
differences$permutation <- row.names(differences)
head(differences)
plot(differences$permutation,differences$sum, xlab="Permutation",
     ylab="Cumulative differences in number of significant peptides by adding an extra permutation") 
# Summary: 200 permutations seem to be sufficient to smooth variation in
# significance due to imputation

# Finding median p-values and median fold changes testing ISM1 against BSA, 
# as well as median imputed values across permutations using MinProp
p.vals <- data.frame()
padj.vals <- data.frame()
fc.vals <- data.frame()
p.vals.ISMvsBSA <- data.frame()
padj.vals.ISMvsBSA <- data.frame()
fc.vals.ISMvsBSA <- data.frame()
pp_name <- data.frame(pp_name=pp_norm@elementMetadata@listData$name)
bsa1 <- data.frame()
bsa2 <- data.frame()
ins1 <- data.frame()
ins2 <- data.frame()
ism1 <- data.frame()
ism2 <- data.frame()
imputed.vals <- data.frame()
permutations <- 200
for (i in 1:permutations){
  data_imp <- impute(pp_norm, fun = "MinProb", q = 0.01)
  data_diff_BSA_vs_ISM1 <- suppressMessages(test_diff(
    data_imp, type = "manual", test = "ISM1_vs_BSA"))
  data_diff_BSA_vs_ISM1 <- add_rejections(data_diff_BSA_vs_ISM1, 
                                          alpha = 0.05, lfc = 0)
  p.val <- data_diff_BSA_vs_ISM1@elementMetadata@listData$ISM1_vs_BSA_p.val
  padj.val <- data_diff_BSA_vs_ISM1@elementMetadata@listData$ISM1_vs_BSA_p.adj
  fc <- data_diff_BSA_vs_ISM1@elementMetadata@listData$ISM1_vs_BSA_diff
  if (i==1){
    p.vals <- cbind(pp_name,p.val)
    padj.vals <- cbind(pp_name,padj.val)
    fc.vals <- cbind(pp_name,fc)
    bsa1 <- as.data.frame(assay(data_diff_BSA_vs_ISM1))$BSA_1
    bsa2 <- as.data.frame(assay(data_diff_BSA_vs_ISM1))$BSA_2
    ins1 <- as.data.frame(assay(data_diff_BSA_vs_ISM1))$INS_1
    ins2 <- as.data.frame(assay(data_diff_BSA_vs_ISM1))$INS_2
    ism1 <- as.data.frame(assay(data_diff_BSA_vs_ISM1))$ISM1_1
    ism2 <- as.data.frame(assay(data_diff_BSA_vs_ISM1))$ISM1_2
  }
  if (i>1){
    p.vals <- cbind(p.vals,p.val)
    padj.vals <- cbind(padj.vals,padj.val)
    fc.vals <- cbind(fc.vals,fc)    
    
    bsa1 <- cbind(bsa1,as.data.frame(assay(data_diff_BSA_vs_ISM1))$BSA_1)
    bsa2 <- cbind(bsa2,as.data.frame(assay(data_diff_BSA_vs_ISM1))$BSA_2)
    ins1 <- cbind(ins1,as.data.frame(assay(data_diff_BSA_vs_ISM1))$INS_1)
    ins2 <- cbind(ins2,as.data.frame(assay(data_diff_BSA_vs_ISM1))$INS_2)
    ism1 <- cbind(ism1,as.data.frame(assay(data_diff_BSA_vs_ISM1))$ISM1_1)
    ism2 <- cbind(ism2,as.data.frame(assay(data_diff_BSA_vs_ISM1))$ISM1_2)
  }  
  if (i %% 10 ==0) print(i)
  if (i==permutations){
    p.vals.ISMvsBSA <- data.frame(pp_name,
      p.vals.median=apply(p.vals[,2:permutations], 1, FUN = median,na.rm=TRUE))
    padj.vals.ISMvsBSA <- data.frame(pp_name,
      padj.vals.median=apply(padj.vals[,2:permutations], 1, 
                             FUN = median,na.rm=TRUE))
    fc.vals.ISMvsBSA <- data.frame(pp_name,
      fc.vals.median=apply(fc.vals[,2:permutations], 1, 
                           FUN = median,na.rm=TRUE))
    imputed.vals <- data.frame(pp_name,
      BSA_1=apply(bsa1[,1:permutations], 1, FUN = median,na.rm=TRUE),
      BSA_2=apply(bsa2[,1:permutations], 1, FUN = median,na.rm=TRUE),
      INS_1=apply(ins1[,1:permutations], 1, FUN = median,na.rm=TRUE),
      INS_2=apply(ins2[,1:permutations], 1, FUN = median,na.rm=TRUE),
      ISM1_1=apply(ism1[,1:permutations], 1, FUN = median,na.rm=TRUE),
      ISM1_2=apply(ism2[,1:permutations], 1, FUN = median,na.rm=TRUE))
    
  }
}
rm(bsa1,bsa2,ins1,ins2,ism1,ism2, p.vals, padj.vals, fc.vals)
rownames(imputed.vals) <- imputed.vals$pp_name
imputed.vals <- imputed.vals[,-1]
head(imputed.vals)
head(p.vals.ISMvsBSA)
table(p.vals.ISMvsBSA < 0.05)
table(padj.vals.ISMvsBSA < 0.05)
head(fc.vals.ISMvsBSA)

# Adding p-values and fold change to normalized data
pp_norm@elementMetadata@listData$p.vals.ISMvsBSA <- 
  p.vals.ISMvsBSA$p.vals.median
pp_norm@elementMetadata@listData$padj.vals.ISMvsBSA <- 
  padj.vals.ISMvsBSA$padj.vals.median
pp_norm@elementMetadata@listData$fc.vals.ISMvsBSA <- 
  fc.vals.ISMvsBSA$fc.vals.median

# Finding median p-values and median fold changes testing INS against BSA
p.vals <- data.frame()
padj.vals <- data.frame()
fc.vals <- data.frame()
p.vals.INSvsBSA <- data.frame()
padj.vals.INSvsBSA <- data.frame()
fc.vals.INSvsBSA <- data.frame()
pp_name <- data.frame(pp_name=pp_norm@elementMetadata@listData$name)
permutations <- 200
for (i in 1:permutations){
  data_imp <- impute(pp_norm, fun = "MinProb", q = 0.01)
  data_diff_BSA_vs_INS <- suppressMessages(test_diff(data_imp, 
    type = "manual", test = "INS_vs_BSA"))
  data_diff_BSA_vs_INS <- add_rejections(data_diff_BSA_vs_INS, 
    alpha = 0.05, lfc = 0)
  p.val <- data_diff_BSA_vs_INS@elementMetadata@listData$INS_vs_BSA_p.val
  padj.val <- data_diff_BSA_vs_INS@elementMetadata@listData$INS_vs_BSA_p.adj
  fc <- data_diff_BSA_vs_INS@elementMetadata@listData$INS_vs_BSA_diff
  if (i==1){
    p.vals <- cbind(pp_name,p.val)
    padj.vals <- cbind(pp_name,padj.val)
    fc.vals <- cbind(pp_name,fc)
  }
  if (i>1){
    p.vals <- cbind(p.vals,p.val)
    padj.vals <- cbind(padj.vals,padj.val)
    fc.vals <- cbind(fc.vals,fc)    
  }  
  if (i %% 10 ==0) print(i)
  if (i==permutations){
    p.vals.INSvsBSA <- data.frame(pp_name,p.vals.median=apply(
      p.vals[,2:permutations], 1, FUN = median,na.rm=TRUE))
    padj.vals.INSvsBSA <- data.frame(
      pp_name,padj.vals.median=apply(padj.vals[,2:permutations], 1, 
                                     FUN = median,na.rm=TRUE))
    fc.vals.INSvsBSA <- data.frame(pp_name,fc.vals.median=apply(
      fc.vals[,2:permutations], 1, FUN = median,na.rm=TRUE))
  }
}
rm(p.vals, padj.vals, fc.vals)
head(p.vals.INSvsBSA)
table(p.vals.INSvsBSA < 0.05)
table(padj.vals.INSvsBSA < 0.05)
head(fc.vals.INSvsBSA)

# Adding p-values and fold change to normalized data
pp_norm@elementMetadata@listData$p.vals.INSvsBSA <- 
  p.vals.INSvsBSA$p.vals.median
pp_norm@elementMetadata@listData$padj.vals.INSvsBSA <- 
  padj.vals.INSvsBSA$padj.vals.median
pp_norm@elementMetadata@listData$fc.vals.INSvsBSA <- 
  fc.vals.INSvsBSA$fc.vals.median

# Finding median p-values and median fold changes testing ISM1 against INS
p.vals <- data.frame()
padj.vals <- data.frame()
fc.vals <- data.frame()
p.vals.ISMvsINS <- data.frame()
padj.vals.ISMvsINS <- data.frame()
fc.vals.ISMvsINS <- data.frame()
pp_name <- data.frame(pp_name=pp_norm@elementMetadata@listData$name)
permutations <- 200
for (i in 1:permutations){
  data_imp <- impute(pp_norm, fun = "MinProb", q = 0.01)
  data_diff_INS_vs_ISM1 <- suppressMessages(test_diff(data_imp, 
    type = "manual", test = "ISM1_vs_INS"))
  data_diff_INS_vs_ISM1 <- add_rejections(data_diff_INS_vs_ISM1, 
    alpha = 0.05, lfc = 0)
  p.val <- data_diff_INS_vs_ISM1@elementMetadata@listData$ISM1_vs_INS_p.val
  padj.val <- data_diff_INS_vs_ISM1@elementMetadata@listData$ISM1_vs_INS_p.val
  fc <- data_diff_INS_vs_ISM1@elementMetadata@listData$ISM1_vs_INS_diff
  
  if (i==1){
    p.vals <- cbind(pp_name,p.val)
    padj.vals <- cbind(pp_name,padj.val)
    fc.vals <- cbind(pp_name,fc)
  }
  if (i>1){
    p.vals <- cbind(p.vals,p.val)
    padj.vals <- cbind(padj.vals,padj.val)
    fc.vals <- cbind(fc.vals,fc)    
  }  
  if (i %% 10 ==0) print(i)
  if (i==permutations){
    p.vals.ISMvsINS <- data.frame(pp_name,p.vals.median=apply(
      p.vals[,2:permutations], 1, FUN = median,na.rm=TRUE))
    padj.vals.ISMvsINS <- data.frame(pp_name,padj.vals.median=apply(
      padj.vals[,2:permutations], 1, FUN = median,na.rm=TRUE))
    fc.vals.ISMvsINS <- data.frame(pp_name,fc.vals.median=apply(
      fc.vals[,2:permutations], 1, FUN = median,na.rm=TRUE))
  }
}
head(p.vals.ISMvsINS)
table(p.vals.ISMvsINS < 0.05)
table(padj.vals.ISMvsINS < 0.05)
head(fc.vals.ISMvsINS)

# Adding p-values and fold change to normalized data
pp_norm@elementMetadata@listData$p.vals.ISMvsINS <- 
  p.vals.ISMvsINS$p.vals.median
pp_norm@elementMetadata@listData$padj.vals.ISMvsINS <- 
  padj.vals.ISMvsINS$padj.vals.median
pp_norm@elementMetadata@listData$fc.vals.ISMvsINS <- 
  fc.vals.ISMvsINS$fc.vals.median

# Finding median p-values and median fold changes testing all conditions
# against each other
all_sigs <- data.frame()
BSA_vs_ISM1_diffs <- data.frame()
BSA_vs_INS_diffs <- data.frame()
INS_vs_ISM1_diffs <- data.frame()
pp_name <- data.frame(pp_name=pp_norm@elementMetadata@listData$name)
permutations <- 200
for (i in 1:permutations){
  data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
  data_diff <- suppressMessages(test_diff(data_imp, type = "all"))
  data_diff <- add_rejections(data_diff)
  all_sig <- data_diff@elementMetadata@listData$significant
  BSA_vs_ISM1_diff <- data_diff@elementMetadata@listData$BSA_vs_ISM1_diff
  BSA_vs_INS_diff <- data_diff@elementMetadata@listData$BSA_vs_INS_diff
  INS_vs_ISM1_diff <- data_diff@elementMetadata@listData$INS_vs_ISM1_diff
  
  if (i==1){
    all_sigs <- cbind(pp_name,all_sig)
    BSA_vs_ISM1_diffs <- cbind(pp_name,BSA_vs_ISM1_diff)
    BSA_vs_INS_diffs <- cbind(pp_name,BSA_vs_INS_diff)
    INS_vs_ISM1_diffs <- cbind(pp_name,INS_vs_ISM1_diff)
  }
  if (i>1){
    all_sigs <- cbind(all_sigs,all_sig)
    BSA_vs_ISM1_diffs <- cbind(BSA_vs_ISM1_diffs,BSA_vs_ISM1_diff)
    BSA_vs_INS_diffs <- cbind(BSA_vs_INS_diffs,BSA_vs_INS_diff)
    INS_vs_ISM1_diffs <- cbind(INS_vs_ISM1_diffs,INS_vs_ISM1_diff)
    
  }  
  if (i %% 10 ==0) print(i)
  if (i==permutations){
    all_sigs <- data.frame(pp_name,
      sig.all.median=apply(all_sigs[,2:permutations], 1, 
               FUN = median,na.rm=TRUE))
      BSA_vs_ISM1_diffs <- data.frame(pp_name,
              BSA_vs_ISM1_diff=apply(BSA_vs_ISM1_diffs[,2:permutations], 1, 
              FUN = median,na.rm=TRUE))
      BSA_vs_INS_diffs <- data.frame(pp_name,
              BSA_vs_INS_diff=apply(BSA_vs_INS_diffs[,2:permutations], 1, 
              FUN = median,na.rm=TRUE))
      INS_vs_ISM1_diffs <- data.frame(pp_name,
              INS_vs_ISM1_diff=apply(INS_vs_ISM1_diffs[,2:permutations], 1, 
              FUN = median,na.rm=TRUE))
  }
}
head(all_sigs)
table(all_sigs$sig.all.median)
head(BSA_vs_ISM1_diffs[all_sigs$sig.all.median,])
head(BSA_vs_INS_diffs[all_sigs$sig.all.median,])
head(INS_vs_ISM1_diffs[all_sigs$sig.all.median,])

# Adding p-values and fold change to normalized data
pp_norm@elementMetadata@listData$significant <- all_sigs$sig.all.median
pp_norm@elementMetadata@listData$BSA_vs_ISM1_diff <- 
  BSA_vs_ISM1_diffs$BSA_vs_ISM1_diff
pp_norm@elementMetadata@listData$BSA_vs_INS_diff <- 
  BSA_vs_INS_diffs$BSA_vs_INS_diff
pp_norm@elementMetadata@listData$INS_vs_ISM1_diff <- 
  INS_vs_ISM1_diffs$INS_vs_ISM1_diff

# Saving data including permuted statistical differences 
saveRDS(pp_norm,file="pp_norm.rds")
