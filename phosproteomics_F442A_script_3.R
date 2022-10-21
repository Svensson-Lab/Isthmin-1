### Data analysis - 3, phosphoproteomics, F442A cells
# Script used for heatmap generation and pathway analysis
# Niels Banhos Danneskiold-Samsoee, Sep 16 2022

# clear workspace
rm(list=ls())

library(Matrix)
library(readxl)
library(DEP)
library(dplyr)
library(SummarizedExperiment)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(stats)
library(dendsort)
library(wesanderson)
library(grid)
library(topGO)
library(rrvgo)
library(MASS)
library(ComplexHeatmap)
library(compositions)
library(dendsort)
library(GO.db)
library(dittoSeq)
library(stringr)
library(limma)

## set working directory
setwd("~/Isthmin phoshoproteomics")

# Loading raw data and experimental design
pp_norm <- readRDS(file="pp_norm_log_filter.rds")
experimental_design <- read.csv("ed.txt", header=T,sep = "\t")

### Drawing heatmap of phosphopeptides with lowest adjusted p-values
# Setup clustering function
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

# finding unique phosphopeptides with lowest 100 p-values for each comparison  
pos_low_p_val <- unique(c(sort(pp_norm@elementMetadata$padj.vals.INSvsBSA,
                               index.return=TRUE)$ix[1:100],
                          sort(pp_norm@elementMetadata$padj.vals.ISMvsBSA,
                               index.return=TRUE)$ix[1:100],
                          sort(pp_norm@elementMetadata$padj.vals.ISMvsINS,
                               index.return=TRUE)$ix[1:100]))

# Calculating center log ratio of normalized intensities
clr_int <- clr(as.data.frame(pp_norm@elementMetadata[,
                          grep("Intensity",colnames(pp_norm@elementMetadata))]))
colnames(clr_int) <- unname(as.data.frame(strsplit(colnames(clr_int),"_"))[2,])
# Cluster phosphopeptides
mat_cluster_rows <- sort_hclust(hclust(dist(clr_int[pos_low_p_val,],
                          method="minkowski"),method="ward.D")) 
# Make heatmap with imputed values instead of non-detectable
ha <- HeatmapAnnotation(condition = c("BSA","BSA","INS","INS","ISM1","ISM1"),
                        col = list(condition= c(BSA=brewer.pal(9, "Set1")[1],
                                                INS=brewer.pal(9, "Set1")[2],
                                                ISM1=brewer.pal(9, "Set1")[3])),
                        annotation_legend_param = list(condition = list(
                          title = "condition",
                          at = c("BSA", "INS","ISM1"),
                          labels = c("BSA", "INS","ISM1"))))
h<- Heatmap(clr_int[pos_low_p_val,]*-1,
        na_col = "grey",clustering_distance_columns ="canberra",
        cluster_rows=mat_cluster_rows,col = circlize::colorRamp2(
          seq(-1, 1, (1/5)),
          rev(RColorBrewer::brewer.pal(11, "RdBu"))),column_km = 3,
        top_annotation = ha,
        row_dend_width=unit(10, "cm"),column_dend_height=unit(2, "cm"),
        show_column_names = FALSE, 
        heatmap_legend_param=list(title = "Log2 ratio"))
h
grid.ls(grid.force()) # "col_annotation" looks like it's the one to edit
grid.gedit("GRID.rect", gp = gpar(col="grey70"))
dev.copy2pdf(file=paste0("review/plots/heatmap_top100_grid_no_na_filtered",paste0(".pdf")), 
             width=6.5, height=12)
dev.off()
h
dev.copy2pdf(file=paste0("review/plots/heatmap_top100_nogrid_no_na",paste0(".pdf")), 
             width=6.5, height=12)
dev.off()

### GO analysis

## Making graphs of GO enrichment for ISM1 vs BSA
# Making vector with adjusted p-values
gl <- pp_norm@elementMetadata$padj.vals.ISMvsBSA
# Investigating cases where genes cannot be uniquely assigned
pp_norm@elementMetadata$Gene_names[grep(";",pp_norm@elementMetadata$Gene_names)]
# Summary: most gene names with multiple gene assignments belong to the same
# family. 
# Restricting results to first gene
names(gl) <- gsub(";\\S+","",pp_norm@elementMetadata$Gene_names)
head(gl)
#gl <- gl[pp_norm@elementMetadata$fc.vals.ISMvsBSA >0] # Testing results using only phosphopeptides with higher phosphorylation
table(gl <0.05)
selection <- function(allScore){ return(allScore < 0.05)} 
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, 
                          mapping="org.Mm.eg.db", ID="symbol")
# Fetch go enrichment data
GOdata <- new("topGOdata", 
                    ontology = "BP", 
                    allGenes = gl, 
                    geneSel = selection, 
                    nodeSize = 10, annot = annFUN.GO2genes, 
                    GO2genes=allGO2genes)
# Testing more stringent enrichment
# Save pathway analysis
saveRDS(GOdata,"GO_data_log_filtered.rds")

Fisher <- runTest(GOdata, statistic = "fisher")
# Fisher test results in 97 GO terms with p<0.01
# generate the result of the fisher test
sig.tab <- GenTable(GOdata, Fis = Fisher, topNodes = 800)
sig.tab <- sig.tab[is.na(as.numeric(sig.tab$Fis)) | 
                     as.numeric(sig.tab$Fis) <0.05,]
# Get GO terms
goterms <- unlist(Term(GOTERM))
sig.tab$Term <- unname(goterms[match(sig.tab$GO.ID,names(goterms))])
KS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
# Stringent KS test results in 25 GO terms with p<0.01
KS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
# Less stringent KS test results in 1586 GO terms with p<0.01
# generate the result of less stringent Kolmogorov-Smirnov test
sig.tab <- GenTable(GOdata, classicFisher = Fisher,
                    classicKS = KS, elimKS = KS.elim,
                    orderBy = "classicKS", 
                    ranksOf = "classicFisher", topNodes = 2500)
sig.tab$KSnumeric <- as.numeric(gsub("<","",sig.tab$classicKS))
sig.tab <- sig.tab[sig.tab$KSnumeric <0.01,]
sig.tab$Term <- unname(goterms[match(sig.tab$GO.ID,names(goterms))])
# Investigate pathways for common linkages to Isthmin-1 and metabolic organs
sig.tab[grep("neuro",sig.tab$Term),] 
dim(sig.tab[grep("neuro",sig.tab$Term),]) #in 28 terms
sig.tab[grep("immune",sig.tab$Term),]  
dim(sig.tab[grep("immune",sig.tab$Term),]) #in 11 terms 
sig.tab[grep("adipo",sig.tab$Term),] 
dim(sig.tab[grep("adipo",sig.tab$Term),]) #in 0 terms
sig.tab[grep("fat cell",sig.tab$Term),] 
dim(sig.tab[grep("fat cell",sig.tab$Term),]) #in 3 terms
sig.tab[grep("liver|hepat",sig.tab$Term),] 
dim(sig.tab[grep("liver|hepat",sig.tab$Term),]) #liver in 1 term, hepat in 1 term
sig.tab[grep("pancrea",sig.tab$Term),] 
dim(sig.tab[grep("pancrea",sig.tab$Term),]) #in 0 terms
sig.tab[grep("muscle",sig.tab$Term),] 
dim(sig.tab[grep("muscle",sig.tab$Term),]) #in 25 terms 
sig.tab[grep("lung|pulmo",sig.tab$Term),] 
dim(sig.tab[grep("lung|pulmo",sig.tab$Term),]) #lung in 1 terms, pulmo in 0 terms
sig.tab[grep("kidney",sig.tab$Term),] 
dim(sig.tab[grep("kidney",sig.tab$Term),]) #in 1 term
sig.tab[grep("spleen",sig.tab$Term),] 
dim(sig.tab[grep("spleen",sig.tab$Term),]) #in 0 terms
sig.tab[grep("intest",sig.tab$Term),] 
dim(sig.tab[grep("intest",sig.tab$Term),]) #in 0 terms
sig.tab[grep("adren",sig.tab$Term),] 
dim(sig.tab[grep("adren",sig.tab$Term),]) #in 0 terms
sig.tab[grep("skin",sig.tab$Term),] 
dim(sig.tab[grep("skin",sig.tab$Term),]) #in 1 terms
sig.tab[grep("bladder",sig.tab$Term),] #in 0 terms
sig.tab[grep("sex",sig.tab$Term),] #in 4 terms
sig.tab[grep("thyroid",sig.tab$Term),] #in 0 terms
sig.tab[grep("thymus",sig.tab$Term),] #in 2 terms

head(sig.tab)
tail(sig.tab)
sig.tab$escore <- sig.tab$Significant/sig.tab$Expected

# saving results of pathway analysis
write.table(sig.tab,file=paste0("review/tables/",
  paste0("2022-09-01_GOs_ISMvsBSA_log_filter",".tsv")),sep="\t")
dim(sig.tab)
# Summary: List includes 1586 GO terms with different phosphorylation in 
# ISM1 vs BSA. Similar results obtained if restricting to genes with higher 
# phosphorylation in ISM1 samples. List includes immune functions consistent
# with previous observations, neural development consistent with its 
# importance for fetal brain development and muscle development, but not
# pathways relating to adipose tissue or the pancreas and only two related to
# liver

## Clustering GO terms by semantic similarity
# Calculating similarity matrix
simMatrix <- calculateSimMatrix(sig.tab$GO.ID,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")
saveRDS(simMatrix,"GOsimilaritymatrix_log_filter.rds")

# Grouping terms based on similarity
scores <- setNames(-log10(sig.tab$KSnumeric), 
                   sig.tab$GO.ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Mm.eg.db")
# Reviewing clustering
heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)
# Summary: some clustering based on parent GO term 

# Calculate PCOA based on matrix and cluster dendrogram
sim_pcoa <- cmdscale(as.matrix(as.dist(1-simMatrix)), eig=TRUE, k=2)
plot(hclust(as.dist(1-simMatrix), method="ward.D"))
dend <- hclust(as.dist(1-simMatrix), method="ward.D2") 
branch_cut <- cutree(dend, k = 4)
table(branch_cut)

# prepare PCOA plot
goplot <- cbind(as.data.frame(sim_pcoa$points),
                reducedTerms[match(rownames(sim_pcoa$points), reducedTerms$go), 
                c("term","cluster", "parent", "parentTerm", "size")])
goplot$GOID <- rownames(sim_pcoa$points)
head(sig.tab)
rownames(goplot) <- goplot$GOID
goplot$escore <- sig.tab$escore[which(sig.tab$GO.ID %in% rownames(goplot))]
goplot$classicKS <- sig.tab$classicKS[
  which(sig.tab$GO.ID %in% rownames(goplot))]
goplot$branch <- branch_cut[which(rownames(goplot) %in% names(branch_cut))]
centroid <- aggregate(.~branch_cut, data=
                        data.frame(goplot$V1,goplot$V2,goplot$branch), mean)
centroid <- centroid[,-1]
colnames(centroid) <- c("clust.x","clust.y","branch")
goplot <- left_join(goplot,as.data.frame(centroid), by="branch")
head(branch_cut)
head(goplot)

write.table(goplot,row.names=FALSE,file=paste0("review/tables/",
                                paste0("GOs_ISMvsBSA_clust_filter",".tsv")),sep="\t")
# Clusters dominated by terms:
# 1: mixed
# 2: metabolic process
# 3: development
# 4: localization/transport

goplot$plabel <- "other"
goplot$plabel[grep("immune",goplot$term)] <- "immune system"
goplot$plabel[grep("neuro",goplot$term)] <- "neural"
goplot$plabel[grep("liver|hepatic",goplot$term)] <- "liver"
goplot$plabel[grep("muscle",goplot$term)] <- "muscle"
goplot$plabel[grep("muscle",goplot$term)] <- "smooth muscle"
goplot$plabel[grep("cardiac muscle",goplot$term)] <- "cardiac muscle"
goplot$plabel[grep("striated muscle",goplot$term)] <- "striated muscle"
table(goplot$plabel)

pm_col <- c("other"="azure3", 
            "muscle"=brewer.pal(9,"Reds")[5],
            "smooth muscle"=brewer.pal(9,"Reds")[3], 
            "cardiac muscle"=brewer.pal(9,"Reds")[7],
            "striated muscle"=brewer.pal(9,"Reds")[9],
            "immune system"=dittoColors()[3],
            "neural"=dittoColors()[2],
            "liver"=dittoColors()[7])

# Plot PCOA of similarities in GO terms
p1 <- ggplot(data=goplot, aes(x=V1,y=V2)) + 
  geom_point(data = goplot %>% filter(plabel == "other"), 
             aes(x = V1, y = V2,color = plabel,size=escore), 
             color="azure3", alpha=.3) +
  labs(size = "Enrichment score") +
  geom_point(data = goplot %>% filter(!str_detect(plabel,"muscle|other")), 
             aes(x = V1, y = V2,color = plabel,size=escore), alpha=0.5) +
  geom_point(data = goplot %>% filter(str_detect(plabel,"muscle")), 
             aes(x = V1, y = V2,color = plabel,size=escore), alpha=0.8) +
  geom_point(data = goplot %>% filter(str_detect(plabel,"smooth muscle")), 
             aes(x = V1, y = V2,color = plabel,size=escore), alpha=0.8) +
  geom_point(data = goplot %>% filter(str_detect(plabel,"cardiac muscle")), 
             aes(x = V1, y = V2,color = plabel,size=escore), alpha=0.8) +
  geom_point(data = goplot %>% filter(str_detect(plabel,"striated muscle")), 
             aes(x = V1, y = V2,color = plabel,size=escore), alpha=0.8) +
  scale_color_manual(values=pm_col, name="GO category") + 
  stat_ellipse(aes(x=V1, y=V2,group=as.factor(branch)),type = "norm",level = 0.8) +
  geom_label(aes(x=clust.x+.1, y=clust.y, label=branch)) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank())
p1
dev.copy2pdf(file=paste0("review/plots/pathways_PCOA_plot_log_filter",paste0(".pdf")), width=5.8, height=4)
dev.off()
p1
ggsave(file=paste0("review/plots/",paste0("pathways_PCOA_plot_log_filter",".svg")),
       width=5.8, height=4)
