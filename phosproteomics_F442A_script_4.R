### Data analysis - 4, phosphoproteomics, F442A cells
# Script to make venn diagrams, distribution diagram,
# and barplots of phosphopeptides
# Niels Banhos Danneskiold-Samsoee, May 16 2022

rm(list=ls())

library(UniprotR)
library(GOfuncR)
library(Mus.musculus)
library(ggvenn)
library(GO.db)
library(circlize)
library(RColorBrewer)
#library(GISTools)
library(viridis)
library(dittoSeq)
library(ComplexHeatmap)
library(DescTools)
library(tidyr)
library(compositions)
library(dendsort)
library(SummarizedExperiment)
library(dplyr)

## set working directory
setwd("~/Isthmin phoshoproteomics")

# Loading normalized phosphoproteomics data and experimental design
pp_norm <- readRDS(file="pp_norm_log_filter.rds")
experimental_design <- read.csv("ed.txt", header=T,sep = "\t")

# extract a named vector of all terms
goterms <- Term(GOTERM)
saveRDS(goterms,"2022_09_01_GOterms.rds")

# finding genes associated with the GO term: insulin receptor signaling pathway
is_pathway_genes <- get_anno_genes("GO:0008286", database = 'Mus.musculus')
table(is_pathway_genes$gene)
# finding genes associated with the GO term: insulin receptor signaling pathway
mTOR_pathway_genes <- get_anno_genes(names(goterms[grep("\\bTOR|MTORC1|MTORC2\\b",goterms)]), 
                                     database = 'Mus.musculus')
table(mTOR_pathway_genes$gene)
# finding genes associated with GO terms including: "muscle"
muscle_pathways_genes <- get_anno_genes(names(goterms[grep("muscle",goterms)]), 
                                        database = 'Mus.musculus')
table(muscle_pathways_genes$gene)

#View(goterms[names(goterms) %in% muscle_pathways_genes$go_id])
length(muscle_pathways_genes$gene)
length(unique(muscle_pathways_genes$gene))

# in case of multiple genes, picking first
genes <- gsub(";\\S+","",pp_norm@elementMetadata$Gene_names)

# Annotate if genes associated with insulin receptor Go term is detected in 
# both replicates for each treatment. 
is_peptide_found <- data.frame(BSA=
                                 apply(pp_norm@assays@data@listData[[1]][
                                   genes %in% is_pathway_genes$gene ,1:2], 
                                   1, function(x) sum(is.na(x))<2),
                               INS=
                                 apply(pp_norm@assays@data@listData[[1]][
                                   genes %in% is_pathway_genes$gene ,3:4], 
                                   1, function(x) sum(is.na(x))<2),
                               ISM1=
                                 apply(pp_norm@assays@data@listData[[1]][
                                   genes %in% is_pathway_genes$gene ,5:6], 
                                   1, function(x) sum(is.na(x))<2))
is_peptide_found[grep("Insr",rownames(is_peptide_found)),]
is_peptide_found$phospep <- rownames(is_peptide_found)

# Make venn diagram of shared and unique phosphopeptides 
ggvenn(is_peptide_found)
ggsave(file=paste0("review/plots/",paste0("is_peptide_shared_venn_0NA_log",".svg")),
       width=8, height=6)
is_peptide_found$pathway <- "IR signalling"

# Annotate if genes associated with TOR GO term is detected in 
# both replicates for each treatment. 
mTOR_peptide_found <- data.frame(BSA=
                                   apply(pp_norm@assays@data@listData[[1]][
                                     genes %in% mTOR_pathway_genes$gene ,1:2], 
                                     1, function(x) sum(is.na(x))<2),
                                 INS=
                                   apply(pp_norm@assays@data@listData[[1]][
                                     genes %in% mTOR_pathway_genes$gene ,3:4], 
                                     1, function(x) sum(is.na(x))<2),
                                 ISM1=
                                   apply(pp_norm@assays@data@listData[[1]][
                                     genes %in% mTOR_pathway_genes$gene ,5:6], 
                                     1, function(x) sum(is.na(x))<2))
mTOR_peptide_found$phospep <- rownames(mTOR_peptide_found)
# Make venn diagram of shared and unique phosphopeptides 
ggvenn(mTOR_peptide_found)
ggsave(file=paste0("review/plots/",paste0("TOR_complex_shared_venn_0NA",".svg")),
       width=8, height=6)
mTOR_peptide_found$pathway <- "TOR complex"

# Annotate if genes associated with muscle GO terms is detected in 
# both replicates for each treatment. 
mu_peptide_found <- data.frame(BSA=
                                 apply(pp_norm@assays@data@listData[[1]][
                                   genes %in% muscle_pathways_genes$gene ,1:2], 
                                   1, function(x) sum(is.na(x))<2),
                               INS=
                                 apply(pp_norm@assays@data@listData[[1]][
                                   genes %in% muscle_pathways_genes$gene ,3:4], 
                                   1, function(x) sum(is.na(x))<2),
                               ISM1=
                                 apply(pp_norm@assays@data@listData[[1]][
                                   genes %in% muscle_pathways_genes$gene ,5:6], 
                                   1, function(x) sum(is.na(x))<2))
mu_peptide_found$phospep <- rownames(mu_peptide_found)
ggvenn(mu_peptide_found)
ggsave(file=paste0("review/plots/",paste0("muscle_peptide_shared_venn_0NA_log",".svg")),
       width=8, height=6)
mu_peptide_found$pathway <- "Muscle"

# Combinging tables of detected genes for pathways
peptide_found <- rbind(is_peptide_found,mTOR_peptide_found,mu_peptide_found)

# saving results of detected phosphopeptides 
p_found <- peptide_found
pos <- match(p_found$phospep,pp_norm@NAMES)
p_found$gene <- pp_norm@elementMetadata$Gene_names[pos]
p_found$AA <- pp_norm@elementMetadata$Amino_acid[pos]
p_found$pos <- pp_norm@elementMetadata$Positions_within_proteins[pos]
p_found$protein <- pp_norm@elementMetadata$Protein[pos]
write.table(p_found,file=paste0("review/tables/",
  paste0("phosphopeptide_goterm_found_log",".tsv")),sep="\t")
rm(p_found)

# Listing and saving specific muscle genes in dataset also in muscle related
# go terms
goterms <- Term(GOTERM)
muscle_pathways_genes <- get_anno_genes(names(goterms[grep("muscle",goterms)]), 
                                        database = 'Mus.musculus')
muscle_genes <- gsub(";\\S+","",
  unique(pp_norm@elementMetadata$Gene_names[pp_norm@NAMES %in% 
                                           (peptide_found %>%
                                            filter(pathway=="Muscle") %>% 
                                            pull(phospep))]))
muscle_genes <- muscle_pathways_genes[muscle_pathways_genes$gene %in% 
                                            muscle_genes,]
muterm <- goterms[names(goterms) %in% muscle_genes$go_id]
muterm <- as.data.frame(muterm)
muterm$go_id <- rownames(muterm)
muscle_genes <- left_join(muscle_genes,muterm,by="go_id")
write.table(muscle_genes,file=paste0("review/tables/",
  paste0("phosphopeptide_muscle_goterm_found_log",".tsv")),sep="\t")
rm(muscle_genes,muterm)

### Making circular distribution plot
# Prepare plot
conditions <- c("BSA","INS","ISM1")
con_lim <- data.frame(xlim1=c(0,0,0),xlim2=c(sum(peptide_found$BSA),
                                             sum(peptide_found$INS),
                                             sum(peptide_found$ISM1)))
rownames(con_lim) <- conditions
path_lim <- peptide_found %>% filter(BSA==TRUE) %>% group_by(pathway) %>% 
  tally(BSA)
path_lim <- rbind(path_lim,peptide_found %>% filter(INS==TRUE) %>% 
                    group_by(pathway) %>% 
                    tally(INS))
path_lim <- rbind(path_lim,peptide_found %>% filter(ISM1==TRUE) %>% 
                    group_by(pathway) %>% 
                    tally(ISM1))
path_lim <- as.data.frame(path_lim)
path_lim$xlim1 <- 0
path_lim$pathway[1:3] <- paste0("BSA ",path_lim$pathway[1:3])
path_lim$pathway[4:6] <- paste0("INS ",path_lim$pathway[4:6])
path_lim$pathway[7:9] <- paste0("ISM1 ",path_lim$pathway[7:9])
rownames(path_lim) <- path_lim$pathway
path_lim <- path_lim[,-1]
path_lim <- path_lim[,c(2,1)]
colnames(path_lim)[2] <- "xlim2"
path_lim

shared_all <- peptide_found %>% group_by(pathway) %>% 
  filter(BSA==TRUE & INS==TRUE & ISM1==TRUE) %>% group_by(pathway)  %>% tally()
shared_BSAINS <- peptide_found %>% group_by(pathway) %>% 
  filter(BSA==TRUE & INS==TRUE & ISM1==FALSE) %>%group_by(pathway) %>% tally()
shared_INSISM1 <- peptide_found %>% group_by(pathway) %>% 
  filter(INS==TRUE & ISM1==TRUE) %>% group_by(pathway) %>% tally()

tally_count <- peptide_found[,-4] %>% group_by_all %>% dplyr::count()
tc_ism <- tally_count %>% filter(BSA==FALSE & INS==FALSE & ISM1==TRUE)
tc_ins <- tally_count %>% filter(BSA==FALSE & INS==TRUE & ISM1==FALSE)
tc_bsa <- tally_count %>% filter(BSA==TRUE & INS==FALSE & ISM1==FALSE)
tc_all <- tally_count %>% filter(BSA==TRUE & INS==TRUE & ISM1==TRUE)
tc_bsaism <- tally_count %>% filter(BSA==TRUE & INS==FALSE & ISM1==TRUE)
tc_ismins <- tally_count %>% filter(BSA==FALSE & INS==TRUE & ISM1==TRUE)
tc_bsains <- tally_count %>% filter(BSA==TRUE & INS==TRUE & ISM1==FALSE)

# Count occurances of each combination where phosphopeptide only detected in 1 treatment
uniq <- peptide_found[,-4] %>% 
  group_by_all %>% 
  dplyr::count() %>% 
  dplyr::filter(BSA+INS+ISM1==1)

# Set bed for circos plot
bed <- data.frame(pathway=c(rownames(path_lim)[7:9],rownames(path_lim)[4:6],
                            rownames(path_lim)[1:3]),
                  start=0,end=uniq$n,value1="text")

# set bed based on genes
bed2 <- data.frame()
for(i in 1:dim(uniq)[1]){
  genes <- peptide_found %>% filter(BSA==uniq$BSA[i] &
                                      INS==uniq$INS[i] &
                                      ISM1==uniq$ISM1[i] &
                                      pathway ==uniq$pathway[i]) %>% 
                                      pull(phospep)
  bed2 <- rbind(bed2,data.frame(pathway=rep(bed$pathway[i],bed$end[i]),start=seq(1,bed$end[i],1),
                                end=seq(1,bed$end[i],1),value1=genes))
}
dim(bed2)
pos <- match(bed2$value1,pp_norm@elementMetadata$name)
bed2$value1 <- paste(gsub(";\\S+","",pp_norm@elementMetadata$Gene_names[pos]),
                     paste0(pp_norm@elementMetadata$Amino_acid[pos],
                       pp_norm@elementMetadata$Positions_within_proteins[pos]))

## Drawing circular plot
plot.new()
par(mar=rep(0,4))
circos.clear()
circos.par(track.margin=c(0,0),canvas.xlim = c(-1.9, 1.9), canvas.ylim = c(-0.1, 0.1))
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.15), start.degree = 90,
           gap.degree =2)
circos.initialize(sectors=conditions,xlim=con_lim,
                  sector.width=con_lim$xlim2)
circos.trackPlotRegion(sectors = conditions, y=con_lim$xlim2,
                       track.height = 0.08,
                       bg.col = brewer.pal(9, "Set1"),
                       #bg.col = c("#EE00007D","#0000EE7D","#00EE007D"),
                       panel.fun = function(x, y) {})
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
              facing = "bending.inside", niceFacing = TRUE, col = "black")
}
circos.clear()
par(new = TRUE) # <- magic
circos.par("canvas.xlim" = c(-1.6, 1.6), "canvas.ylim" = c(-0.01, 0.01))
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.20), start.degree = 90,
           gap.degree =2)
circos.initialize(sectors=rownames(path_lim),xlim=path_lim,
                  sector.width=path_lim[,2])
circos.trackPlotRegion(sectors = rownames(path_lim), y=path_lim[,2],
                       track.height = 0.07,
                       bg.col = rep(dittoColors()[1:38][c(8,13,1)],3),
                       bg.border = "gray",
                       track.index = 1,
                       panel.fun = function(x, y) {
                         circos.axis(labels=NULL,
                                     major.at = c(0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750),
                                     major.tick=TRUE,
                                     minor.ticks = 10,
                                     labels.facing="clockwise",
                                     h = "bottom",
                                     direction ="outside",major.tick.length=400,
                                     labels.cex = 1, col="white")
                       })

pos <- 0.64
# Linking all shared muscle - ism to ins
circos.link(sector.index1=paste("ISM1",tc_all$pathway[2],sep=" "), 
            c(tc_ism$n[2]+tc_ismins$n[2],
              tc_ism$n[2]+tc_ismins$n[2]+tc_all$n[2]), 
            sector.index2=paste("INS",tc_all$pathway[2],sep=" "), 
            c(tc_ins$n[2]+tc_bsains$n[2],
              tc_ins$n[2]+tc_bsains$n[2]+tc_all$n[2]),
            col="grey90", border = NA,rou = pos)
# Linking all shared - bsa to ism1
circos.link(sector.index1=paste("BSA",tc_all$pathway[2],sep=" "), 
            c(tc_bsa$n[2]+tc_bsaism$n[2],
              tc_bsa$n[2]+tc_bsaism$n[2]+tc_all$n[2]), 
            sector.index2=paste("ISM1",tc_all$pathway[2],sep=" "), 
            c(tc_ism$n[2]+tc_ismins$n[2],
              tc_ism$n[2]+tc_ismins$n[2]+tc_all$n[2]),
            col="grey90", border = NA,rou = pos)
# Linking all shared - bsa to ins
circos.link(sector.index1=paste("BSA",tc_all$pathway[2],sep=" "), 
            c(tc_bsa$n[2]+tc_bsaism$n[2],
              tc_bsa$n[2]+tc_bsaism$n[2]+tc_all$n[2]), 
            sector.index2=paste("INS",tc_all$pathway[2],sep=" "), 
            c(tc_ins$n[2]+tc_bsains$n[2],
              tc_ins$n[2]+tc_bsains$n[2]+tc_all$n[2]),
            col="grey90", border = NA,rou = pos)

#col2rgb("grey50")
# Linking all shared insulin receptor pathway - ism to ins
circos.link(sector.index1=paste("ISM1",tc_all$pathway[1],sep=" "), 
            c(tc_ism$n[1]+tc_ismins$n[1],
              tc_ism$n[1]+tc_ismins$n[1]+tc_all$n[1]), 
            sector.index2=paste("INS",tc_all$pathway[1],sep=" "), 
            c(tc_ins$n[1]+tc_bsains$n[1],
              tc_ins$n[1]+tc_bsains$n[1]+tc_all$n[1]),
            col="grey80", border = NA,rou = pos)
# Linking all shared - bsa to ism1
circos.link(sector.index1=paste("BSA",tc_all$pathway[1],sep=" "), 
            c(tc_bsa$n[1]+tc_bsaism$n[1],
              tc_bsa$n[1]+tc_bsaism$n[1]+tc_all$n[1]), 
            sector.index2=paste("ISM1",tc_all$pathway[1],sep=" "), 
            c(tc_ism$n[1]+tc_ismins$n[1],
              tc_ism$n[1]+tc_ismins$n[1]+tc_all$n[1]),
            col="grey80", border = NA,rou = pos)
# Linking all shared - bsa to ins
circos.link(sector.index1=paste("BSA",tc_all$pathway[1],sep=" "), 
            c(tc_bsa$n[1]+tc_bsaism$n[1],
              tc_bsa$n[1]+tc_bsaism$n[1]+tc_all$n[1]), 
            sector.index2=paste("INS",tc_all$pathway[1],sep=" "), 
            c(tc_ins$n[1]+tc_bsains$n[1],
              tc_ins$n[1]+tc_bsains$n[1]+tc_all$n[1]),
            col="grey80", border = NA,rou = pos)

# Linking all shared TOR complex - ism to ins
circos.link(sector.index1=paste("ISM1",tc_all$pathway[3],sep=" "), 
            c(tc_ism$n[3],
              tc_ism$n[3]+tc_all$n[3]), 
            sector.index2=paste("INS",tc_all$pathway[3],sep=" "), 
            c(tc_ins$n[3]+0,
              tc_ins$n[3]+0+tc_all$n[3]),
            col="grey70", border = NA,rou = pos)
# Linking all shared - bsa to ism1
circos.link(sector.index1=paste("BSA",tc_all$pathway[3],sep=" "), 
            c(tc_bsa$n[3]+tc_bsaism$n[3],
              tc_bsa$n[3]+tc_bsaism$n[3]+tc_all$n[3]), 
            sector.index2=paste("ISM1",tc_all$pathway[3],sep=" "), 
            c(tc_ism$n[3]+0,
              tc_ism$n[3]+0+tc_all$n[3]),
            col="grey70", border = NA,rou = pos)
# Linking all shared - bsa to ins
circos.link(sector.index1=paste("BSA",tc_all$pathway[3],sep=" "), 
            c(0+tc_bsaism$n[3],
              0+tc_bsaism$n[3]+tc_all$n[3]), 
            sector.index2=paste("INS",tc_all$pathway[3],sep=" "), 
            c(tc_ins$n[3]+0,
              tc_ins$n[3]+0+tc_all$n[3]),
            col="grey70", border = NA,rou = pos)

# Linking ins+ism shared only
circos.link(sector.index1=paste("INS",tc_all$pathway[2],sep=" "), 
            c(tc_ins$n[2]+tc_bsains$n[2]+tc_all$n[2],
              tc_ins$n[2]+tc_bsains$n[2]+tc_all$n[2]+tc_ismins$n[2]), 
            sector.index2=paste("ISM1",tc_all$pathway[2],sep=" "), 
            c(tc_ism$n[2],tc_ism$n[2]+tc_ismins$n[2]),
            col=viridis(10,alpha = 0.4)[2], border = NA,rou = pos)
# Linking bsa+ism shared only
circos.link(sector.index1=paste("BSA",tc_all$pathway[2],sep=" "), 
            c(tc_bsa$n[2],tc_bsa$n[2]+tc_bsaism$n[2]), 
            sector.index2=paste("ISM1",tc_all$pathway[2],sep=" "), 
            c(tc_ism$n[2]+tc_ismins$n[2]+tc_all$n[2],
              tc_ism$n[2]+tc_ismins$n[2]+tc_all$n[2]+tc_bsaism$n[2]),
            col=viridis(10,alpha = 0.4)[4], border = NA,rou = pos)
# Linking bsa+ins only
circos.link(sector.index1=paste("BSA",tc_all$pathway[2],sep=" "), 
            c(tc_bsa$n[2]+tc_bsaism$n[2]+tc_all$n[2],
              tc_bsa$n[2]+tc_bsaism$n[2]+tc_all$n[2]+tc_bsains$n[2]), 
            sector.index2=paste("INS",tc_all$pathway[2],sep=" "), 
            c(tc_ins$n[2],
              tc_ins$n[2]+tc_bsains$n[2]),
            col=viridis(10,alpha = 0.4)[7], border = NA,rou = pos)

# Linking ins+ism shared only
circos.link(sector.index1=paste("INS",tc_all$pathway[1],sep=" "), 
            c(tc_ins$n[1]+tc_bsains$n[1]+tc_all$n[1],
              tc_ins$n[1]+tc_bsains$n[1]+tc_all$n[1]+tc_ismins$n[1]), 
            sector.index2=paste("ISM1",tc_all$pathway[1],sep=" "), 
            c(tc_ism$n[1],tc_ism$n[1]+tc_ismins$n[1]),
            col=viridis(10,alpha = 0.4)[2], border = NA,rou = pos)
# Linking bsa+ism shared only
circos.link(sector.index1=paste("BSA",tc_all$pathway[1],sep=" "), 
            c(tc_bsa$n[1],tc_bsa$n[1]+tc_bsaism$n[1]), 
            sector.index2=paste("ISM1",tc_all$pathway[1],sep=" "), 
            c(tc_ism$n[1]+tc_ismins$n[1]+tc_all$n[1],
              tc_ism$n[1]+tc_ismins$n[1]+tc_all$n[1]+tc_bsaism$n[1]),
            col=viridis(10,alpha = 0.4)[4], border = NA,rou = pos)
# Linking bsa+ins only
circos.link(sector.index1=paste("BSA",tc_all$pathway[1],sep=" "), 
            c(tc_bsa$n[1]+tc_bsaism$n[1]+tc_all$n[1],
              tc_bsa$n[1]+tc_bsaism$n[1]+tc_all$n[1]+tc_bsains$n[1]), 
            sector.index2=paste("INS",tc_all$pathway[1],sep=" "), 
            c(tc_ins$n[1],
              tc_ins$n[1]+tc_bsains$n[1]),
            col=viridis(10,alpha = 0.4)[7], border = NA,rou = pos)

# Linking ins+ism shared only
circos.link(sector.index1=paste("INS",tc_all$pathway[3],sep=" "), 
            c(tc_ins$n[3]+0+tc_all$n[3],
              tc_ins$n[3]+0+tc_all$n[3]+tc_ismins$n[3]), 
            sector.index2=paste("ISM1",tc_all$pathway[3],sep=" "), 
            c(tc_ism$n[3],tc_ism$n[3]+tc_ismins$n[3]),
            col=viridis(10,alpha = 0.4)[2], border = NA,rou = pos)
# Linking bsa+ism shared only
circos.link(sector.index1=paste("BSA",tc_all$pathway[3],sep=" "), 
            c(0,tc_bsa$n[3]+tc_bsaism$n[3]), 
            sector.index2=paste("ISM1",tc_all$pathway[3],sep=" "), 
            c(tc_ism$n[3]+0+tc_all$n[3],
              tc_ism$n[3]+0+tc_all$n[3]+tc_bsaism$n[3]),
            col=viridis(10,alpha = 0.4)[4], border = NA,rou = pos)
# Linking bsa+ins only
circos.link(sector.index1=paste("BSA",tc_all$pathway[3],sep=" "), 
            c(0+tc_bsaism$n[3]+tc_all$n[3],
              0+tc_bsaism$n[3]+tc_all$n[3]), 
            sector.index2=paste("INS",tc_all$pathway[3],sep=" "), 
            c(tc_ins$n[3],
              tc_ins$n[3]+0),
            col=viridis(10,alpha = 0.4)[7], border = NA,rou = pos)

circos.clear()
par(new = TRUE) # <- magic
circos.par("canvas.xlim" = c(-2.9, 2.9), "canvas.ylim" = c(-0.01, 0.01))
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.20), start.degree = 90,
           gap.degree =2)
circos.initialize(sectors=rownames(path_lim),xlim=path_lim,
                  sector.width=path_lim[,2])
circos.genomicTrackPlotRegion(bed2,ylim = c(0, 1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicText(region, value, y = 22, labels.column = 1, #y=9.5
                                                   facing = "clockwise", adj = c(0, 0),
                                                   cex = 0.7,niceFacing = TRUE, # Previous cex=0.7
                                                   posTransform = posTransform.text,
                                                   padding=0.2, # Previous padding -0.17
                                                   extend = 0.0,col="grey20") # Previous extend=0
                              }, track.height = 0.1, bg.border = NA)
circos.clear()
par(new = TRUE) # <- magic
circos.par("canvas.xlim" = c(-1.77, 1.77), "canvas.ylim" = c(-0.01, 0.01))
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.20), start.degree = 90,
           gap.degree =2)
circos.initialize(sectors=rownames(path_lim),xlim=path_lim,
                  sector.width=path_lim[,2])
circos.genomicPosTransformLines(bed2, track.height = 0.12,
                                posTransform = function(region, value) 
                                  posTransform.text(region, y = 0,
                                                    labels = value[[1]], 
                                                    cex = 0.7, padding = -0.27),
                                direction = "outside")

# Making circular plot legend
lgd_lines1 = Legend(at = c(" TOR complex", 
                           " All incl. term ?muscle?",
                           " Insulin receptor signalling"), type = "lines", 
                    legend_gp = gpar(col = dittoColors()[1:38][c(8,13,1)], lwd = 10), 
                    title_position = "topleft", 
                    title = "GO Pathway category",
                    labels_gp = gpar(fontsize = 9),
                    title_gp = gpar(fontsize = 9))
lgd_lines2 = Legend(at = c(" BSA", " INS"," ISM1"), type = "lines", 
                    legend_gp = gpar(col=brewer.pal(9, "Set1")[1:3], lwd = 10), 
                    title_position = "topleft", 
                    title = "Condition",
                    labels_gp = gpar(fontsize = 9),
                    title_gp = gpar(fontsize = 9))
lgd_lines3 = Legend(at = c(" all conditions, muscle related pathways", 
                           " all conditions, IR signalling pathway",
                           " all conditions, TOR complex pathway",
                           " not in BSA",
                           " not in INS",
                           " not in ISM1"), type = "lines", 
                    legend_gp = gpar(col = c(ColToHex("grey90"),
                                             ColToHex("grey80"),
                                             ColToHex("grey70"),
                                             viridis(10)[2],
                                             viridis(10)[4],
                                             viridis(10)[7]), lwd = 10), 
                    title_position = "topleft", 
                    title = "Shared phosphopeptides",
                    labels_gp = gpar(fontsize = 9),
                    title_gp = gpar(fontsize = 9))
lgd_list_vertical = packLegend(lgd_lines1, lgd_lines2, lgd_lines3)
lgd_list_vertical
draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
dev.copy2pdf(file=paste0("plots/circular_plot_0NA_extraTOR",paste0(".pdf")), width=25, height=25)
dev.off()
svg("plots/2022-09-01_circular_extraTOR.svg")

### Making bar plots of normalized intensities for selected phosphopeptides
# and phosphopeptides belonging to genes related to the GO terms: insulin
# receptor signalling, TOR pathway and any term with muscle.

# Using filtered data
pp_norm <- readRDS(file="pp_norm_log_filter.rds")

# Making bar plot of insulin receptor phosphorylation between conditions
Insr <- pp_norm@assays@data@listData[[1]][grep("Insr",
                                       rownames(pp_norm@assays@data@listData[[1]])),]
Insr <- as.data.frame(Insr)
Insr$condition <- gsub("_\\S","",rownames(Insr))
min_norm_val <- min(pp_norm@assays@data[[1]],na.rm=TRUE)
Insr$Insr <- Insr$Insr - min_norm_val
# If missing set value to zero
Insr[is.na(Insr)] <- 0
Insr$Insr <- Insr$Insr/max(Insr$Insr,na.rm=TRUE)
Insr <- Insr[c(1,2,5,6,3,4),]
Insr$condition <- factor(Insr$condition,levels = c("BSA", "ISM1", "INS"))
cond <- c("BSA","BSA","ISM1","ISM1","INS","INS")
ggplot(Insr, aes(x=condition,y=Insr,fill=condition,width=.8)) + 
  geom_bar(position = 'dodge', stat = 'summary') +
  scale_fill_manual(values = brewer.pal(9, "Set1")[c(1,3,2)]) +
  geom_dotplot(binaxis='y', stackdir='center') +
  ylab("Fraction of maximum normalized intensity") +
  theme_classic() + 
  theme(axis.title.x=element_blank(),
                     axis.text.x = element_text(angle = 45, 
                                                vjust = 1, hjust=1))
ggsave(paste0("review/plots/",paste0("Insr_bar_percentmax_filter",".jpeg")),
       units="cm", width=8, height=10, dpi=600)
ggsave(file=paste0("review/plots/",paste0("Insr_bar_percentmax_filter",".svg")),
       width=8, height=10)

# Making barplot of Rps6ka5 phosphosite 375
phospos <- 375
as.data.frame(pp_norm@elementMetadata[grep("Rps6ka5",
  rownames(pp_norm@assays@data@listData[[1]])),]) %>% 
  filter(Positions_within_proteins==phospos)
Rps6ka5 <- pp_norm@assays@data@listData[[1]][pp_norm@elementMetadata$Positions_within_proteins == phospos,]
Rps6ka5 <- Rps6ka5[grep("Rps6ka5",rownames(Rps6ka5)),]
Rps6ka5 <- as.data.frame(Rps6ka5)
Rps6ka5$condition <- gsub("_\\S","",rownames(Rps6ka5))
min_norm_val <- min(pp_norm@assays@data[[1]],na.rm=TRUE)
Rps6ka5$Rps6ka5 <- Rps6ka5$Rps6ka5 - min_norm_val
# If missing set value to zero
Rps6ka5[is.na(Rps6ka5)] <- 0
Rps6ka5$Rps6ka5 <- Rps6ka5$Rps6ka5/max(Rps6ka5$Rps6ka5,na.rm=TRUE)
Rps6ka5 <- Rps6ka5[c(1,2,5,6,3,4),]
Rps6ka5$condition <- factor(Rps6ka5$condition,levels = c("BSA", "ISM1", "INS"))
cond <- c("BSA","BSA","ISM1","ISM1","INS","INS")
ggplot(Rps6ka5, aes(x=condition,y=Rps6ka5,fill=condition,width=.8)) + 
  geom_bar(position = 'dodge', stat = 'summary') +
  scale_fill_manual(values = brewer.pal(9, "Set1")[c(1,3,2)]) +
  geom_dotplot(binaxis='y', stackdir='center') +
  ylab("Fraction of maximum normalized intensity") +
  ggtitle(paste("Rps6ka5",phospos,sep=" "))+
  theme_classic() + 
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   vjust = 1, hjust=1))
ggsave(paste0("review/plots/",paste0("Rps6ka5_bar_percentmax_filter",".jpeg")),
       units="cm", width=8, height=10, dpi=600)
ggsave(file=paste0("review/plots/",paste0("Rps6ka5_bar_percentmax_filter",".svg")),
       width=8, height=10)

# Making barplot of Rps3a phosphosite 237
phospos <- 237
as.data.frame(pp_norm@elementMetadata[grep("Rps3a",
  rownames(pp_norm@assays@data@listData[[1]])),]) %>% 
  filter(Positions_within_proteins==phospos)
Rps3a <- pp_norm@assays@data@listData[[1]][pp_norm@elementMetadata$Positions_within_proteins == 237,]
Rps3a <- Rps3a[grep("Rps3a",rownames(Rps3a)),]
Rps3a <- as.data.frame(Rps3a)
Rps3a$condition <- gsub("_\\S","",rownames(Rps3a))
min_norm_val <- min(pp_norm@assays@data[[1]],na.rm=TRUE)
Rps3a$Rps3a <- Rps3a$Rps3a - min_norm_val
# If missing set value to zero
Rps3a[is.na(Rps3a)] <- 0
Rps3a$Rps3a <- Rps3a$Rps3a/max(Rps3a$Rps3a,na.rm=TRUE)
Rps3a <- Rps3a[c(1,2,5,6,3,4),]
Rps3a$condition <- factor(Rps3a$condition,levels = c("BSA", "ISM1", "INS"))
cond <- c("BSA","BSA","ISM1","ISM1","INS","INS")
ggplot(Rps3a, aes(x=condition,y=Rps3a,fill=condition,width=.8)) + 
  geom_bar(position = 'dodge', stat = 'summary') +
  scale_fill_manual(values = brewer.pal(9, "Set1")[c(1,3,2)]) +
  geom_dotplot(binaxis='y', stackdir='center') +
  ylab("Fraction of maximum normalized intensity") +
  ggtitle(paste("Rps3a",phospos,sep=" "))+
  theme_classic() + 
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   vjust = 1, hjust=1))
ggsave(paste0("review/plots/",paste0("Rps3a_237_bar_percentmax_filter",".jpeg")),
       units="cm", width=8, height=10, dpi=600)
ggsave(file=paste0("review/plots/",paste0("Rps3a_237_bar_percentmax_filter",".svg")),
       width=8, height=10)

# Making barplot of Rps3a phosphosite 238
phospos <- 238
as.data.frame(pp_norm@elementMetadata[grep("Rps3a",
  rownames(pp_norm@assays@data@listData[[1]])),]) %>% 
  filter(Positions_within_proteins==phospos)
Rps3a <- pp_norm@assays@data@listData[[1]][pp_norm@elementMetadata$Positions_within_proteins == phospos,]
Rps3a <- Rps3a[grep("Rps3a",rownames(Rps3a)),]
Rps3a <- as.data.frame(Rps3a)
Rps3a$condition <- gsub("_\\S","",rownames(Rps3a))
min_norm_val <- min(pp_norm@assays@data[[1]],na.rm=TRUE)
Rps3a$Rps3a <- Rps3a$Rps3a - min_norm_val
# If missing set value to zero
Rps3a[is.na(Rps3a)] <- 0
Rps3a$Rps3a <- Rps3a$Rps3a/max(Rps3a$Rps3a,na.rm=TRUE)
Rps3a <- Rps3a[c(1,2,5,6,3,4),]
Rps3a$condition <- factor(Rps3a$condition,levels = c("BSA", "ISM1", "INS"))
cond <- c("BSA","BSA","ISM1","ISM1","INS","INS")
ggplot(Rps3a, aes(x=condition,y=Rps3a,fill=condition,width=.8)) + 
  geom_bar(position = 'dodge', stat = 'summary') +
  scale_fill_manual(values = brewer.pal(9, "Set1")[c(1,3,2)]) +
  geom_dotplot(binaxis='y', stackdir='center') +
  ylab("Fraction of maximum normalized intensity") +
  ggtitle(paste("Rps3a",phospos,sep=" "))+
  theme_classic() + 
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   vjust = 1, hjust=1))
ggsave(paste0("review/plots/",paste0("Rps3a_238_bar_percentmax_filter",".jpeg")),
       units="cm", width=8, height=10, dpi=600)
ggsave(file=paste0("review/plots/",paste0("Rps3a_238_bar_percentmax_filter",".svg")),
       width=8, height=10)

# Making barplot of Rpl13 phosphosite 77
phospos <- 77
as.data.frame(pp_norm@elementMetadata[grep("Rpl13",
  rownames(pp_norm@assays@data@listData[[1]])),]) %>% 
  filter(Positions_within_proteins==phospos)
Rpl13 <- pp_norm@assays@data@listData[[1]][pp_norm@elementMetadata$Positions_within_proteins == phospos,]
Rpl13 <- Rpl13[grep("Rpl13",rownames(Rpl13)),]
Rpl13 <- as.data.frame(Rpl13)
Rpl13$condition <- gsub("_\\S","",rownames(Rpl13))
min_norm_val <- min(pp_norm@assays@data[[1]],na.rm=TRUE)
Rpl13$Rpl13 <- Rpl13$Rpl13 - min_norm_val
# If missing set value to zero
Rpl13[is.na(Rpl13)] <- 0
Rpl13$Rpl13 <- Rpl13$Rpl13/max(Rpl13$Rpl13,na.rm=TRUE)
Rpl13 <- Rpl13[c(1,2,5,6,3,4),]
Rpl13$condition <- factor(Rpl13$condition,levels = c("BSA", "ISM1", "INS"))
cond <- c("BSA","BSA","ISM1","ISM1","INS","INS")
ggplot(Rpl13, aes(x=condition,y=Rpl13,fill=condition,width=.8)) + 
  geom_bar(position = 'dodge', stat = 'summary') +
  scale_fill_manual(values = brewer.pal(9, "Set1")[c(1,3,2)]) +
  geom_dotplot(binaxis='y', stackdir='center') +
  ylab("Fraction of maximum normalized intensity") +
  ggtitle(paste("Rpl13",phospos,sep=" "))+
  theme_classic() + 
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   vjust = 1, hjust=1))
ggsave(paste0("review/plots/",paste0("Rpl13_77_bar_percentmax_filter",".jpeg")),
       units="cm", width=8, height=10, dpi=600)
ggsave(file=paste0("review/plots/",paste0("Rpl13_77_bar_percentmax_filter",".svg")),
       width=8, height=10)

# Making barplot of Eif3 phosphosite 584
phospos <- 584
as.data.frame(pp_norm@elementMetadata[grep("Eif3",
  rownames(pp_norm@assays@data@listData[[1]])),]) %>% 
  filter(Positions_within_proteins==phospos)
Eif3a <- pp_norm@assays@data@listData[[1]][pp_norm@elementMetadata$Positions_within_proteins == phospos,]
Eif3a <- Eif3a[grep("Eif3a",rownames(Eif3a)),]
Eif3a <- as.data.frame(Eif3a)
Eif3a$condition <- gsub("_\\S","",rownames(Eif3a))
min_norm_val <- min(pp_norm@assays@data[[1]],na.rm=TRUE)
Eif3a$Eif3a <- Eif3a$Eif3a - min_norm_val
# If missing set value to zero
Eif3a[is.na(Eif3a)] <- 0
Eif3a$Eif3a <- Eif3a$Eif3a/max(Eif3a$Eif3a,na.rm=TRUE)
Eif3a <- Eif3a[c(1,2,5,6,3,4),]
Eif3a$condition <- factor(Eif3a$condition,levels = c("BSA", "ISM1", "INS"))
cond <- c("BSA","BSA","ISM1","ISM1","INS","INS")
ggplot(Eif3a, aes(x=condition,y=Eif3a,fill=condition,width=.8)) + 
  geom_bar(position = 'dodge', stat = 'summary') +
  scale_fill_manual(values = brewer.pal(9, "Set1")[c(1,3,2)]) +
  geom_dotplot(binaxis='y', stackdir='center') +
  ylab("Fraction of maximum normalized intensity") +
  ggtitle(paste("Eif3a",phospos,sep=" "))+
  theme_classic() + 
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   vjust = 1, hjust=1))
ggsave(paste0("review/plots/",paste0("Eif3a_584_bar_percentmax_filter",".jpeg")),
       units="cm", width=8, height=10, dpi=600)
ggsave(file=paste0("review/plots/",paste0("Eif3a_584_bar_percentmax_filter",".svg")),
       width=8, height=10)

###  Making heatmap of insulin receptor substrate phosphorylation
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

# Calculating center log ratio of normalized intensities
clr_int <- clr(as.data.frame(pp_norm@assays@data[[1]]))
colnames(clr_int) <- unname(as.data.frame(strsplit(colnames(clr_int),"_"))[1,])
rn <- paste(gsub(";\\S+","",pp_norm@elementMetadata$Gene_names),
            paste0(pp_norm@elementMetadata$Amino_acid,
                   gsub(";\\S+","",
                        pp_norm@elementMetadata$Positions_within_proteins)))
irs <- grep("Irs",pp_norm@elementMetadata$Gene_names)

d <- clr_int[irs,3:6][is.na(clr_int[irs,3:6])] <-0
# Make heatmap of all insulin-receptor substrate phosphorylation between 
# INS and ISM1. Set non-detectables to grey 
mat_cluster_rows <- sort_hclust(hclust(dist(clr_int[irs,3:6],
                      method="manhattan"),method="ward.D"))
clr_int[is.na(assay(pp_norm))] <- NA
ha <- HeatmapAnnotation(condition = c("INS","INS","ISM1","ISM1"),
                        col = list(condition= c(
                                                INS=brewer.pal(9, "Set1")[2],
                                                ISM1=brewer.pal(9, "Set1")[3])),
                        annotation_legend_param = list(condition = list(
                          title = "condition",
                          at = c("INS","ISM1"),
                          labels = c("INS","ISM1"))))
h<- Heatmap(clr_int[irs,3:6],
        na_col = "grey",#clustering_distance_columns ="canberra",
        cluster_rows=mat_cluster_rows,
        col = circlize::colorRamp2(
          seq(-0.1, 0.1, (0.1/5)),
          rev(RColorBrewer::brewer.pal(11, "RdBu"))),#column_km = 2,
        top_annotation = ha,
        row_dend_width=unit(10, "cm"),column_dend_height=unit(2, "cm"),
        show_column_names = FALSE,
        heatmap_legend_param=list(title = "Log2 ratio"),
        row_labels =rn[irs])
h
grid.ls(grid.force()) # "col_annotation" looks like it's the one to edit
grid.gedit("GRID.rect", gp = gpar(col="grey70"))
dev.copy2pdf(file=paste0("review/plots/irs_heatmap_grid_na_filter",
                         paste0(".pdf")), width=7.2, height=4)
dev.off()
h
dev.copy2pdf(file=paste0("review/plots/irs_heatmap_nogrid_no_na_filter",
                         paste0(".pdf")), width=7.2, height=4)
dev.off()

# Making heatmap of phosphopeptides belonging to the 
# insulin receptor signalling pathway that are present in at least one 
# ISM1 and INS replicate, and where there is a significant difference 
# between the two conditions
#is_phospep <- peptide_found %>% 
#  filter(INS==TRUE & ISM1==TRUE) %>% 
#  filter(pathway=="IR signalling") %>%
#  pull(phospep)

is_phospep <- peptide_found$phospep[peptide_found$pathway=="IR signalling"]
is_norm_intens <- pp_norm@assays@data[[1]][rownames(pp_norm@assays@data[[1]]) %in% 
  is_phospep,]
is_sig <- pp_norm@NAMES[pp_norm@elementMetadata$padj.vals.ISMvsINS < 0.05 & rownames(pp_norm@assays@data[[1]]) %in% 
  rownames(is_norm_intens)]
is_norm_intens <- is_norm_intens[rownames(is_norm_intens) %in% is_sig,]
is_norm_intens <- is_norm_intens-min_norm_val
# If missing set value to zero
is_norm_intens[is.na(is_norm_intens)] <- 0
is_norm_intens <- apply(is_norm_intens,1, function(x) {x/max(x)})
pos <- match(colnames(is_norm_intens),pp_norm@elementMetadata$name)
# Tag phosphosite
phospos <- paste(gsub(";\\S+","",pp_norm@elementMetadata$Gene_names[pos]),
                     paste0(pp_norm@elementMetadata$Amino_acid[pos],
                            pp_norm@elementMetadata$Positions_within_proteins[pos]))
is_norm_intens <- as.data.frame(is_norm_intens)
is_norm_intens$condition <- gsub("_\\S","",rownames(is_norm_intens))
is_norm_intens$condition <- factor(is_norm_intens$condition,levels = c("BSA", "ISM1", "INS"))

# Loop that makes bar graphs of significant phosphopeptides present in 
# insulin receptor signalling GO term
for(i in 1:(ncol(is_norm_intens)-1)){
  p <- ggplot(is_norm_intens, aes(x=condition,y=is_norm_intens[,i],fill=condition,width=.8)) + 
    geom_bar(position = 'dodge', stat = 'summary') +
    scale_fill_manual(values = brewer.pal(9, "Set1")[c(1,3,2)]) +
    geom_dotplot(binaxis='y', stackdir='center',binwidth = 0.05) +
    ylab("Fraction of maximum normalized intensity") +
    ggtitle(phospos[i])+
    theme_classic() + 
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 45, 
                                     vjust = 1, hjust=1))
  ggsave(p,file=paste0("review/plots/bargraphs/filtered/is/",paste0(paste0(colnames(is_norm_intens)[i]),"_bar_percentmax"),".jpeg"),
         units="cm", width=8, height=10, dpi=600)
  ggsave(p,file=paste0("review/plots/bargraphs/filtered/is/",paste0(paste0(colnames(is_norm_intens)[i]),"_bar_percentmax"),".svg"),
         width=8, height=10)
}

# Making heatmap of phosphopeptides belonging to 
# muscle pathways that are present in at least one 
# ISM1 and INS replicate, and where there is a significant difference 
# between the two conditions
#mu_phospep <- peptide_found %>% 
#  filter(INS==TRUE & ISM1==TRUE) %>% 
#  filter(pathway=="Muscle") %>%
#  pull(phospep)

mu_phospep <- peptide_found$phospep[peptide_found$pathway=="Muscle"]
mu_norm_intens <- pp_norm@assays@data[[1]][rownames(pp_norm@assays@data[[1]]) %in% mu_phospep,]
mu_sig <- pp_norm@NAMES[pp_norm@elementMetadata$padj.vals.ISMvsINS < 0.05 &  rownames(pp_norm@assays@data[[1]]) %in% rownames(mu_norm_intens)]
mu_norm_intens <- mu_norm_intens[rownames(mu_norm_intens) %in% mu_sig,]
mu_norm_intens <- mu_norm_intens-min_norm_val
# If missing set value to zero
mu_norm_intens[is.na(mu_norm_intens)] <- 0
mu_norm_intens <- apply(mu_norm_intens,1, function(x) {x/max(x)})
pos <- match(colnames(mu_norm_intens),pp_norm@elementMetadata$name)
# Tag phosphosite
phospos <- paste(gsub(";\\S+","",pp_norm@elementMetadata$Gene_names[pos]),
                 paste0(pp_norm@elementMetadata$Amino_acid[pos],
                        pp_norm@elementMetadata$Positions_within_proteins[pos]))
mu_norm_intens <- as.data.frame(mu_norm_intens)
mu_norm_intens$condition <- gsub("_\\S","",rownames(mu_norm_intens))
mu_norm_intens$condition <- factor(mu_norm_intens$condition,levels = c("BSA", "ISM1", "INS"))

# Loop that makes bar graphs of significant phosphopeptides for genes present 
# in GO terms with "muscle" included
for(i in 1:(ncol(mu_norm_intens)-1)){
  p <- ggplot(mu_norm_intens, aes(x=condition,y=mu_norm_intens[,i],fill=condition,width=.8)) + 
    geom_bar(position = 'dodge', stat = 'summary') +
    scale_fill_manual(values = brewer.pal(9, "Set1")[c(1,3,2)]) +
    geom_dotplot(binaxis='y', stackdir='center',binwidth = 0.05) +
    ylab("Fraction of maximum normalized intensity") +
    ggtitle(phospos[i])+
    theme_classic() + 
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 45, 
                                     vjust = 1, hjust=1))
  ggsave(p,file=paste0("review/plots/bargraphs/filtered/mu/ism_vs_ins/",paste0(paste0(colnames(mu_norm_intens)[i]),"_bar_percentmax"),".jpeg"),
         units="cm", width=8, height=10, dpi=600)
  ggsave(p,file=paste0("review/plots/bargraphs/filtered/mu/ism_vs_ins/",paste0(paste0(colnames(mu_norm_intens)[i]),"_bar_percentmax"),".svg"),
         width=8, height=10)
}

# Making heatmap of phosphopeptides belonging to 
# muscle pathways that are present in at least one 
# ISM1 and BSA replicate, and where there is a significant difference 
# between the two conditions
#mu_phospep <- peptide_found %>% 
#  filter(BSA==TRUE & ISM1==TRUE) %>% 
#  filter(pathway=="Muscle") %>%
#  pull(phospep)

mu_phospep <- peptide_found$phospep[peptide_found$pathway=="Muscle"]
mu_norm_intens <- pp_norm@assays@data[[1]][rownames(pp_norm@assays@data[[1]]) %in% mu_phospep,]
mu_sig <- pp_norm@NAMES[pp_norm@elementMetadata$padj.vals.ISMvsBSA < 0.05 &  rownames(pp_norm@assays@data[[1]]) %in% rownames(mu_norm_intens)]
mu_norm_intens <- mu_norm_intens[rownames(mu_norm_intens) %in% mu_sig,]
mu_norm_intens <- mu_norm_intens-min_norm_val
# If missing set value to zero
mu_norm_intens[is.na(mu_norm_intens)] <- 0
mu_norm_intens <- apply(mu_norm_intens,1, function(x) {x/max(x)})
pos <- match(colnames(mu_norm_intens),pp_norm@elementMetadata$name)
# Tag phosphosite
phospos <- paste(gsub(";\\S+","",pp_norm@elementMetadata$Gene_names[pos]),
                 paste0(pp_norm@elementMetadata$Amino_acid[pos],
                        pp_norm@elementMetadata$Positions_within_proteins[pos]))
mu_norm_intens <- as.data.frame(mu_norm_intens)
mu_norm_intens$condition <- gsub("_\\S","",rownames(mu_norm_intens))
mu_norm_intens$condition <- factor(mu_norm_intens$condition,levels = c("BSA", "ISM1", "INS"))

# Loop that makes bar graphs of significant phosphopeptides for genes present 
# in GO terms with "muscle" included
for(i in 1:(ncol(mu_norm_intens)-1)){
  p <- ggplot(mu_norm_intens, aes(x=condition,y=mu_norm_intens[,i],fill=condition,width=.8)) + 
    geom_bar(position = 'dodge', stat = 'summary') +
    scale_fill_manual(values = brewer.pal(9, "Set1")[c(1,3,2)]) +
    geom_dotplot(binaxis='y', stackdir='center',binwidth = 0.05) +
    ylab("Fraction of maximum normalized intensity") +
    ggtitle(phospos[i])+
    theme_classic() + 
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 45, 
                                     vjust = 1, hjust=1))
  ggsave(p,file=paste0("review/plots/bargraphs/filtered/mu/ism_vs_bsa/",paste0(paste0(colnames(mu_norm_intens)[i]),"_bar_percentmax"),".jpeg"),
         units="cm", width=8, height=10, dpi=600)
  ggsave(p,file=paste0("review/plots/bargraphs/filtered/mu/ism_vs_bsa/",paste0(paste0(colnames(mu_norm_intens)[i]),"_bar_percentmax"),".svg"),
         width=8, height=10)
}

# Making heatmap of phosphopeptides belonging to the
# TOR complex pathway that are present in at least one 
# ISM1 and INS replicate, and where there is a significant difference 
# between the two conditions
#tor_phospep <- peptide_found %>% 
#  filter(INS==TRUE & ISM1==TRUE) %>% 
#  filter(pathway=="TOR complex") %>%
#  pull(phospep)

tor_phospep <- peptide_found$phospep[peptide_found$pathway=="TOR complex"]
tor_norm_intens <- pp_norm@assays@data[[1]][rownames(pp_norm@assays@data[[1]]) %in% tor_phospep,]
tor_sig <- pp_norm@NAMES[pp_norm@elementMetadata$padj.vals.ISMvsINS < 0.05 &  rownames(pp_norm@assays@data[[1]]) %in% rownames(tor_norm_intens)]
tor_norm_intens <- tor_norm_intens[rownames(tor_norm_intens) %in% tor_sig,]
tor_norm_intens <- tor_norm_intens-min_norm_val
# If missing set value to zero
tor_norm_intens[is.na(tor_norm_intens)] <- 0
tor_norm_intens <- apply(tor_norm_intens,1, function(x) {x/max(x)})
pos <- match(colnames(tor_norm_intens),pp_norm@elementMetadata$name)
# Tag phosphosite
phospos <- paste(gsub(";\\S+","",pp_norm@elementMetadata$Gene_names[pos]),
                 paste0(pp_norm@elementMetadata$Amino_acid[pos],
                        pp_norm@elementMetadata$Positions_within_proteins[pos]))
tor_norm_intens <- as.data.frame(tor_norm_intens)
tor_norm_intens$condition <- gsub("_\\S","",rownames(tor_norm_intens))
tor_norm_intens$condition <- factor(tor_norm_intens$condition,levels = c("BSA", "ISM1", "INS"))

# Loop that makes bar graphs of significant phosphopeptides for genes present 
# in GO terms with "muscle" included
for(i in 1:(ncol(tor_norm_intens)-1)){
  p <- ggplot(tor_norm_intens, aes(x=condition,y=tor_norm_intens[,i],fill=condition,width=.8)) + 
    geom_bar(position = 'dodge', stat = 'summary') +
    scale_fill_manual(values = brewer.pal(9, "Set1")[c(1,3,2)]) +
    geom_dotplot(binaxis='y', stackdir='center',binwidth = 0.05) +
    ylab("Fraction of maximum normalized intensity") +
    ggtitle(phospos[i])+
    theme_classic() + 
    #scale_y_continuous(expand = c(0, 0))+
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 45, 
                                     vjust = 1, hjust=1))
  ggsave(p,file=paste0("review/plots/bargraphs/filtered/tor/ism_vs_ins/",paste0(paste0(colnames(tor_norm_intens)[i]),"_bar_percentmax"),".jpeg"),
         units="cm", width=8, height=10, dpi=600)
  ggsave(p,file=paste0("review/plots/bargraphs/filtered/tor/ism_vs_ins/",paste0(paste0(colnames(tor_norm_intens)[i]),"_bar_percentmax"),".svg"),
         width=8, height=10)
}

# Making heatmap of phosphopeptides belonging to the
# TOR complex pathway that are present in at least one 
# ISM1 and BSA replicate, and where there is a significant difference 
# between the two conditions
#tor_phospep <- peptide_found %>% 
#  filter(BSA==TRUE & ISM1==TRUE) %>% 
#  filter(pathway=="TOR complex") %>%
#  pull(phospep)

tor_phospep <- peptide_found$phospep[peptide_found$pathway=="TOR complex"]
tor_norm_intens <- pp_norm@assays@data[[1]][rownames(pp_norm@assays@data[[1]]) %in% tor_phospep,]
tor_sig <- pp_norm@NAMES[pp_norm@elementMetadata$padj.vals.ISMvsBSA < 0.05 &  rownames(pp_norm@assays@data[[1]]) %in% rownames(tor_norm_intens)]
tor_norm_intens <- tor_norm_intens[rownames(tor_norm_intens) %in% tor_sig,]
tor_norm_intens <- tor_norm_intens-min_norm_val
# If missing set value to zero
tor_norm_intens[is.na(tor_norm_intens)] <- 0
tor_norm_intens <- apply(tor_norm_intens,1, function(x) {x/max(x)})
pos <- match(colnames(tor_norm_intens),pp_norm@elementMetadata$name)
# Tag phosphosite
phospos <- paste(gsub(";\\S+","",pp_norm@elementMetadata$Gene_names[pos]),
                 paste0(pp_norm@elementMetadata$Amino_acid[pos],
                        pp_norm@elementMetadata$Positions_within_proteins[pos]))
tor_norm_intens <- as.data.frame(tor_norm_intens)
tor_norm_intens$condition <- gsub("_\\S","",rownames(tor_norm_intens))
tor_norm_intens$condition <- factor(tor_norm_intens$condition,levels = c("BSA", "ISM1", "INS"))

# Loop that makes bar graphs of significant phosphopeptides for genes present 
# in the GO term TOR complex 
for(i in 1:(ncol(tor_norm_intens)-1)){
  p <- ggplot(tor_norm_intens, aes(x=condition,y=tor_norm_intens[,i],fill=condition,width=.8)) + 
    geom_bar(position = 'dodge', stat = 'summary') +
    scale_fill_manual(values = brewer.pal(9, "Set1")[c(1,3,2)]) +
    geom_dotplot(binaxis='y', stackdir='center',binwidth = 0.05) +
    ylab("Fraction of maximum normalized intensity") +
    ggtitle(phospos[i])+
    theme_classic() + 
    #scale_y_continuous(expand = c(0, 0))+
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 45, 
                                     vjust = 1, hjust=1))
  ggsave(p,file=paste0("review/plots/bargraphs/tor/ism_vs_bsa/",paste0(paste0(colnames(tor_norm_intens)[i]),"_bar_percentmax"),".jpeg"),
         units="cm", width=8, height=10, dpi=600)
  ggsave(p,file=paste0("review/plots/bargraphs/tor/ism_vs_bsa/",paste0(paste0(colnames(tor_norm_intens)[i]),"_bar_percentmax"),".svg"),
         width=8, height=10)
}

norm_intens <- as.data.frame(pp_norm@assays@data[[1]])
norm_intens <- norm_intens-min_norm_val
norm_intens[is.na(norm_intens)] <- 0
norm_intens <- t(apply(norm_intens,1, function(x) {x/max(x)}))

# Saving tables of significant differences between conditions
sig_table <- data.frame(proteins=pp_norm@elementMetadata$Proteins,
                        protein=pp_norm@elementMetadata$Protein,
                        prot_name=pp_norm@elementMetadata$Protein_names,
                        gene=pp_norm@elementMetadata$Gene_names,
                        aa=pp_norm@elementMetadata$Amino_acid,
                        pos=pp_norm@elementMetadata$Positions_within_proteins,
                        pval=pp_norm@elementMetadata$p.vals.ISMvsBSA,
                        padjval=pp_norm@elementMetadata$padj.vals.ISMvsBSA,
                        fc=pp_norm@elementMetadata$fc.vals.ISMvsBSA)
sig_table <- cbind(sig_table,norm_intens)
write.table(sig_table,row.names=FALSE,file=paste0("review/tables/",
  paste0("sig_ISMvsBSA_log_filter",".tsv")),sep="\t") 

sig_table <- data.frame(proteins=pp_norm@elementMetadata$Proteins,
                        protein=pp_norm@elementMetadata$Protein,
                        prot_name=pp_norm@elementMetadata$Protein_names,
                        gene=pp_norm@elementMetadata$Gene_names,
                        aa=pp_norm@elementMetadata$Amino_acid,
                        pos=pp_norm@elementMetadata$Positions_within_proteins,
                        pval=pp_norm@elementMetadata$p.vals.INSvsBSA,
                        padjval=pp_norm@elementMetadata$padj.vals.INSvsBSA,
                        fc=pp_norm@elementMetadata$fc.vals.INSvsBSA)
sig_table <- cbind(sig_table,norm_intens)
write.table(sig_table,row.names=FALSE,file=paste0("review/tables/",
  paste0("sig_INSvsBSA_log_filter",".tsv")),sep="\t")

sig_table <- data.frame(proteins=pp_norm@elementMetadata$Proteins,
                        protein=pp_norm@elementMetadata$Protein,
                        prot_name=pp_norm@elementMetadata$Protein_names,
                        gene=pp_norm@elementMetadata$Gene_names,
                        aa=pp_norm@elementMetadata$Amino_acid,
                        pos=pp_norm@elementMetadata$Positions_within_proteins,
                        pval=pp_norm@elementMetadata$p.vals.ISMvsINS,
                        padjval=pp_norm@elementMetadata$padj.vals.ISMvsINS,
                        fc=pp_norm@elementMetadata$fc.vals.ISMvsINS)
sig_table <- cbind(sig_table,norm_intens)
write.table(sig_table,row.names=FALSE,file=paste0("review/tables/",
  paste0("sig_ISMvsINS_log_filter",".tsv")),sep="\t") 
