rm(list=ls())

######################
### Load libraries ###
######################

library(DESeq2)
library(easylabel)
library(ComplexHeatmap)
library(enrichR)
library(ggplot2)
library(reshape2)

#################
### Load data ###
#################

load(url("https://raw.githubusercontent.com/cankutcubuk/QMUL-msc-bioinfo/main/msc_peac_data.RData"))
# load("./msc_peac_data.RData")
ls()
str(syn_metadata)
str(rawcounts)
rawcounts[1:5,1:5]
syn_metadata[1:5,]
all(colnames(rawcounts)==syn_metadata$SampleID)
rawcounts <- rawcounts[,match(syn_metadata$SampleID, colnames(rawcounts))]
all(colnames(rawcounts)==syn_metadata$SampleID)

##################
### Run DESeq2 ###
##################

dds  <- DESeqDataSetFromMatrix(rawcounts, colData = syn_metadata, design=~Pathotype)
dds <- DESeq(dds) 

head(results(dds, contrast=c("Pathotype","Myeloid","Lymphoid")))
head(results(dds, contrast=c("Pathotype","Lymphoid","Myeloid"))) 

dssDF <- results(dds, contrast=c("Pathotype","Myeloid","Lymphoid")) 
dssDF <- as.data.frame(dssDF)
head(dssDF)

boxplot(split(log2(rawcounts["IGLC2",]+1), syn_metadata$Pathotype))
dssDF["IGLC2",]
sapply(split(rawcounts["IGLC2",], syn_metadata$Pathotype), sort)

boxplot(split(log2(rawcounts["MMRN1",]+1), syn_metadata$Pathotype))
dssDF["MMRN1",]
sapply(split(rawcounts["MMRN1",], syn_metadata$Pathotype), sort)

dssDF$Direction1 <- ifelse(dssDF$log2FoldChange > 0, "Up in Myeloid", "Up in Lymphoid")
dssDF$Direction2 <- ifelse(dssDF$log2FoldChange > 0, "Down in Lymphoid", "Down in Myeloid")

dssDF["IGLC2",]

####################
### Volcano plot ###
####################

easyVolcano(data=dssDF, x = "log2FoldChange", y = "pvalue", output_shiny = F, startLabels = c("IGLC2", "MMRN1"),
            fccut = 1, fdrcutoff = 0.05, colScheme=c('darkgrey', 'blue', 'lightblue', 'orange', 'red'), outline_col = NA)
easyVolcano(data=dssDF, x = "log2FoldChange", y = "pvalue", padj = "padj", output_shiny = F, startLabels = c("IGLC2", "MMRN1","CEMIP"),
            fccut = 1, fdrcutoff = 0.05, colScheme=c('darkgrey', 'blue', 'lightblue', 'orange', 'red'), outline_col = NA, ylim=c(0,25))

dssDF <- dssDF[order(dssDF$padj, decreasing = F),]
dssDF_p005 <- dssDF[which(dssDF$padj < 0.05),]
summary(dssDF_p005$log2FoldChange)

mygenes_up <- rownames(dssDF_p005)[which(dssDF_p005$log2FoldChange > 0)][1:20]
mygenes_down <- rownames(dssDF_p005)[which(dssDF_p005$log2FoldChange < 0)][1:20]
mygenes_all <- c(mygenes_down, mygenes_up)

###########################
### Count normalization ###
###########################

par(mfrow=c(3,1))
vstCounts <- DESeq2::vst(dds, blind=T)
vstCounts <- assay(vstCounts)
normCounts <- log2(counts(dds, normalized=TRUE) + 1)

lowexpressedGenes <- which(apply(counts(dds), 1, function(x) all(x<=100)))
boxplot(vstCounts[-lowexpressedGenes,], main="VST normalized")
boxplot(normCounts[-lowexpressedGenes,], main="Normalized by size factors")
boxplot(log2(counts(dds)[-lowexpressedGenes,] + 1), main="Log2(raw counts)", las=2)

df_dens_vst <- melt(vstCounts[-lowexpressedGenes,])
colnames(df_dens_vst) <- c("Gene","Sample","Expression")
p_vst <- ggplot(data = df_dens_vst, aes(x=Expression, color=Sample)) + geom_density(alpha=0.4) + theme_classic() + theme(legend.position = "none") + ylim(0, 0.3) + ggtitle("VST normalized")
p_vst

df_dens_norm <- melt(normCounts[-lowexpressedGenes,])
colnames(df_dens_norm) <- c("Gene","Sample","Expression")
p_norm <- ggplot(data = df_dens_norm, aes(x=Expression, color=Sample)) + geom_density(alpha=0.4) + theme_classic() + theme(legend.position = "none") + ylim(0, 0.3) + ggtitle("Normalized by size factors")
p_norm

df_dens_raw <- melt(log2(counts(dds)[-lowexpressedGenes,] + 1))
colnames(df_dens_raw) <- c("Gene","Sample","Expression")
p_raw <- ggplot(data = df_dens_raw,aes(x=Expression, color=Sample)) + geom_density(alpha=0.4) + theme_classic() + theme(legend.position = "none") + ylim(0, 0.3) + ggtitle("Log2(raw counts)")
p_raw

p_raw + p_vst + p_norm

####################
### Heatmap plot ###
####################

ScaledvstCounts <- t(scale(t(vstCounts[mygenes_all,])))

colBatch <- c("Batch1" = "#c7e9b6", "Batch2" = "#22318d")
colPatho <- c("Fibroid" = "#4bd950", "Myeloid" = "#d03435", "Lymphoid" = "#0004f3", "Ungraded" = "#bfc0bd")

column_ha = HeatmapAnnotation("Pathotype" = syn_metadata$Pathotype,
                              "Batch"=syn_metadata$Batch,
                              "CD20"=as.numeric(syn_metadata$CD20.max),
                              "CD138"=as.numeric(syn_metadata$CD68L.max),
                              col = list("Pathotype" = colPatho,
                                         "Batch"=colBatch))

set.seed(123)
HM_gene <- Heatmap(matrix = ScaledvstCounts, 
                     name = "z-score", 
                     cluster_rows = T, cluster_columns = T, 
                     column_km = 2, row_km = 3, 
                     col=c("blue","gray90","red"),
                     row_names_gp = gpar(fontsize = 8), 
                     top_annotation = column_ha
)

draw(HM_gene)

###################################
### Geneset enrichment analysis ###
###################################

listEnrichrSites()
setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
mydbs <- as.data.frame(dbs)
head(mydbs)
grep("human|sapiens|reactome|kegg",mydbs$libraryName, value = T, ignore.case = T)
selectedDatabases <- c("KEGG_2021_Human","Reactome_2022","WikiPathway_2023_Human")

enriched <- enrichr(genes = rownames(dssDF_p005), databases = selectedDatabases)
enriched$MergedDB <- rbind(enriched$KEGG_2021_Human, enriched$Reactome_2022, enriched$WikiPathway_2023_Human) # do.call("rbind", enriched[selectedDatabases])
enriched$MergedDB$GeneCount <- as.numeric(sapply(strsplit(enriched$MergedDB$Overlap,split = "/"), function(x) x[1]))
dim(enriched$MergedDB[enriched$MergedDB$P.value < 0.05,])
enriched$MergedDB <- enriched$MergedDB[enriched$MergedDB$P.value < 0.05,]
plotEnriched <- plotEnrich(enriched$MergedDB, showTerms = 34, numChar = 70, y = "Count", orderBy = "P.value")
plotEnriched <- plotEnriched + ggtitle("") + theme_bw(base_size = 14) + theme(text=element_text(color="black"),axis.text=element_text(color="black"))
plotEnriched

directionSelect <- "Up in Lymphoid"
enriched <- enrichr(genes = rownames(dssDF_p005)[dssDF_p005$Direction1 == directionSelect], databases = selectedDatabases)
enriched$MergedDB <- rbind(enriched$KEGG_2021_Human, enriched$Reactome_2022,  enriched$WikiPathway_2023_Human)
enriched$MergedDB$GeneCount <- as.numeric(sapply(strsplit(enriched$MergedDB$Overlap,split = "/"), function(x) x[1]))
dim(enriched$MergedDB[enriched$MergedDB$P.value < 0.05,])
enriched$MergedDB <- enriched$MergedDB[enriched$MergedDB$P.value < 0.05,]
plotEnriched_LM <- plotEnrich(enriched$MergedDB, showTerms = 26, numChar = 70, y = "Count", orderBy = "P.value")
plotEnriched_LM <- plotEnriched_LM + ggtitle(directionSelect) + theme_bw(base_size = 14) + theme(text=element_text(color="black"),axis.text=element_text(color="black"))
plotEnriched_LM

directionSelect <- "Up in Myeloid"
enriched <- enrichr(genes = rownames(dssDF_p005)[dssDF_p005$Direction1 == directionSelect], databases = selectedDatabases)
enriched$MergedDB <- rbind(enriched$KEGG_2021_Human, enriched$Reactome_2022,  enriched$WikiPathway_2023_Human)
enriched$MergedDB$GeneCount <- as.numeric(sapply(strsplit(enriched$MergedDB$Overlap,split = "/"), function(x) x[1]))
dim(enriched$MergedDB[enriched$MergedDB$P.value < 0.05,])
enriched$MergedDB <- enriched$MergedDB[enriched$MergedDB$P.value < 0.05,]
plotEnriched_ML <- plotEnrich(enriched$MergedDB, showTerms = 50, numChar = 70, y = "Count", orderBy = "P.value")
plotEnriched_ML <- plotEnriched_ML + ggtitle(directionSelect) + theme_bw(base_size = 14) + theme(text=element_text(color="black"),axis.text=element_text(color="black"))
plotEnriched_ML

##################
### PCA for QC ###
##################

dim(vstCounts)
vst.pca <- vstCounts[which(!apply(counts(dds), 1, function(x) all(x<=50))),]
dim(vst.pca)
pc <- prcomp(t(vst.pca), scale. = T)
varExp <- round(pc$sdev^2/sum(pc$sdev^2),digits = 2) * 100
pc_df <- as.data.frame(pc$x[,c("PC1","PC2")])
pc_df$Batch <- as.character(syn_metadata$Batch)
pc_df$Pathotype <- as.character(syn_metadata$Pathotype)

p <- ggplot(pc_df,aes(x=PC1,y=PC2, color=Pathotype))
p <- p+geom_point(size=3, alpha = 0.5) 
p <- p + theme_classic(base_size = 14) + guides(colour = guide_legend(override.aes = list(size=6)))
p <- p + xlab(paste0("PC1 variance %", varExp[1]))
p1 <- p + ylab(paste0("PC2 variance %", varExp[2]))

p <- ggplot(pc_df,aes(x=PC1,y=PC2, color=Batch))
p <- p+geom_point(size=3, alpha = 0.5) 
p <- p + theme_classic(base_size = 14) + guides(colour = guide_legend(override.aes = list(size=6)))
p <- p + xlab(paste0("PC1 variance %", varExp[1]))
p2 <- p + ylab(paste0("PC2 variance %", varExp[2]))

p1 + p2
table(syn_metadata$Pathotype, syn_metadata$Batch)

