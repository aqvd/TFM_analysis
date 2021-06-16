#######################################################
####             DEG detection
##
##   From the RNAseq unnormalized counts 
##   
#######################################################


setwd("/Users/aqo/Desktop/cluster/MCF10A/RNA-seq")

library(tidyverse)

##################################################
####              COUNT MATRIX                ####
##################################################
counts <- read.table("/Users/aqo/Desktop/cluster/MCF10A/RNA-seq/data/rna_seq_counts_all_samples.txt",
 sep="\t", header = TRUE, row.names = 1) %>% as.matrix()
counts%>%head()

condition <- colnames(counts) %>% stringr::str_replace("_[0-9]+","") %>%
  as.factor() 

condition 
sampleData <- data.frame(condition=condition)
rownames(sampleData) <- colnames(counts)
sampleData$exp <- c(rep("exp2", 9), rep("exp1",11)) %>% as.factor()
all(rownames(sampleData) == colnames(counts))
sampleData


##################################################
####              DESeq2 Annalysis            ####
##################################################
library(DESeq2)

## =============== DESeq NORMALIZATION ================= ##
## Construct a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sampleData,
                              design = ~ condition)
dds

## Add annotation data
featureData <- data.frame(gene=rownames(counts))
mcols(dds) <- data.frame(mcols(dds), featureData)
mcols(dds)
rowRanges(dds) <- rRangesDDS
## Set WT as reference for comparisons (altough contrasts argument could do
## more accurately)
dds$condition <- relevel(dds$condition, ref = "WT1")

## 1-Run analysis
dds <- DESeq(dds)
resultsNames(dds)

## =============== FPKM NORMALIZATION ================= ##
bedGenes <-  AnnotationDbi::select(edb, # database
                keys = rownames(dds),  # data to use for retrieval
        columns = c("GENEID", "EXONSEQSTART", "EXONSEQEND"), #GENEID == ENSEMBL
        keytype = "SYMBOL") %>% # search 'keys' argument 
      dplyr::filter(SYMBOL %in% rownames(dds))

table(bedGenes$SYMBOL %in% rownames(dds))

## Fill in rowRanges SLOW!!! Load dds pre annotated
load("/Users/aqo/Desktop/cluster/MCF10A/RNA-seq/dds-RNAsMCF10-A.RData")
rRangesDDS <- rowRanges(dds)
rRangesDDS
# for ( g in rownames(dds) ) {
#   if ( ! any(bedGenes$SYMBOL == g)) {
#     rowRanges(dds)[g] <- GRangesList(GRanges("chr1", IRanges(0,1000)))
#   } else {
#     start<-bedGenes$EXONSEQSTART[bedGenes$SYMBOL == g] %>% as.numeric()
#     end<-bedGenes$EXONSEQEND[bedGenes$SYMBOL == g] %>% as.numeric()
#     rowRanges(dds)[g] <- GRangesList(GRanges("chr1", IRanges(start,end)))
#   }
# }

fpkm_counts <- fpkm(dds) %>% as_tibble() %>% 
  mutate(gene=mcols(dds)$gene); fpkm_counts; fpkm_counts$gene %>% head()

####################################################
#### PREANNALYSIS: QC, PCA, Heatmap, clustering #### 
####################################################

### PCA ###
### Transform counts for data visualization via regulariz log transformation
rld <- rlog(dds, blind=TRUE)

rld_mat <- assay(rld) ## assay() is function from the "SummarizedExperiment" 
pca <- prcomp(t(rld_mat), scale=TRUE)

## Calculate proportion of variance
eigs <- pca$sdev^2 
percentVar <- round((eigs/sum(eigs)) * 100, 2)

# Create data frame with metadata and PC3 and PC4 values for input to ggplot
pca_df <- cbind(pca$x, sampleData) %>% as_tibble()

## select which components we wish to plot, change aes to be the same too!!
comp <- c(1,2)
ggplot(pca_df, aes(x=PC1, y=PC2, color=condition)) + 
  geom_point(size=2) +
  xlab(paste0("PC1: ", percentVar[comp[1]], "% variance")) +
    ylab(paste0("PC2: ", percentVar[comp[2]], "% variance")) +
    coord_fixed() +
    ggtitle("PCA with VST data") +
      theme(plot.title= element_text(size=rel(1.3), hjust=0.5))
  


### CORRELATION HEATMAP ###
library(pheatmap)
### Extract the rlog matrix from the object
rld_mat <- assay(rld)    
### Compute pairwise correlation values
rld_cor <- cor(rld_mat, method="pearson") ## cor() is a base R function
pheatmap(rld_cor, 
    annotation = sampleData, fontsize = 10, 
      fontsize_row = 10, height=20, 
      main="Pearson Correlation",
      show_colnames = FALSE)

### DISTANCE HEATMAP ###
rld_dist <- dist(t(rld_mat), method='manhattan')
rld_dist_mat <- as.matrix( rld_dist )
rownames(rld_dist_mat) <- colnames(counts) %>% 
              stringr::str_replace("_[0-9]+","")
colnames(rld_dist_mat) <- NULL
## check the output of cor(), make note of the rownames and colnames
head(rld_cor); head(rld_dist_mat)   

### Plot heatmap
pheatmap(rld_dist_mat,
    clustering_distance_rows = rld_dist,
    clustering_distance_cols = rld_dist,
    fontsize = 10, 
      fontsize_row = 10, height=15, 
      main="Manhattan Distance")


#######################################################
######                      DEG                  ######
#######################################################

## DF with all posible contrasts
df_contrasts <- dds$condition %>% as.character() %>% 
  grep("WT.", x=., invert = TRUE) %>%
  dds$condition[.] %>% as.character() %>% 
  expand.grid(cond="condition",treat= ., cont=c("WT1","WT2")) %>% 
  unique() %>% data.frame() %>% 
  mutate(temp=paste(treat, cont, sep="_vs_")) %>%
  tibble::column_to_rownames(var="temp"); df_contrasts

## Run results for every contrast
res <- apply(df_contrasts, 1, FUN = function(x) {
  results(dds,contrast=c(x[1],x[2],x[3]), alpha = 0.05)
})
names(res) 
## Select proper contrasts. WT1 and WT2 are controls of different experiment
contrasts <- names(res)[c(5,3,4,6,7)] ## CHANGE IT if new contrasts !!!!!!!!!!!!!!!!!!!!!!!
res <- res[contrasts]; names(res)

####################################
####  LFCShrink WT1 as reference ####
#####################################

dds$condition <- relevel(dds$condition, ref = "WT1")
## 1-Run analysis
dds <- DESeq(dds)
resultsNames(dds)

## =========== lfcShrink ========== ##
resWT1 <- rownames(df_contrasts) %>% 
  grep("WT1", .) %>% rownames(df_contrasts)[.]
resWT1 <- resWT1[grep("(siSA1)|(siSA2)|(siCTCF)", resWT1)]

resWT1shrink <- sapply(resWT1, FUN = function(x) {
  coeff <- paste0("condition_", x)
  lfcShrink(dds=dds, res=res[[x]], coef = coeff, type = 'apeglm') 
  })
## Summaryze number of up or downregulated genes
sapply(resWT1shrink, function(x) {
  return(summary(x))
})

#####################################
####  LFCShrink WT2 as reference ####
#####################################
dds$condition <- relevel(dds$condition, ref = "WT2")
## 1-Run analysis
dds <- DESeq(dds)
resultsNames(dds)

## =========== lfcShrink ========== ##
resWT2 <- rownames(df_contrasts) %>% 
  grep("WT2", .) %>% rownames(df_contrasts)[.]
resWT2 <- resWT2[grep("(siNIPBL)|(siWAPL)", resWT2)]

resWT2shrink <- sapply(resWT2, FUN = function(x) {
  coeff <- paste0("condition_", x)
  lfcShrink(dds=dds, res=res[[x]], coef = coeff, type = 'apeglm') 
  })
## Summaryze number of up or downregulated genes
sapply(resWT2shrink, function(x) {
  return(summary(x))
})

res <- c(resWT1shrink, resWT2shrink)

## MA plots after shrinking log2FC
par(mfrow=c(2,3))
for(i in names(res)){
  DESeq2::plotMA(res[[i]], main=i)
}

## ========== FPKM GREATER THAN 3  ========== ##
## fpkm > 3 in any of the 3 contrasts!
## if mean(FPKM) > 3 in any experiment, get that gene (row index)
ix_fpkm3 <- sapply(levels(condition), function(x){
  ixc <- grep(x, colnames(fpkm_counts))
  ixr <- which(rowMeans(fpkm_counts[,ixc]) > 3)
  return(ixr)
})
class(ix_fpkm3)
## Concatenate the list of rows and get rid of duplicates
ix_genesFPKM_3 <- do.call("c", ix_fpkm3) %>% .[!duplicated(.)]
length(ix_genesFPKM_3) ## from 14053 to 9570 genes remaining

## ===================== FILTER DEG ======================== ##
## ====== padj <= 0.05 & abs(LFC) >= log2(1.5) & FPKM > 3 FILTER ====== ##

## Set thresholds
padj_cutoff <- 0.05
# 1 is no change, 1.5 is 50% (up or under)
abs_fold_change_percent <-  1.5 
(lfc_cutoff <- log2(abs_fold_change_percent))

resFilter <- plyr::llply(res, function(x) {
  x %>% data.frame() %>%
  .[ix_genesFPKM_3, ] %>%
  dplyr::filter(padj <= padj_cutoff & abs(log2FoldChange) >= lfc_cutoff) %>% 
  rownames_to_column(var="gene") %>%
  arrange(padj) %>% as_tibble()

})
resFilter
plyr::llply(resFilter,function(x){
  t <- as.data.frame(x) %>% pull(log2FoldChange)
  return(table(t>0))
})

## List of all DEG genes between conditions
degGenes <- unlist(sapply(resFilter, FUN = function(x) {
  x$gene
})) %>% unname() %>% unique()
length(degGenes) ## 2265 unique DEG genes in all exps


#######################################################################
####                     EULER PLOT DEGs                           ####
#######################################################################
resFilter

DEGs <- vector("list", length = length(resFilter)) %>%
  `names<-`(names(resFilter))
UPs <- vector("list", length = length(resFilter)) %>%
  `names<-`(names(resFilter))
DOWNs <- vector("list", length = length(resFilter)) %>%
  `names<-`(names(resFilter))

for(name in names(resFilter)){
  res <- resFilter[[name]]
  DEGs[[name]] <- pull(res,gene)
  UPs[[name]] <- res %>% dplyr::filter(log2FoldChange > 0) %>% pull(gene)
  DOWNs[[name]] <- res %>% dplyr::filter(log2FoldChange < 0) %>% pull(gene)
}
cols <- c(siSA1_vs_WT1="firebrick3", siSA2_vs_WT1="dodgerblue3", siCTCF_vs_WT1="goldenrod1",
	siNIPBL_vs_WT2="forestgreen", siWAPL_vs_WT2="ghostwhite")

labels <- gsub("_"," ",names(resFilter)) %>% `names<-`(names(resFilter))

ix <- c("siSA2_vs_WT1","siSA1_vs_WT1", "siNIPBL_vs_WT2")

pdf("DEGs_SA1_SA2_NIPBL.pdf", width = 8.83*0.5, height = 9.19*0.5)
plot(eulerr::euler(DEGs[ix],
  shape='ellipse'),
  quantities=list(
  	cex=1,
  	type="counts"), 
  fill=list(
  	fill=cols[ix],
  	alpha=0.8), 
  edges=FALSE,
  labels=list(
  	labels="",
  	font=2,
  	col="black",
  	cex=1.1), 
  legend=list(
  	labels=labels[ix],
  	cex=1.1), 
  main=list(
  	label="DEGs",
  	font=4,
  	cex=1.6))
dev.off()

 
##################################################
####              ANNOTATION                  ####
##################################################
library(AnnotationDbi)

## /* NOT RUN: To check version of ensembl for hg19/GRCh37
# library(AnnotationHub) 
# ah <- AnnotationHub()
# query(ah, c("GRCh37", "ensembl")) 
##*/ NOT RUN

## Load Ensembl version for our assembly
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
columns(edb) # possible columns to select

# Return the Ensembl IDs for a set of genes
annotations_hg19 <- AnnotationDbi::select(edb, # database
                           keys = rownames(dds),  # data to use for retrieval
                           columns = c("GENEID", "ENTREZID"), #GENEID == ENSEMBL
                           keytype = "SYMBOL") # search 'keys' argument
dim(annotations_hg19)

## Remove NA values for GENEID (ensemblID) or ENTREZID, select unique simbols
annotations_hg19 <- annotations_hg19 %>%  
  dplyr::filter(!is.na(GENEID) & !is.na(ENTREZID)) %>%
  dplyr::distinct()

dim(annotations_hg19)
head(annotations_hg19)

## Left Join normalized counts and annotations
normalized_counts <- counts(dds, normalized=T) %>%
   data.frame() %>%
   rownames_to_column(var="gene") %>%
   as_tibble() %>%
   dplyr::left_join(annotations_hg19, by=c("gene" = "SYMBOL")) 


##################################################
####                    SVA                   ####
##################################################
library(sva)

## Our data are DSeq2 normalized counts
colData(dds)$exp <- as.factor(sampleData$exp)

## Build model matrices
mod = model.matrix( ~ exp, data=colData(dds))
mod0 = model.matrix(~ 1, data=colData(dds))
countsBat <- sva::ComBat(counts(dds, normalized = TRUE),
   batch=colData(dds)$exp)

##################################################
####                 HEATMAP                  ####
##################################################
library(pheatmap)
library(RColorBrewer)

## ================= ALL Differentialy Expresed Genes =============== ##
## Extract normalized expresion of significant genes. Column 1 is gene
## Duplicated genes are summarized by mean
datHeatCombat <- countsBat %>% data.frame() %>% 
  rownames_to_column(var="gene") %>%
  dplyr::filter(gene %in% degGenes) %>%

  dplyr::group_by(gene) %>%
  dplyr::summarise_all(mean) %>%
  column_to_rownames(var="gene") %>% 
  as.matrix()

head(datHeatCombat); dim(datHeatCombat)

### Set a color palette
heat_colors <- colorRampPalette(brewer.pal(11, "RdBu"))(255)

### Run pheatmap using the metadata data frame for the annotation
heat_annCol <- list(
  condition=c(siNIPBL="forestgreen",siSA1="firebrick3",siSA2="dodgerblue3",
    siWAPL="ghostwhite",WT1="grey20",WT2="grey50", siCTCF="goldenrod1"),
  exp = c(exp1="black",exp2="orange"))

heatName <- "~/Desktop/MCF10A/RNA-seq/heatmap_allDEG_ComBat_annoColors.tiff"
pheatmap(datHeatCombat, 
    color = rev(heat_colors), 
    cluster_rows = TRUE,
    show_rownames = FALSE,
    scale = "row", 
    annotation_col = sampleData,
    annotation_colors = heat_annCol,
    border_color = NA, 
    main = "All DEGs (No FPKM fitered)",
    fontsize = 10, 
    fontsize_col = 8,
    fontsize_row = 8,
    legend=TRUE,
    legend.cex=0.1, 
    height = 10,
    width= 9, 
    filename = heatName)

