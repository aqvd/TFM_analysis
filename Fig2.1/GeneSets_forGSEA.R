#######################################################
##     Create gene sets to use in "GSEA.R"
##
##    fgesea needs gene sets in a list object. 
##    This will be "m_t2g", named after "map_term2gene"
##
##    The names of the list elements are gene sets and the
##    elements are the ENTREZIDs of genes
##
##    

library(tidyverse)
options(dplyr.width = Inf)
options(dplyr.print_max = 30)

#######################################################
##      	Gene signatures from RNAseq analysis       
##
##   extract DEGs from DESeq2 results in variable "res"
#######################################################


setwd("/Users/aqo/Desktop/MCF10A/RNA-seq")

## Log2FC sorted names
genesL2FCsort_l <- plyr::llply(res, function(resTb){
	geneList <- resTb %>% as.data.frame() %>%
		rownames_to_column(var="SYMBOL") %>%
		dplyr::select(log2FoldChange,SYMBOL) %>%
		dplyr::filter(!is.na(log2FoldChange)) %>% 
		dplyr::arrange(desc(log2FoldChange)) %>% as_tibble() 
	return(geneList)
})
genesL2FCsort_l[[3]][3254:3299, ]
genesL2FCsort_l$siNIPBL_vs_WT2 %>% tail

## ========= Map symbols to ENSEMBL ========== ##
library(AnnotationDbi)

## Load Ensembl version for our assembly
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75 ##hg19 version of ensembl

# Return the Ensembl IDs for a set of genes
annotations_hg19 <- AnnotationDbi::select(edb, # database
                           keys = genesL2FCsort_l[[1]]$SYMBOL,  # data to use for retrieval
                           columns = c("GENEID", "ENTREZID"), #GENEID == ENSEMBL
                           keytype = "SYMBOL") # search 'keys' argument
dim(annotations_hg19)
annotations_hg19 <- as_tibble(annotations_hg19)

genesL2FCannot_l <- plyr::llply(genesL2FCsort_l, function(genes){
	genesAnnot <- dplyr::full_join(genes, annotations_hg19,
		by=c("SYMBOL"="SYMBOL")) %>%
		dplyr::filter(! is.na(GENEID)) %>%
		unique()
	return(genesAnnot)
})
genesL2FCannot_l


#######################################################
##   Gene sets from MsigDB retrieved using msigdbr
##  
##  Also add Hallmark (H), canonical pathways (CP)
##  and Chemical And Genetic Perturbations (CGP)
#######################################################

library(msigdbr)

## helper funcitons
msigdbr_species()
msigdbr_collections() %>% data.frame()

## retrieve speciffic set
cpg_gene_sets <- msigdbr(species = "Homo sapiens")
unique(cpg_gene_sets$gs_subcat)

## term 2 gene: df to map entrez to gene sets
categories <- c("H")
subcat <- c("CP", "CGP")

m_t2g <- cpg_gene_sets %>% 
	dplyr::filter(gs_cat %in% categories | gs_subcat %in% subcat ) %>%
	dplyr::select(gs_name, entrez_gene) %>%
	dplyr::mutate_at("gs_name", as.factor)
summary(m_t2g)

## >>>> Add our DEG genes to m_t2g <<<<<
names(resFilter)==names(genesL2FCannot_l)

## m_t2g has 2 columns: gs_name (endung in _DN or _UP) and entrez_gene
our_gs <- matrix(nrow=0, ncol=ncol(m_t2g)) %>%
	as.data.frame() %>% `colnames<-`(colnames(m_t2g))

for ( n in names(resFilter) ){
	## inner join with resFilter (FPKM, padj & l2fc filters) to get all DEG genes
	joined <- dplyr::inner_join(genesL2FCannot_l[[n]], resFilter[[n]],
		by=c("SYMBOL"="gene"))

	## divide into up or down deg
	up <- joined  %>% 
		dplyr::filter(log2FoldChange.y > 0) %>% 
		dplyr::mutate(gs_name=paste0("===>>>", n,"_UP")) %>%
		dplyr::rename(entrez_gene="ENTREZID") %>% 
		dplyr::select(gs_name,entrez_gene) %>%
		unique()
	our_gs <- rbind(our_gs, up)
	
	down <- joined %>% 
		dplyr::filter(log2FoldChange.y < 0) %>% 
		dplyr::mutate(gs_name=paste0("===>>>",n,"_DN")) %>%
		dplyr::rename(entrez_gene="ENTREZID") %>% 
		dplyr::select(gs_name,entrez_gene) %>%
		unique()

	our_gs <- rbind(our_gs, down)
}	
our_gs <- as_tibble(our_gs)

unique(our_gs$gs_name)

## rbind gith mapping term to gene (m_t2g)
m_t2g <- rbind(m_t2g, our_gs)
head(m_t2g);tail(m_t2g)
dim(m_t2g)

#######################################################
##         Cornelia de Lange signatures 
##
## Add signatures from CdLS from Liu et al. PLoS Biol (2009)
#######################################################

cdlSet_df <- read.table("/Users/aqo/Desktop/MCF10A/RNA-seq/CdLs_up_down_geneSets.txt", header = TRUE,
	sep="\t", stringsAsFactors = FALSE, fill=TRUE) 

cdlUP <- cdlSet_df %>% dplyr::select(CDLS_UP) %>%
	dplyr::inner_join(.,annotations_hg19, by=c("CDLS_UP"="SYMBOL")) %>%
	dplyr::select(ENTREZID) %>% mutate(gs_name="###CDLS_UP") %>%
	dplyr::rename(entrez_gene="ENTREZID")
m_t2g <- rbind(m_t2g, cdlUP)

cdlDOWN <- cdlSet_df %>% dplyr::select(CDLS_DOWN) %>%
	dplyr::inner_join(.,annotations_hg19, by=c("CDLS_DOWN"="SYMBOL")) %>%
	dplyr::select(ENTREZID) %>% mutate(gs_name="###CDLS_DOWN") %>%
	dplyr::rename(entrez_gene="ENTREZID")
m_t2g <- rbind(m_t2g, cdlDOWN)
dim(m_t2g)
tail(m_t2g)

## save gene sets for later use
save(m_t2g, file="geneSets-C2_H_CdLS_MCF10ArnaSeq@m_t2g.RData")