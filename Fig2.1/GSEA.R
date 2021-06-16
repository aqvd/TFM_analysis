#######################################################
## 			Gene Set ENrichment Analysis
##  
##  Genes sorted by log2FC from DESesq2 analysis
##  and Gene sets from "GeneSets_forGSEA.R" passed to
##  fgesa() and then  barplots with top deregulated sets
##  are generated
#######################################################

library(tidyverse)
library(fgsea)

setwd("/Users/aqo/Desktop/MCF10A/RNA-seq")

## Log2FC sorted genes from DESeq2 results ENTREZID is also present as a column
genesL2FCsort_l <- plyr::llply(res, function(resTb){
	geneList <- resTb %>% as.data.frame() %>%
		rownames_to_column(var="SYMBOL") %>%
		dplyr::select(log2FoldChange,SYMBOL) %>%
		dplyr::filter(!is.na(log2FoldChange)) %>% 
		dplyr::arrange(desc(log2FoldChange)) %>% as_tibble() 
	return(geneList)
})

## load gene sets
load("/Users/aqo/Desktop/MCF10A/RNA-seq/geneSets-C2_H_CdLS_MCF10ArnaSeq@m_t2g.RData")

## List of gene sets name containing ENTREZIDs. 
msigdbr_list <- split(x = m_t2g$entrez_gene, f = m_t2g$gs_name)

## run fgsea
fgsea_l <- plyr::llply(genesL2FCannot_l, function(geneDF){
	## deframe sets the first selected variable as names of the vector
	ranks <- geneDF %>% dplyr::select(ENTREZID, log2FoldChange) %>% 
		dplyr::filter(!duplicated(ENTREZID)) %>% deframe()

	fg <- fgsea(pathways = msigdbr_list, 
		   stats = ranks,
		   eps = 0,
		   sampleSize = 101,
		   minSize = 10,
		   maxSize = Inf,
		   nPermSimple = 1000)

	return(fg)
})
fgsea_l


## Barplots of Top dereguñated Pathways
n <- names(fgsea_l)[1] 
for ( n in names(fgsea_l) ){
	ranks <- genesL2FCannot_l[[n]] %>%
			 dplyr::select(ENTREZID, log2FoldChange) %>% 
			 dplyr::filter(!duplicated(ENTREZID)) %>% deframe()

	fgseaRes <- fgsea_l[[n]]

	topPathwaysUp <- fgseaRes[ES > 0 & padj <= 0.05][head(order(padj), n=20),
		c("pathway","NES","ES", "padj") ]
	topPathwaysDown <- fgseaRes[ES < 0 & padj <= 0.05][head(order(padj), n=20),
		c("pathway","NES","ES", "padj") ]

	toGG <- rbind(topPathwaysUp, topPathwaysDown) %>% as_tibble() %>%
		dplyr::mutate_at("pathway", ~stringr::str_trunc(.,50,"right"))

	plotName <- paste0("barplotGSEA-H_C2CP_Shrink_", n, ".png")
	plotTitle <- paste0("Top 20 pathways GSEA\n", n)
	
	ggplot(toGG, aes(reorder(pathway, NES), NES)) +
		geom_col(aes(fill=as.factor(sign(NES)))) +
		theme_minimal() + 
		scale_fill_manual(values=c("dodgerblue3","firebrick2"),
			limits=c("-1","1")) +
		guides(fill="none") +
		coord_flip() +
		labs(x="", y="Normalized Enrichment Score",
			title=plotTitle) + 
		ggsave(plotName)
}
