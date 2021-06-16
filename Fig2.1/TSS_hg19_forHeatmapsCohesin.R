##################################################
##  get TSS coordinates of DEGs and 5000 random
##  no DEGs to draw heatmaps of cohesin distribution
##  around TSS
##
##  DESeq results "res" and "resFilter" variables
##  are needed in order to know which genes to select 
##################################################

library(tidyverse)

wd <- "/Users/aqo/Desktop/cluster/MCF10A/RNA-seq/"

dir.create(file.path(wd, "DEGs_TSS_bed"), showWarnings = F)
setwd(file.path(wd, "DEGs_TSS_bed"))

##Function to write .bed from a list of GENEIDs contained in gtf
write_TSS_bed <- function(geneIDs, filename){
	AnnotationDbi::select(gtf, # database
		                  keys = geneIDs,  
						  columns = c("GENEID","TXCHROM","TXSTART","TXEND"), 
						  keytype = "GENEID") %>% 
	dplyr::group_by(GENEID) %>%
	dplyr::filter(row_number()==1) %>%
	dplyr::ungroup() %>%
	dplyr::transmute(`#chr`=TXCHROM, start=TXSTART, end=as.integer(TXSTART+1), id=GENEID) %>% 
	write.table(file=filename, sep = "\t", col.names = F, row.names = F,
		quote=F)
}


## GTF used for RNAseq gene counting
gtf <- GenomicFeatures:::makeTxDbFromGFF('/Users/aqo/Desktop/cluster/ewingCells/ewing_cellsRNAseq/DEGs_noParental/genes_hg19.gtf')

GenomeInfoDb::seqlevels(gtf) <- GenomeInfoDb::seqlevels(gtf)[1:24] # Only standard chromosomes

## UP GENES. resFilter contains only genes differentialy expressed
up <- resFilter$siNIPBL_vs_WT2 %>%
	dplyr::filter(log2FoldChange > 0) %>% 
	dplyr::pull(gene)
length(up) # 937 genes UP

write_TSS_bed(up, "siNIPBL_up_TSS.bed")

## DOWN genes
down <- resFilter$siNIPBL_vs_WT2 %>%
	dplyr::filter(log2FoldChange < 0) %>% 
	dplyr::pull(gene)
length(down) # 670

write_TSS_bed(down, "siNIPBL_down_TSS.bed")

## Sample 1000 genes no DEGs as a control in heatmaps
res$siNIPBL_vs_WT2 %>% data.frame() %>% 
	rownames() %>% .[! . %in% resFilter$siNIPBL_vs_WT2$gene] %>% 
	sample(., 5000, replace = FALSE) %>% 
	write_TSS_bed(., "siNIPBL_noDEG_TSS.bed")

