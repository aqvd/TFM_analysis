#######################################################
##
##   Volcano plot of genes annalyzed in 
##   "DEG_analysys_DESeq2_heatmap.R"
##
#######################################################


library(tidyverse)


##################################################
####        CLASIFY genes as DownNotUp        ####
##################################################

# ======== Annotate genes as Down/Up/Unchanged ========== #
DownNotUp <- function(resTable, fc) {
  if ( resTable["padj"] > 0.05 || is.na(resTable["padj"])) {
    return("NoDEG")
  } else {
    if ( abs(resTable["log2FoldChange"] ) < log2(fc) ) {
      return("LowFC") 
    } else if ( resTable["log2FoldChange"] <= -log2(fc) ) {
      return("Down")
    } else {
      return("Up")
    }
  }
}

## Add column with annotated change to res results list. 
## Better add at the end so that apply does not coerce row into character
resChange <- plyr::llply(res, function(resTb){  
   change <- apply(resTb, 1, function(row) {
                    DownNotUp(row, 1.5)
   })
   resTb <- resTb %>% data.frame() %>% rownames_to_column(var="gene") %>%
   dplyr::mutate(change=change)
   resTb$FPKMgt3[ix_genesFPKM_3] <- "yes"
   resTb$FPKMgt3[is.na(resTb$FPKMgt3)] <- "no" 
   resTb <-resTb %>% mutate_at(c("change","FPKMgt3"), as.factor)
   return(as_tibble(resTb))
})

resChange[[1]] %>% head # check the generated tibbles
resChange[[1]] %>% summary()

##################################################
####                 VOLCANO                  ####
##################################################
library(ggrepel)

## Iterate through every contrast. resChange is a list with DESeq2 results and the type of
## change annotated. It is not filtered by padj nor log2Fc
cont <- "siNIPBL_vs_WT2"
for (cont in names(resChange)) {
  resTb <- resChange[[cont]]

  head(resTb)
  ## colum genelabels contains symbols of DEG genes according to our filtering
  ## criteria. Also remove padj == 0
  toGG <- resTb %>%  dplyr::select(gene, padj, log2FoldChange) %>% 
    dplyr::filter(padj <= 0.05 & abs(log2FoldChange) >= lfc_cutoff) %>% 
    dplyr::mutate(genelabels= gene) %>%
    dplyr::select(genelabels, gene) %>%
    dplyr::right_join(resTb, by=c("gene"="gene")) %>% 
    dplyr::select(gene, genelabels, log2FoldChange, padj, change, FPKMgt3) %>%
    dplyr::filter(FPKMgt3 == "yes") %>%
    dplyr::arrange(padj) %>%
    dplyr::mutate(padjVolc=ifelse(padj==0, 1e-300, padj))

  head(toGG)
  tail(toGG)
  ## DEG genes have its SYMBOl in genelabels, No DEG have NAs
  n_DEG <- sum(!is.na(toGG$genelabels))
  n_DEG <- paste0(" Number of DEG genes: ", n_DEG)
  ## Keep 20 most significant and interest genes manualy introduced
  ix <- which(toGG$gene %in% c("NIPBL", "CTCF", "STAG1", "STAG2", "WAPAL", "PDS5", "ESCO1"))
  toGG$genelabels[20:nrow(toGG)] <- ""
  toGG$genelabels[ix] <- toGG[ix, ]$gene

  ggplot(toGG, aes(x = log2FoldChange, y = -log10(padjVolc))) +
    theme_classic() +
      geom_point(aes(colour = change), size=1) +
      scale_color_manual(values=c('Down'='deepskyblue2', 'LowFC'='gray50',  
                    'NoDEG'='gray90','Up' = 'firebrick'))+
      geom_text_repel(aes(
        label = genelabels, segment.curvature= -.1, segment.angle=30 ,
        segment.ncp= 3, segment.shape=-.5, segment.square= TRUE, segment.inflect= TRUE, 
        ),
       arrow = arrow(length = unit(0.015, "npc")),
       max.overlaps= Inf, max.time= 1, size = 2.5,  min.segment.length=3,
        
        ) +
      scale_y_continuous(expand = expansion(mult = 0.1)) +
      ggtitle(paste(cont,"FPKM > 3 filter.\n" ,n_DEG, sep=" ")) +
      xlab("log2 fold change") + 
      ylab("-log10 adjusted p-value") +
      theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
            axis.title = element_text(size = rel(1.25))) 
  filename <- paste0("~/Desktop/cluster/MCF10A/RNA-seq/volcano_",cont,"_FPKMgt3.png" )
  ggsave(filename, width=7, height=7)
}

#######################################################
######                SAVE RESULTS               ######
#######################################################
save(resChange, 
  file="/Users/aqo/Desktop/MCF10A/RNA-seq/DESeqResults@resChange.RData")
 