#####################################################
##    ODDS RATIO DEGs in Super Enhancers
##
## To measure if there is an association between 
## genes in Super Enhancers and being
## deregulated by NIPBL, SA1 or SA2 depletion 
library(tidyverse)

wd <- "/Users/aqo/Desktop/cluster/MCF10A/RNA-seq"
setwd(wd)

## DESeq results to classify genes in UP, NoDEG and DOPWN
load("/Users/aqo/Desktop/cluster/MCF10A/RNA-seq/DESeqResults@resChange.RData")
resChange
load("/Users/aqo/Desktop/cluster/MCF10A/RNA-seq/dds-RNAsMCF10-A.RData")
dds$condition

##################################################
####       SUPER ENHANCER ANNOTAION         ####
##################################################
## Read enhancer table
SE <- read.table("/Users/aqo/Desktop/cluster/MCF10A/enhancers-MCF10A/super-enhancers_MCF10A-TableS1-Zhang2019-BreastCancerRes.bed",
  header=FALSE, sep = "\t") %>% as_tibble() %>%
  `colnames<-`(c("chr", "start", "end"))
dim(SE) ## 343 super enhancers
median(SE$end - SE$start)
SE_ranges <-GenomicRanges::makeGRangesFromDataFrame(SE, 
	keep.extra.columns = TRUE)

SE_hiC <- read.table("/Users/aqo/Desktop/cluster/MCF10A/enhancers-MCF10A/super_genes_5000.bed",
  header=FALSE, sep = "\t") %>% as_tibble() %>%
  `colnames<-`(c("chr", "start", "end", "Hi-C", "Hi-C_chr","V6","V7","GeneInteraction"))
dim(SE_hiC)
length(unique(SE_hiC$GeneInteraction)) ## 93 unique enh:hene interaction

SE_hiCRanges <-GenomicRanges::makeGRangesFromDataFrame(SE_hiC, 
	keep.extra.columns = TRUE)

#####################################################
####  LOAD ANNOTATION hg19 to map genes into SE  ####
#####################################################
library(AnnotationDbi)

## Load Ensembl version for our assembly
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
columns(edb) # possible columns to select

# Return the Ensembl IDs for a set of genes
annotations_hg19 <- AnnotationDbi::select(edb, # database
                           keys = resChange[[1]]$gene,  # data to use for retrieval
                           columns = c("SEQNAME","GENESEQSTART","GENESEQEND",
                           	"GENEID", "ENTREZID"), #GENEID == ENSEMBL
                           keytype = "SYMBOL") # search 'keys' argument
dim(annotations_hg19)

## Remove NA values for GENEID (ensemblID) or ENTREZID, select unique simbols
annotations_hg19 <- annotations_hg19 %>%  
  dplyr::filter(!is.na(GENEID) & !is.na(ENTREZID)) %>%
  dplyr::distinct() %>% 
  dplyr::mutate(SEQNAME=paste0("chr",SEQNAME)) %>%
  dplyr::select(SEQNAME, GENESEQSTART, GENESEQEND, SYMBOL, GENEID, ENTREZID) %>%
  as_tibble()

names(annotations_hg19) <- c("chr", "start", "end", "symbol", "ensID", "entrezID")
dim(annotations_hg19)
head(annotations_hg19)

hg19 <- GenomicRanges::makeGRangesFromDataFrame(annotations_hg19, 
	keep.extra.columns = TRUE)

#####################################################
#### 		 MAP genes into Super Enhancers      ####
#####################################################
library(annotatr)

## >>>>>>>>> DEFINITION OF SUPER ENHANCER ANNOTATION: >>>>>>>>>>> ##
## Get all genes that lie inside the SE plus the genes 
## where a Hi-C significative interaction is detected
## <<<<<<<<< DEFINITION OF SUPER ENHANCER ANNOTATION: <<<<<<<<<<< ##

## INERSECT: some intersect once, others more than one time
GenomicRanges::countOverlaps(SE_ranges,hg19)
## intersect hg19 annotations with SE, pull unique symbols and 
## concatenate wuth SE_hiC$GeneInteraction
genesInSE <- IRanges::subsetByOverlaps(hg19,SE_ranges) %>% 
	data.frame() %>%
	pull(symbol) %>% unique() %>% 
	c(.,as.character(SE_hiC$GeneInteraction)) %>%
	unique() %>%
	tibble(SEgene=., Enhancer_ID=1:length(.))

## 360 uniquie genes annotated in super enhancers 
## combinig intersect (279) and Hi-C (93)
dim(genesInSE) 

## tabulate the number of genes that we can annotate
plyr::llply(resChange, function(resTb) {
  table(genesInSE$SEgene %in% resTb$gene)
})
# FALSE  TRUE 
#    57   303

## ## Do right joining genesInSE - resChange to keep all genes and some 
## annotated as with Enhancer ID
enhChange <- plyr::llply(resChange, function(resTb) {
  genesInSE %>% dplyr::select(SEgene, Enhancer_ID) %>% as_tibble() %>%
  	dplyr::right_join(resTb, by = c("SEgene"="gene"))
})
table(is.na(enhChange[[1]]$Enhancer_ID)) ## Check 303 are not NA

#################################################
####     ODDS RATIO: Up,Down,LowFC, NoDEG    ####
#################################################

## Count how many genes change for each results table for enhancer/NoEnhancer
tt <- plyr::llply(enhChange, function(resTb) {
  e <- resTb %>% dplyr::filter(!is.na(Enhancer_ID)) %>%
    group_by(change) %>% summarize(n=n()) 
  ne <- resTb %>% dplyr::filter(is.na(Enhancer_ID)) %>%
    group_by(change) %>% summarize(n=n())

  return(list(enhancer=e,noEnhancer= ne))
})
tt

contigTable <- function(v_chang, Tb, filterBy) {
  resL <- vector("list", length(v_chang)); names(resL) <- v_chang
  for (ch in v_chang) {
    ## Extract nº genes per change class for enhancer / NoEnhancer.
    ## Sum(.) to transfor cases where filter returns integer(0) to 0
    yes <- Tb$enhancer %>% dplyr::filter( .data[[filterBy]] == ch ) %>% pull(n) %>% sum(.)
    restY <- Tb$enhancer %>% dplyr::filter(.data[[filterBy]]!= ch ) %>% pull(n) %>% sum(.)
    no <- Tb$noEnhancer %>% dplyr::filter(.data[[filterBy]] == ch ) %>% pull(n) %>% sum(.)
    restN <- Tb$noEnhancer %>% dplyr::filter(.data[[filterBy]] != ch ) %>% pull(n) %>% sum(.)


    ## Create matrix
    dnms <- list(filterBy=c(ch, "-"), "Enhancer"=c("Yes", "No"))
    names(dnms)[1] <- filterBy 
    mat <- matrix(data =c(yes,restY,no,restN), byrow = FALSE, ncol = 2, 
            dimnames = dnms)
    resL[[ch]] <- mat
  }
  return(resL)
}
changes <- levels(enhChange[[1]]$change)
contrasts <- names(enhChange)

#############################################
#### 			     ODDS RATIO 	CI			     ####
#############################################

## Function for OR with CI plotting
plotOR <- function(contrasts, changes, contTbl) {
  or_df <-  matrix(1:(length(contrasts) * length(changes) * 6),
    ncol = 6) %>%   data.frame() %>%  
    `colnames<-`(c("contrast", "change", "OR", "lower", "upper","p.val")) %>%
    mutate(contrast=rep(contrasts, length(changes)), 
      change=rep(changes,length(contrasts))) 

  for (cont in names(resChange) ){
    for (ch in changes ) {
      or <-  exact2x2::exact2x2(contTbL[[cont]][[ch]], tsmethod = "minlike")
      ix <- with (or_df, change == ch & contrast == cont)
      or_df[ix, "OR"] <- or$estimate
      or_df[ix, "lower"] <- or$conf.int[1]
      or_df[ix, "upper"] <- or$conf.int[2]
      or_df[ix, "p.val"] <- or$p.value
    }
  }
  or_df <- dplyr::mutate(or_df, signif=ifelse(or_df$p.val < 0.05, "yes", "no")) %>%
    mutate_at(c("OR","upper","lower"), ~.+1e-7) %>%
    mutate_at(c("OR","upper","lower"), log2)

  p <- ggplot(data=or_df, mapping=aes(x=OR, y=change)) +
    geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") +
    facet_wrap(~contrast, scales="free") +
    geom_errorbarh(aes(xmax=upper, xmin=lower), size=0.5, 
      height = 0.2, color = "gray50")+ 
    geom_point(aes(x=OR, y=change, color=signif), shape=18, size=2, ) +
    geom_text(aes(x=OR, label=round(OR,2)), size=1.5,
          vjust=-0.9, position=position_dodge2(width=0.9)) +
    scale_color_manual(
      values=c("orange","grey79"), 
      limits=c("yes","no"),
      name="Significative") +
    scale_x_continuous( limits=c(-5,5),
                        breaks=seq(-4,4, by=2)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    ylab('') +
    xlab('Log2 Odds Ratio') +
    ggtitle('Odds ratio Super Enhancer Genes') +
    theme(plot.title=element_text(hjust=0.5),
      axis.text.x=element_text(
      angle=30, hjust=0.5, vjust= 0.5,
      color= "black"),
      axis.text.y=element_text(
        color= "black"),  
      panel.grid.major= element_blank(),
      strip.background = element_rect(
      fill="black"
      ),
      strip.text = element_text(
        face= "bold", color= "white"
        )
      )

  return(p)
}

## /* UP DOWN NoDEG /*
## after applying FPKM filter and creation of IsDEG column, split Up and Down
enhChange <- plyr::llply(enhChange, function(resTb) {
  resTb$change[resTb$FPKMgt3 == "no"] <- "NoDEG"
  resTb <- dplyr::mutate(resTb, IsDEG=factor(IsDEG, levels = c("Down", "Up", "NoDEG")))
  resTb$IsDEG[resTb$change=="Up"] <- "Up"
  resTb$IsDEG[resTb$change=="Down"] <- "Down"
  return(resTb)
})
summary(enhChange[[1]])
## */  UP DOWN NoDEG */

## Count number of genes per with enhancer and no enhancer for DEG NoDEG
tt <- plyr::llply(enhChange, function(resTb) {
  e <- resTb %>% dplyr::filter(!is.na(Enhancer_ID)) %>%
  group_by(IsDEG) %>% summarize(n=n()) 
  ne <- resTb %>% dplyr::filter(is.na(Enhancer_ID)) %>%
  group_by(IsDEG) %>% summarize(n=n())

  return(list(enhancer=e,noEnhancer= ne))
})

changes <- levels(enhChange[[1]]$IsDEG)
## Create contig tables
contTbL <- lapply(tt, function(tb) {
  contigTable(changes, tb, "IsDEG")
})

p <- plotOR(contrasts, changes, contTbl)
print(p)
ggsave("oddsRatioCI-UpNoDown-SuperEnhancers.png", height=4, width=6.5)
