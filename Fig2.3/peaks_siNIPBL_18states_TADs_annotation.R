###############################################################
##			Annotate peaks with States and TADs
##  
##  Script used to annotate called peaks from siNIPBL ChIP
##  with 18 states chromatin states and also annotate them
##  inside the TAD they lie in.
##
##  Then save an .RData object to load for different analysis
###############################################################

library(tidyverse)


## Peaks directory and dense.bed states
## Durectory where peaks pre condition .bed files are
PEAKSDIR <- "/Users/aqo/Desktop/cluster/MCF10A/siNipbl/Peaks/merged_replicates_callpeak/unique"

## Macs directory with summits.bed
MACSDIR <- "/Users/aqo/Desktop/cluster/MCF10A/siNipbl/Peaks/merged_replicates_callpeak"

##Chromain states model Dense .bed file
CHMMSTATES <- "/Users/aqo/Desktop/cluster/MCF10A/chromHMM-MCF10A/model-18states/MCF-10A_model18states_dense.bed"

## table to map state numbers to annotations
STATEANNO <- "/Users/aqo/Desktop/cluster/MCF10A/siNipbl/log2ratio_CHRMMstates/state_annotations_byHand.txt"

## tads hg19
TADS <- "/Users/aqo/Desktop/cluster/MCF10A/TADS-MCF10A/TAD_mcf10a_hg19_middle_and_length.sorted.bed"

## create relational table to organize files
configTab <- data.frame(PeakFile=list.files(PEAKSDIR, pattern="uniquePeaks.bed",
						full.names=T)) %>%
							mutate_if(is.factor, as.character)


configTab <- configTab %>%
	mutate(Sample=stringr::str_replace(basename(PeakFile),"_uniquePeaks.bed","")) %>%
	mutate(temp=Sample) %>%
	separate(temp, sep = "_", into = c("Prot", "Condition", "OnlyIn", "siRNA"))

## Use summits instead of peaks
summits <- list.files(MACSDIR, pattern="merged_summits.bed")
## Map each peak file with their corresponding summits
configTab$summits <- apply(configTab, 1, function(peakInfo){
	summitsProt <- summits[grep(peakInfo["Prot"], summits)]
	cond <- ifelse(peakInfo["Condition"]=="Treat", "siNipbl","siC")
	summitsToInter <- summitsProt[grep(cond, summitsProt)]
	return(summitsToInter)
})

configTab <- mutate(configTab, summits=file.path(MACSDIR, summits))

table(configTab$Prot, configTab$OnlyIn)
###############################################################
####		READ PEAKS AND CHROMATIN STATES  	 		   ####
###############################################################
peaksList <- apply(configTab, 1, function(peakInfo){
	peaks <- plyranges::read_bed(peakInfo[["PeakFile"]])
	summits <-  plyranges::read_bed(peakInfo[["summits"]])

	plyranges::join_overlap_inner(summits, peaks) %>%
		plyranges::select(name.y)
})

names(peaksList) <- configTab$Sample

## CHROMATIN STATES
states <- rtracklayer::import(CHMMSTATES, format="BED")
## =========== Add state annotations ============= ##
anno <- read.table(STATEANNO, 
	header=FALSE, sep="\t") %>% 
	`names<-`(c("state","anno", "rgb")) %>%
	dplyr::mutate_at("state", as.factor) %>% 
	dplyr::select(state, anno) %>%
	dplyr::mutate(anno=factor(anno, levels=anno)) %>%
	as_tibble()

statesAnno <- states %>% as.data.frame() %>%
		dplyr::left_join(.,anno, by=c("name"="state")) %>%
		GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

###############################################################
####				INTERSECT PEAKS AND STATES 	 		   ####
###############################################################
## Check multiple overlaps
cout_ol <- plyr::llply(peaksList, function(peaks){
	n_overlaps <- mutate(peaks, n_overlaps = plyranges::count_overlaps(peaks, statesAnno)) %>%
		as.data.frame() %>% pull(n_overlaps)
	sum(n_overlaps > 1)
})
## Intersect and join. If 1 peak intessect with more than 1 state,
## 2 or more tiems the peak returned
peakStates <- plyr::llply(peaksList, function(peaks){
	plyranges::join_overlap_inner(peaks, statesAnno) %>% 
		as.data.frame() %>%
		dplyr::rename(ID = "name.y", state = "name") %>% 
		dplyr::select(seqnames,start,end,width,strand,ID,itemRgb,anno) %>%
		dplyr::mutate(exp=gsub("_[0-9]+$","",ID)) %>%
		as_tibble()
})
summary(peakStates[[3]])
###############################################################
####				INTERSECT PEAKS AND TADS 	 		   ####
###############################################################
tads <- read.table(TADS, sep="\t", header=FALSE, stringsAsFactors = FALSE, 
			col.names = c("chr", "start", "end", "tadID", "middle", "length"))
tadRanges <- GenomicRanges::makeGRangesFromDataFrame(tads, 
				keep.extra.columns = TRUE) 

## Add TADS to states in peak summits
peakStatesTAD <- plyr::llply(peakStates, function(peaks){
	peaks <- GenomicRanges::makeGRangesFromDataFrame(peaks, 
				keep.extra.columns = TRUE) 
	plyranges::join_overlap_inner(peaks, tadRanges) %>%
		as.data.frame() %>%	
		dplyr::select(seqnames,start,end,width,strand,ID,tadID,itemRgb,anno) %>%
		dplyr::mutate(exp=gsub("_[0-9]+$","",ID)) %>%
		as_tibble() %>%
		dplyr::inner_join(.,tads, by=c("tadID"="tadID")) %>%
		dplyr::rename(start="start.x",end="end.x",
			startTAD="start.y", endTAD="end.y",
			middleTAD="middle", lengthTAD="length") %>%
		dplyr::select(-c(chr, width))

})

saveDir <- "/Users/aqo/Desktop/cluster/MCF10A/siNipbl/Peaks/merged_replicates_callpeak/unique"
savePath <- file.path(saveDir, "uniquePeaks_siNIPBL_intersected_18states_TADS@peakStatesTAD.RData")