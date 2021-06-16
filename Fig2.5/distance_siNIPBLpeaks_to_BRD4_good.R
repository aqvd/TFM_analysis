###########################################################
####	  	  log2FC cohesin ~ distance BRD4           ####
##
##	CohesinSA2 is removed from BRD4 sites, so what happens
##  with its FC with respect to its distance to BRD4 peak?
##	
##	I have previously obtained log2FC from log2Ratio bw
##  using "BWsocres_inRegions_sameSample.py" so that each
##  protein gets their scores from their merged peaks
##
##  I have used the same strategy in "TEST_intersect_TADSandSCORES.R"
##  in order to plot FC arond the 6 states implified model.
##  So, starting from that point, i just have to join_nearest
##  with BRD4 peaks.
###########################################################

library(tidyverse)

options(dplyr.width=Inf)

## SET WORKING DIRECTORY
basedir <- "/Users/aqo/Desktop/cluster/MCF10A"
wd <- file.path(basedir, "gainedLost_peaks_distance_BRD4")
dir.create(wd, showWarnings = F)
setwd(wd)


## ====== MAKE GRANGES FROM BRD4 PEAKS SORTED BY INTENSITY AND 10groups =========
BRD4_peaks_file <- "/Users/aqo/Desktop/cluster/Other_Experiments_UCSC/BRD4_PRJNA295270/BRD4_merged_GSE72931_peaks_mergePeaks5000bp.bed"

BRD4_10groups_dir <- "/Users/aqo/Desktop/cluster/Other_Experiments_UCSC/BRD4_PRJNA295270/peaksBRD4_sorted"
BRD4_peaks_file_l <- list.files(BRD4_10groups_dir, pattern="10_groups.bed$") %>%
	setNames(., .)

## To plot deeptools sorted regions divided in 10 groups
BRD4_10groupsGR <- plyr::llply(BRD4_peaks_file_l, function(BRD4_peaks_file){
	pos <- gsub("BRD4.+DescBRD4_([^.]+)[.]bed.+", "\\1", BRD4_peaks_file )
	file <- file.path(BRD4_10groups_dir, BRD4_peaks_file)
	BRD4_peaks <- read.table(file, sep="\t", header=T) %>%
		as_tibble() %>%
		`colnames<-`(c("chr", "start", "end", "group")) %>%
		dplyr::mutate(toSEP = paste(chr, start, end, sep="_"),
			positions=pos) %>%
		GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)	
})

BRD4_peaks <- BRD4_10groupsGR[[1]]

## ====== GET LOG2FC SCORES FOR EACH PROTEIN "TEST_intersect_TADSandSCORES.R" ======== ##
scoresFile <- "/Users/aqo/Desktop/cluster/MCF10A/siNipbl-siSA_combined/log2ratio_siNipbl_siSA_ChIPs_newPseudocount/sorted_regions_log2desc/data_log2scores_siNIPBL_sameSample.tsv"

tadsDir <- "/Users/aqo/Desktop/cluster/MCF10A/TADS-MCF10A"
CHMM_states <- "/Users/aqo/Desktop/cluster/MCF10A/chromHMM-MCF10A/model-18states/MCF-10A_model18states_dense.bed"

#####################################################
####            SCORES, TADS, STATES             ####
#####################################################
tads_gr <- read.table(file.path(tadsDir,"TADborders_mcf10a_hg19_middle_and_length.sorted.bed"),
 		sep="\t", header=FALSE, stringsAsFactors = FALSE) %>% 
 		as_tibble() %>%
 		`colnames<-`(c("chr","start","end","tadID", "TADmiddle", "TADlength")) %>%
 		GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

states_gr <- read.table(CHMM_states, sep="\t", skip= 1, header= FALSE) %>% 
	dplyr::select(V1, V2, V3, V4) %>%
	`colnames<-`(c("chr", "start", "end", "state")) %>%
	GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

myIntersect <- function(data, states_gr) {
	data_gr <-  data %>% dplyr::select(chr, start, end, score) %>%
		dplyr::mutate(start= start + 50 , end= end - 50) %>%
		GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns=TRUE)

	res <- plyranges::join_overlap_inner(data_gr, states_gr)

	return(res)
}

# Read scores, pivot longer and group by regions (.bed) and experiment (.bw) to NEST
# then find Nearest tad borders and I have the scores to rank and compute correlations
nested <- read.delim(scoresFile, sep="\t", header=T) %>% 
	tidyr::separate(regions, sep="@", into= c("exp", "regions")) %>%
	dplyr::group_by(regions, exp) %>%
	tidyr::nest() %>%
	dplyr::mutate(expID=gsub("^(.+)_(.+)_(.*)bw_score$", "\\1_\\2", exp)) %>% 
	dplyr::mutate(intersected=purrr::map(data, 
						~myIntersect(., .GlobalEnv$states_gr)
						)
					) %>%
	dplyr::mutate(intersected=purrr::map(intersected,
					 	~plyranges::join_nearest(., .GlobalEnv$tads_gr) 
					)
				)

#Function to compute relative distance to TAD middle
distToMiddle_rel <- function(GR) {
	TADmid <- GenomicRanges::mcols(GR) %>% .$TADmiddle
	TADlen <- GenomicRanges::mcols(GR) %>% .$TADlength
	start <-  GenomicRanges::ranges(GR) %>% IRanges::start(.)

	dist_mid <- abs(TADmid - start)
	return(dist_mid / TADlen)

}

## Compute distance to TAD middle
nested <- nested %>% 
	dplyr::mutate( distTADmiddle=purrr::map(intersected, ~distToMiddle_rel(.)) )

## Compute distance to BRD4 peaks. join_nearest and parse toSEP col
nested <- nested %>% 
	dplyr::mutate( intersected=purrr::map(intersected, 
										  ~plyranges::join_nearest(., .GlobalEnv$BRD4_peaks) 
										  )
				  )
distBRD4_abs <- function(GR) {
	## get BRDcoordenates from toSEP column
	toSEP <- GenomicRanges::mcols(GR) %>% .$toSEP
	ch_st_nd <- stringr::str_split(toSEP, "_")
	midBRD4 <- lapply(ch_st_nd, function(x){
		( ( as.integer(x[2]) + as.integer(x[3]) ) / 2 )  %>% floor(.)
	}) %>% unlist()

	## get chip protein peak summit start coordinate
	start <-  GenomicRanges::ranges(GR) %>% IRanges::start(.)

	## absolute ditanve BRD4
	absBRD4=abs(midBRD4 - start)

	return(absBRD4)
}
## function to add amny vector as metadata of GRanges
add_mCol <- function(GR, vector, colName) {
	GenomicRanges::mcols(GR)[[colName]] <- vector

	return(GR)
}

nested <- nested %>% 
	dplyr::mutate( absBRD4=purrr::map(intersected, ~distBRD4_abs(.)) ) %>%
	dplyr::mutate( intersected= purrr::map2(intersected, absBRD4,
				 							~add_mCol(.x, .y, "absDist_BRD4")
				 							)
				 )


## ======= PREPARE FOR PLOTTING ======== ##
toName <- paste(dplyr::pull(nested, expID), dplyr::pull(nested, regions), sep = "_AT_" )
intersectedGR <- nested %>% dplyr::pull(intersected) %>% `names<-`(toName)

## Add simplified model anotations
annoDir <- "/Users/aqo/Desktop/cluster/MCF10A/siNipbl-siSA_combined/log2ratio_siNipbl_siSA_ChIPs_newPseudocount/sorted_regions_log2desc/10groups/CHMM_18states"

annotations <- read.table(file.path(annoDir,"state_annotations_simplified.txt"), 
	header=FALSE, sep="\t") %>% 
	`names<-`(c("state","anno","rgb")) %>%
	dplyr::mutate_at(c("anno"), as.factor) %>% 
	as_tibble()

cols <- colnames(intersectedGR[[1]] %>% data.frame())
cols <- c(cols, "exp")
toGG <- matrix(nrow=0, ncol= length(cols) ) %>% `colnames<-`(cols)

for (name in names(intersectedGR)) {
	GR <- intersectedGR[[name]]
	toGG <- GR %>% as.data.frame() %>% as_tibble() %>%
		dplyr::mutate(exp=name) %>% rbind(toGG, .)
}

toGG <- dplyr::left_join(toGG, annotations, by = c('state'="state"))
toGG <- toGG %>% tidyr::separate(exp, into=c("cond", "positions"),
	sep="_AT_") %>%
	dplyr::mutate(positions=gsub("siNipbl_merged_.+","peaks", positions),
		cond= gsub("_log2Ratio.*", "", cond)) %>%
	dplyr::mutate_at(c("cond", "positions"), as.factor) %>% 
	dplyr::filter(!absDist_BRD4 < 0, !is.na(absDist_BRD4), !is.infinite(absDist_BRD4)) %>%
	dplyr::mutate(cut_absDist_BRD4= 
						cut(absDist_BRD4, 
							breaks=c(seq(0,10000, length.out= 5), Inf),
							labels= c(1:5) 
							) 
				  )

toGG <- toGG %>% dplyr::filter(!is.infinite(absDist_BRD4))
summary(toGG)


## PLOT score ~ absDist_BRD4
toGG %>%
	dplyr::filter(anno %in% "insulator", group %in% 1:10) %>% 
	dplyr::group_by(cond) %>% tally %>% data.frame()


require(scales)

toGG %>% dplyr::group_by(cond) %>%
	dplyr::summarize(quantile( absDist_BRD4, 0:10/10) ) %>% 
	data.frame()


dummy_y <- toGG %>% dplyr::select(cond, score) %>%
	dplyr::mutate_at("cond", as.character) %>%
	dplyr::mutate(x= case_when(
								cond == "SA1_siNipbl" ~ -0.2,
								cond == "SA2_siNipbl" ~ -0.4,
								cond == "SMC1_siNipbl" ~ -.6
								) 
				  )
dummy_x <- toGG %>% dplyr::select(cond, absDist_BRD4) %>%
	dplyr::mutate_at("cond", as.character) %>%
	dplyr::mutate(y= case_when(
								cond == "SA1_siNipbl" ~ -3 -.05,
								cond == "SA2_siNipbl" ~ -3 -0.15,
								cond == "SMC1_siNipbl" ~ -3 -.25
								) 
				  )
myPlot <- function(dataGG, Condition, Groups, valueFill, Xlim, Ylim) {
	p <- ggplot(dataGG %>% dplyr::filter(group %in% Groups, cond %in% Condition), 
		aes(
				y = score,
			    x = absDist_BRD4 # natural log
	   		)
	   ) +
	theme_bw() +
	coord_cartesian(ylim=Ylim,
					xlim=Xlim
					) +
	geom_smooth( aes(color=cond) ) +
	facet_wrap(~cond, scales="free") +
	geom_rug( aes(color=cond), size=.1, alpha=.2, ) +
	scale_color_manual(breaks=Condition,
					   values=valueFill) + 
	scale_x_continuous(labels = comma) +
	guides(fill=FALSE) +
	labs(x="", y="") +
	theme(axis.text.x=element_text(size=rel(1.2), color="black", angle=45, hjust=1),
		axis.text.y=element_text(size=rel(1.2), color="black", angle=0, hjust=1),
		legend.position="none")

	return(p)
}

p1 <- myPlot(toGG, "SA1_siNipbl", 1:10, "firebrick3", c(0,39196), c(0,1))
p2 <- myPlot(toGG, "SA2_siNipbl", 1:10, "dodgerblue3", c(0,60425), c(-2.5,-1.5))
p3 <- myPlot(toGG, "SMC1_siNipbl", 1:10, "forestgreen", c(0,55301), c(-1.5,-0.5))

library(gridExtra)
pdf("logFC_distBRD4_linear_groups1-5_grid.pdf", width=8, height=3 )
	grid.arrange(p3, p2, p1, ncol =3)
dev.off()
 

## ========== OTEHR PLOTS =========== ##
ggplot(toGG %>% dplyr::filter(group %in% 1:5), 
		aes(
				y = score,
			    x = absDist_BRD4 # natural log
	   		)
	   ) +
	theme_bw() +
	coord_cartesian(ylim=c(-3, 1),
					xlim=c(200, 90000)
					) +
	geom_smooth( aes(color=cond) ) +
	facet_wrap(~cond, scales="free") +
	geom_rug( aes(color=cond), size=.1, alpha=.2, ) +
	scale_color_manual(breaks=c("SA1_siNipbl", "SA2_siNipbl", "SMC1_siNipbl"),
					   values=c("firebrick3","dodgerblue3","forestgreen")) + 
	scale_x_continuous(labels = comma) +
	theme(axis.text.x=element_text(size=rel(.9), color="black", angle=45, hjust=1))

ggsave("log2FC_distBRD4_log_allPeaks_STATES.pdf", height=7, width=7)

ggplot(toGG %>% dplyr::filter(group %in% 1:5), 
		aes(
				y = score,
			   x = absDist_BRD4 # natural log
	   		)
	   ) +
	theme_bw() +
	coord_cartesian(ylim=c(-3.5,3.5)) +
	geom_rug( aes(color=cond), size=.3, alpha=.4, ) +
	geom_smooth( aes(color=cond) ) +
	facet_wrap(~anno, scales="free") +
	scale_x_continuous(labels = comma) +
	theme(axis.text.x=element_text(size=rel(.9), color="black", angle=45, hjust=1))

ggsave("log2FC_distBRD4_linear_halfPeaks_STATES.pdf", height=7, width=7)

ggplot(toGG %>% dplyr::filter(group %in% 1), 
		aes(
				y = score,
			    x = absDist_BRD4 # natural log
	   		)
	   ) +
	theme_bw() +
	coord_cartesian(ylim=c(-3.5,3.5)) +
	geom_rug( aes(color=cond), size=.3, alpha=.4, ) +
	geom_smooth( aes(color=cond) ) +
	facet_wrap(~anno, scales="free") +
	scale_x_continuous(labels = comma) +
	theme(axis.text.x=element_text(size=rel(.9), color="black", angle=45, hjust=1))

ggsave("log2FC_distBRD4_linear_strongerPeaks_STATES.pdf", height=7, width=7)


##################################################

ggplot(toGG, 
		aes(
				y = absDist_BRD4,
			    x = cut(score, 10)
	   		)
	   ) +
	theme_classic() +
	geom_boxplot(
				aes( ),
				weight=20,
				outlier.size=.1,
				outlier.shape=12
				) +
	facet_wrap(~cond) +
	theme(axis.text.x=element_text(size=rel(.9), color="black", angle=45, hjust=1))