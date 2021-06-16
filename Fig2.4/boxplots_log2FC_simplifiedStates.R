###########################################################################
##     log2FC in NIPBL KD in the groupped 6 states model
##
##  For Fig 2.4B
##
##  Using "BWsocres_inRegions_sameSample.py" get log2fold change in
##  summits +-50bp.
##
##  Then, intersect with chromatin states and TADs, as we keep the original
##  summits positions
##
##  Finally, draw boxplots
###########################################################################


library(tidyverse)

scoresFile <- "/Users/aqo/Desktop/cluster/MCF10A/siNipbl-siSA_combined/log2ratio_siNipbl_siSA_ChIPs_newPseudocount/sorted_regions_log2desc/data_log2scores_siNIPBL_sameSample.tsv"

tadsDir <- "/Users/aqo/Desktop/cluster/MCF10A/TADS-MCF10A"
CHMM_states <- "/Users/aqo/Desktop/cluster/MCF10A/chromHMM-MCF10A/model-18states/MCF-10A_model18states_dense.bed"

wd <- "/Users/aqo/Desktop/cluster/MCF10A/siNipbl-siSA_combined/log2ratio_siNipbl_siSA_ChIPs_newPseudocount/sorted_regions_log2desc/"
setwd(wd)

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

# Function to compute relative distance to TAD middle
distToMiddle_rel <- function(GR) {
	TADmid <- GenomicRanges::mcols(GR) %>% .$TADmiddle
	TADlen <- GenomicRanges::mcols(GR) %>% .$TADlength
	start <-  GenomicRanges::ranges(GR) %>% IRanges::start(.)

	dist_mid <- abs(TADmid - start)
	return(dist_mid / TADlen)

}

nested <- nested %>% 
	dplyr::mutate( distTADmiddle=purrr::map(intersected, ~distToMiddle_rel(.)) )

# Correlation bwtween socre and dist TADmiddle

corGRanges <- function(GR, column, corMethod ){
	score <-  GenomicRanges::mcols(GR) %>% .$score
	corr <- cor.test(column, score, method = corMethod)

	return(corr)
}


nested <- nested %>% 
	dplyr::mutate( spearmanTest=purrr::pmap(
						list(ranges=intersected, dist=distTADmiddle),
						function(ranges, dist){
								corGRanges(ranges, dist, "spearman")
							}
						)
					)
warnings()
nested <- nested %>% 
	dplyr::mutate(
		rho = purrr::map_dbl(spearmanTest, ~.$estimate),
		pVal = purrr::map_dbl(spearmanTest, ~.$p.value)
		)

#####################################################
####            MEAN SCORE PER STATE             ####
#####################################################
toName <- paste(dplyr::pull(nested, expID), dplyr::pull(nested, regions), sep = "_AT_" )
scoresTADs <- nested %>% dplyr::pull(intersected) %>% `names<-`(toName)

## =========== Add state annotations ============= ##
annoDir <- "/Users/aqo/Desktop/cluster/MCF10A/siNipbl-siSA_combined/log2ratio_siNipbl_siSA_ChIPs_newPseudocount/sorted_regions_log2desc/10groups/CHMM_18states"

annotations <- read.table(file.path(annoDir,"state_annotations_simplified.txt"), 
	header=FALSE, sep="\t") %>% 
	`names<-`(c("state","anno","rgb")) %>%
	dplyr::mutate_at(c("anno"), as.factor) %>% 
	as_tibble()

cols <- colnames(scoresTADs[[1]] %>% data.frame())
cols <- c(cols, "exp")
toGG <- matrix(nrow=0, ncol= length(cols) ) %>% `colnames<-`(cols)

for (name in names(scoresTADs)) {
	GR <- scoresTADs[[name]]
	toGG <- GR %>% as.data.frame() %>% as_tibble() %>%
		dplyr::mutate(exp=name) %>% rbind(toGG, .)
}

toGG <- dplyr::left_join(toGG, annotations, by = c('state'="state"))
toGG <- toGG %>% tidyr::separate(exp, into=c("cond", "positions"),
	sep="_AT_") %>%
	dplyr::mutate(positions=gsub("siNipbl_merged_.+","peaks", positions),
		cond= gsub("_log2Ratio.*", "", cond)) %>%
	dplyr::mutate_at(c("cond", "positions"), as.factor) 

summary(toGG)

## create mapping between states and colors
char_toRGB <- function(rgb_str) {
	sapply(strsplit(rgb_str, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue = 255)
)}
color <- unique(annotations[, c(2,3)]) %>% dplyr::mutate(fill=char_toRGB(rgb))

data <- toGG 

labeller <- data %>% group_by(cond) %>% tally %>%
	dplyr::mutate(label=paste0(cond, "\n(n = ", n, ")")) %>%
	dplyr::select(!n) %>%
	deframe

per_cond_anno <- data %>% dplyr::group_by(cond, anno, rgb) %>% tally()

p1 <- ggplot(data, aes(x=anno, y=score)) + 
	theme_bw() + 
	stat_boxplot(geom = "errorbar", width = 0.5) +
	geom_boxplot(aes(fill = anno), position= "dodge", width= 0.5,
		outlier.shape= 1, outlier.size= 0.3, outlier.alpha=0.2, notch= TRUE,
		varwidth= FALSE, width=0.2) +
	stat_summary(fun= mean, geom= "point", size= rel(0.75), alpha= .8) +
	scale_y_continuous(limits = c(-3, 5)) +
	geom_label(data=per_cond_anno, 
		aes(x=anno, label= n, fill=anno),
		y = 4.5,
		hjust=0.5, vjust= 0.5, size = 2.2,
		label.padding = unit(0.2, "lines"),
		label.r = unit(0.15, "lines"),
		label.size = 0.4,
		colour= "white"
		) +
	scale_fill_manual(
		limits=color$anno,
		values=color$fill, 
		labels=color$anno) +
	facet_wrap(~cond, ncol= length(unique(data$cond)),
		labeller = labeller(cond = labeller)) +
	labs(x="", y="Log2Ratio score", fill="State") +
	theme(axis.text.x=element_blank(),
		panel.grid= element_blank()) +
	ggsave("boxplot_scores_per_state_siNIPBL@siNIPBL_peaks.png", width= 7.3, height= 4.5)

