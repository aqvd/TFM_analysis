############################################################
##			BOXPLOTS Score Dynamics
##		
## log2FC soreted regions siNIPBL ChIP are intersected with chromatin states
## and divided in 10 groups
##
## This script is used to plot the nuber of peaks at each state from group 1 to 10
## Uses the complete chromatin states model with 18 states.	
##
############################################################


library(tidyverse)

wd<-"/Users/aqo/Desktop/cluster/MCF10A/siNipbl-siSA_combined/log2ratio_siNipbl_siSA_ChIPs_newPseudocount/sorted_regions_log2desc/10groups/CHMM_18states"
setwd(wd)

plotDir <- file.path(wd, "plots_fullModel")
dir.create(plotDir, showWarnings = FALSE)

NGROUPS <- 10 # needed for plot legends

## List files
log2files <- list.files(wd, pattern = "intersected.bed")

## read files into tables
bed_l <- plyr::llply(log2files, function(file){ 
	read.table(file, sep="\t", header = TRUE, comment.char="") %>% as_tibble() %>%
		dplyr::mutate_at(c("state","group"), as.factor)
})

## Name list with the name of the protein where log2Reatio comes from
names(bed_l) <- gsub("_log2Ratio.+","",log2files)

## =========== Add state annotations ============= ##
annotations <- read.table(file.path(wd,"state_annotations_byHand.txt"), 
	header=FALSE, sep="\t") %>% 
	`names<-`(c("state","anno","rgb")) %>%
	dplyr::mutate_at(c("state","anno"), as.factor) %>% 
	as_tibble()

bed_l <- lapply(bed_l, function(bed){
	bed %>% dplyr::left_join(.,annotations, by=c("state"="state"))
})
bed_l[[1]] 

#####################################################
####       BARPLOTS 18states COLOr STRIPS        ####
#####################################################

char_toRGB <- function(rgb_str) {
	sapply(strsplit(rgb_str, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue = 255)
)}

gtable_select <- function (x, ...) {
  matches <- c(...)
  x$layout <- x$layout[matches, , drop = FALSE]
  x$grobs <- x$grobs[matches]
  x
}

gtable_stack <- function(g1, g2){
  g1$grobs <- c(g1$grobs, g2$grobs)
  g1$layout <- rbind(g1$layout, g2$layout)
  g1
}
## =================================================================== ##

for (bedName in names(bed_l)){
	## get total number of peaks per group
	nPerGroup <- bed_l[[bedName]] %>% group_by(group) %>%
		dplyr::summarize(n=n())  %>% .[1,2]
	
	## Name to plot in title and legend
	name<-gsub("_"," ", bedName)
	v_name <- str_split(bedName,"_") %>% unlist()

	legendName <- paste("Log2Ratio in",v_name[1],v_name[2], sep=" ")
	plotTitle <- paste("Log2Ratio distribution of peaks in",
		v_name[1],v_name[2],
		"\n Peaks per group:", nPerGroup, sep=" ")

	## Filename
	filename <- paste0(bedName,"_number_peaks_18chromatinStates.pdf")

	data <- bed_l[[bedName]] %>% dplyr::mutate(facet_fill_color=char_toRGB(rgb))
	print(head(data))
	
	## Main plot
	mainP <- ggplot(data, aes(x=group,y=as.numeric(state)/as.numeric(state))) +
		geom_col(aes(fill=group), alpha=1) + 
		theme_bw() + 
		scale_fill_manual(
			values= RColorBrewer::brewer.pal(NGROUPS,"RdBu"),
			breaks=factor(1:NGROUPS),
			limits=factor(1:NGROUPS),
			name=legendName,
			labels=c("+++",rep("",NGROUPS-2),"---")) +
		scale_y_continuous(expand = expansion(mult = c(0, .03))) + 
		scale_x_discrete(
			breaks=1:10,
			limits=factor(1:NGROUPS),
			labels=c("+++",rep("",NGROUPS-2),"---")) +
		labs(x="", y="Number of Peaks") + 
		theme(legend.position = "none",
			panel.grid.major.x = element_blank(),
			panel.grid.major.y = element_line(linetype=2, color="grey80")) +
		ggtitle(plotTitle) + 
		facet_wrap(~anno, scales="free", ncol=6) +
	    theme(strip.background=element_blank())
		
	## create mapping between states and colors
	color <- unique(annotations[, c(2,3)]) %>% dplyr::mutate(fill=char_toRGB(rgb))

	## create dummy data frame, filled with color ~ state
	dummy <- ggplot(data, aes(x=anno)) +
		facet_wrap(~anno, scales="free", ncol=6) +
		geom_rect(aes(fill = anno),
			xmin=-Inf, ymin=-Inf, xmax=Inf, ymax=Inf) +
		scale_fill_manual(limits=color$anno,
			values=color$fill, 
			labels=color$anno) +
		theme_minimal()
	
	## Transform plots into grid tables
	library(gtable)
	library(grid)

	g1 <- ggplotGrob(mainP)
	g2 <- ggplotGrob(dummy)

	## Move strip panel 1 unit down, where plot is filled with one colour
	panels <- grepl(pattern="panel", g2$layout$name)
	strips <- grepl(pattern="strip-t", g2$layout$name)
	g2$layout$t[panels] <- g2$layout$t[panels] - 1
	g2$layout$b[panels] <- g2$layout$b[panels] - 1

	## Draw the main plot and over it the strips (facet labels) that have
	## the colour we want
	new_strips <- gtable_select(g2, panels | strips)
	new_plot <- gtable_stack(g1, new_strips)

	pdf(file.path(plotDir,filename), width=12, height=8, )
		grid.draw(new_plot)
	dev.off()
}	

