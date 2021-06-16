###############################################################
####						TAD BORDERS  				   ####
##
##  Intersect peaks divided in 10 groups sorted according to
##  desc log2FC with TADs at 1bp resolution. 
##
##  Then intersect wach peak with the TADs and move to 
##  "plot_TADcentreDistance.R"
##  
##
###############################################################

###############################################################
####						TAD MIDDLE  				   ####
###############################################################
RESDIR="/Users/aqo/Desktop/cluster/MCF10A/siNipbl-siSA_combined/log2ratio_siNipbl_siSA_ChIPs_newPseudocount/TAD_borders"
mkdir -p ${RESDIR}

TADS="/Users/aqo/Desktop/cluster/MCF10A/TADS-MCF10A/TAD_mcf10a_hg19_middle_and_length.sorted.bed"

head ${TADS}

WORKDIR="/Users/aqo/Desktop/cluster/MCF10A/siNipbl-siSA_combined/log2ratio_siNipbl_siSA_ChIPs_newPseudocount/sorted_regions_log2desc/10groups/CHMM_18states"
cd ${WORKDIR}
 
head *.SortedRegions.bed_10_groups.bed_model18_intersected.bed

for f in ${RESDIR}/*10_groups_sorted.bed
do
	echo -e ">> $(wc -l ${f})"

	bedtools intersect -wa -wb -a ${f} -b ${TADS} >\
	${f/10_groups_sorted.bed/INTERSECT_TADmiddle.bed} && 

	echo -e "-- $(wc -l ${f/10_groups_sorted.bed/INTERSECT_TADmiddle.bed}) \n"
done

