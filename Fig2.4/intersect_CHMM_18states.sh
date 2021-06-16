######################################################
## siNIPBL and siSATAGs CHIP
## Intersect Log2FC sorted, 10 groups divided, with chromatin States
##
##
######################################################
RATIODIR="/Users/aqo/Desktop/cluster/MCF10A/siNipbl-siSA_combined/log2ratio_siNipbl_siSA_ChIPs_newPseudocount/sorted_regions_log2desc/10groups"
WORKDIR="${RATIODIR}/CHMM_18states"

mkdir -p ${WORKDIR}
cd ${WORKDIR}

## Dense .BED chromatin 18 states
CHMM="/Users/aqo/Desktop/cluster/MCF10A/chromHMM-MCF10A/model-18states/MCF-10A_model18states_dense.bed"
head ${CHMM}
tail ${CHMM}

## Log2Ratio bed
ratios=($(ls ${RATIODIR}/*10_groups.bed | grep -Ev "/NiPBL"))

echo "${ratios[@]}"
## number of files
echo ${#ratios[@]}

for r in "${ratios[@]}"
do
	#r=$(head -n 2000 ${r})
	head ${r}
	echo "=============="
	head ${CHMM}
	echo "============="

	bedtools intersect -wa -wb -a ${r} -b ${CHMM} | \
	awk -F "\t" 'BEGIN{ OFS="\t"; print "#chr","start","end","group","state";};
	{ print $1,$2,$3,$4,$8}'> \
	${WORKDIR}/${r##*/}_model18_intersected.bed && echo "Done..."
done

## check files generated
ls -lh ${WORKDIR}
wc -l ${WORKDIR}/* | sed -E 's|/.+/||g'
head ${WORKDIR}/*

## Check the same number of regions present
for r in "${ratios[@]}"
do
	wc -l ${r} ##ratios files have header, so 1 line more!!!
done

## Count the number of states present
for f in ${WORKDIR}/*
do
	echo ">> ${f##*/}"
	cut -d$'\t' -f 4 ${f} | sort -n | uniq -c
done

## FIMISHED, move to R
