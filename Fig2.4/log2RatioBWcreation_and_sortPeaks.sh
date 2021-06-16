##############################################################
##    Sort descending log2FC combined summits ChIP siNIPBL
##
##  To generate the plots of score dynamics
##
##  1- Use this script to sort genomic positions in a .bed file
##  2- use "1"0groups_log2FC_sortedRegions.sh" to  divide .bed file
##	   in 10 equally sized groups
##  3- Use "intersect_CHMM_18states.sh" to annotate groups with 
##     states
##  4- Finally, plot using "boxplots_scoreDynamics.R"


resdir="/storage/scratch01/users/aquevedo/log2ratio_siSTAG_siNIPBL_ChIPs_NewPseudocount"
mkdir -p ${resdir}/log

###########################################
#####            siNIPBL              #####
###########################################
## SA1
cd /storage/scratch01/users/aquevedo/6Nov-siNipbl-newSnak/res/bw/
  
sbatch -o ${resdir}/log/SA1_siNipbl_log2ratio.log \
-e ${resdir}/log/SA1_siNipbl_log2ratio.txt \
-J bwCompare \
-N 1 \
-c 11 \
--mem=10G \
-t 30 \
--wrap "
bigwigCompare \
-b1 SA1_siNipbl_S5_RPKM_scaled.bw \
-b2 SA1_siC_S1_RPKM_scaled.bw  \
--pseudocount 0.005 \
--operation log2 \
-p 11 --outFileFormat bigwig \
--outFileName ${resdir}/SA1_siNipbl_log2ratio.bw "

## SA2
ls SA2*

sbatch -o ${resdir}/log/SA2_siNipbl_log2ratio.log \
-e ${resdir}/log/SA2_siNipbl_log2ratio.txt \
-J bwCompare \
-N 1 \
-c 11 \
--mem=10G \
-t 30 \
--wrap "
bigwigCompare \
-b1 SA2_siNipbl_S6_RPKM_scaled.bw \
-b2 SA2_siC_mean.bw  \
--pseudocount 0.005 \
--operation log2 \
-p 11 --outFileFormat bigwig \
--outFileName ${resdir}/SA2_siNipbl_log2ratio.bw "

squeue -u aquevedo

## SMC1 
sbatch -o ${resdir}/log/SMC1_ratio.log \
-e ${resdir}/log/SMC1_ratio.txt \
-J ratioBW \
-N 1 \
-c 13 \
--mem=10G \
-t 30 \
--wrap "
bigwigCompare \
-b1 SMC1_siNipbl_S7_RPKM_scaled.bw \
-b2 SMC1_siC_mean.bw  \
--pseudocount 0.005 \
--operation log2 \
-p 11 --outFileFormat bigwig \
--outFileName ${resdir}/SMC1_siNipbl_log2ratio.bw "

squeue -u aquevedo
ls -lh ${resdir}

for f in /storage/scratch01/users/aquevedo/log2ratio_siSA/*_merged_unique_summits.bed; do
	echo ${f}
	cp ${f} ${resdir}
done

for f in /storage/scratch01/users/aquevedo/log2ratio_siNipbl/*_merged_unique_summits.bed; do
	echo ${f}
	cp ${f} ${resdir}
done


##########################################################################
## MATRIX of Log2Ratio Sorted descending with unique summits per prot   ##
##########################################################################
##################### siNiPBL ################################
cd ${resdir}

prots=("SA1" "SA2" "SMC1")

for p in "${prots[@]}" 
do
	echo ">>=========== $p ==========<<"
	## Get .bw
	bw="${p}_siNipbl_log2ratio.bw"
	echo "BIG WIG: 
	${bw}"
	ls -l ${bw}
	## Get summits per prot
	reg="${p}_siNipbl_merged_unique_summits.bed"
	echo "REGIONS: 
	${reg}"
	ls -l ${reg}

	## Generate command to compute matrix
	mat="Matrix-${p}_siNipbl_log2Ratio@${p}_siNipbl_merged_unique_summits-SortedDown.gz"
	sortReg="${p}_siNipbl_log2Ratio@${p}_siNipbl_merged_unique_summits-SortedRegions.bed"
	echo -e "\n>> MATRIX:\n${mat##*/}"

	comm="\
	computeMatrix reference-point \
	--referencePoint center \
	--scoreFileName ${bw} \
	--sortRegion descend \
	--sortUsing mean \
	--samplesLabel ${p} \
	--regionsFileName ${reg} \
	--outFileName ${mat} \
	--outFileSortedRegions ${sortReg} \
	--upstream 100 --downstream 100 \
	--binSize 50 \
	--missingDataAsZero \
	--numberOfProcessors 5"

	echo -e "\n>>COMMAND MATRIX:\n${comm}\n"

	# ============== MATRIX comment submission command once matrix is computed
	if [[ ! -e ${mat} ]]
	then
		echo ">> MATRIX MUST BE COMPUTED <<"
		sbatch \
		-J mat${p} \
		-N 1 \
		-c 5 \
		--mem=10G \
		-t 60 \
		--wrap "${comm} |& tee log/${mat##*/}.log"
	else
		echo ">> MATRIX EXISTS, heatmap will be plotted <<"
	fi

	## ============== HEATMAP
	heatFilename="heatmap-${p}_siNipbl_log2Ratio@${p}_siNipbl_merged_unique_summits-1.png"
	heat="${heatFilename}"

	bwsLab="\
	\"${p}-siNipbl Log2Ratio\" \
	" 

	regLab="\
	\"${p} CHIP siNipbl merged unique peaks\" \
	"
	echo $regLab

	comm="\
	plotHeatmap \
	--matrixFile ${mat} \
	--samplesLabel ${bwsLab} \
	--colorMap bwr \
	--zMax 5 --zMin -5 \
	--regionsLabel ${regLab} \
	--outFileName ${heat} \
	--legendLocation none \
	--plotType lines \
	--heatmapHeight 15 \
	--heatmapWidth 4 \
	--xAxisLabel \"\" \
	--refPointLabel \"\" |& tee log/${heatFilename}.log \
"
	echo -e "\n>>COMMAND HEATMAP:\n${comm}\n"
	
	if [[ -e ${mat} ]] ## Si la matriz existe
	then
		echo ">> PLOTTING HEATMAP <<"
		# Comment until matrix is computed
		sbatch \
		-J heat${p} \
		-N 1 \
		-c 1 \
		--mem=8G \
		-t 30 \
		--wrap "${comm}"
	else
		echo ">> CAN'T PLOT HEATMAP WITHOUT MATRIX"
	fi
done

squeue -u aquevedo

