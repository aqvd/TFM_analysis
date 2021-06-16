#####################################################################
###             HEATMAP Cohesin at SUPER ENHANCERS 				 ####
##
##  Plot cohesin signal at super enhancer regions from Zhang et al.
#####################################################################

bw="\
/storage/scratch01/users/aquevedo/6Nov-siNipbl-newSnak/res/bw/SMC1_siC_mean.bw \
/storage/scratch01/users/aquevedo/6Nov-siNipbl-newSnak/res/bw/SMC1_siNipbl_S7_RPKM_scaled.bw \
/storage/scratch01/users/aquevedo/6Nov-siNipbl-newSnak/res/bw/SA2_siC_mean.bw \
/storage/scratch01/users/aquevedo/6Nov-siNipbl-newSnak/res/bw/SA2_siNipbl_S6_RPKM_scaled.bw \
/storage/scratch01/users/aquevedo/9Nov-siSA1_siSA2/res/bw/SA2_siC_mean.bw \
/storage/scratch01/users/aquevedo/9Nov-siSA1_siSA2/res/bw/SA2_siSA1_mean.bw \
/storage/scratch01/users/aquevedo/9Nov-siSA1_siSA2/res/bw/SA2_siSA2_mean.bw \
/storage/scratch01/users/aquevedo/dani_files/GSE107176_H3K27ac_final.bw \
"

reg="/storage/scratch01/users/aquevedo/MCF10A/super-enhancers_MCF10A-TableS1-Zhang2019-BreastCancerRes.bed"
mat="/storage/scratch01/users/aquevedo/snakemakeTest/res/deeptools/matrix_enhancers-siNipblSE.gz"
conda activate deeptools

medianSE=19000
binSE=95

sbatch -o ~/matrixSE.log \
-e ~/matrixSE.txt \
-J matrixSE \
-N 1 \
-c 12 \
--mem=12G \
-t 40 \
--wrap "computeMatrix scale-regions \
--regionBodyLength ${medianSE} \
--startLabel 5' --endLabel 3' \
--scoreFileName ${bw} \
--regionsFileName ${reg} \
--upstream 1805 --downstream 1805 \
--binSize ${binSE} \
--missingDataAsZero \
--smartLabels \
--numberOfProcessors 12 \
--outFileName ${mat} \
"

squeue -u aquevedo
less matrixSE.txt

## ==================== HEATMAP ====================== ##
workdir="/Users/aqo/Desktop/MCF10A/enhancers-MCF10A/"
matSE="/Users/aqo/Desktop/MCF10A/enhancers-MCF10A/matrix_enhancers-siNipblSE.gz"
sampLab="NIPBL_siC NIPBL_siNipbl SA1_siC SA1_siNipbl SA2_siC SA2_siC-2 SA2_siNipbl SMC1_siC SMC1_siC-2 SMC1_siNipbl"

heat="/Users/aqo/Desktop/MCF10A/enhancers-MCF10A/heatmap-siNIPBLsiSA_merged_replicates_chips@SuperEnhancersMCF10A-1.png"
echo ${heat/?.png/clusteredRegions.bed}

ls ${matSE}

plotHeatmap \
--matrixFile ${matSE} \
--outFileName ${heat} \
--kmeans 7 \
--startLabel "5'" --endLabel "3'" \
--colorMap Reds \
--plotType lines \
--heatmapHeight 36 \
--heatmapWidth 6 \
--xAxisLabel "distance (kb)" \
--refPointLabel "Summit" \
--outFileSortedRegions ${workdir}${heat/?.png/clusteredRegions.bed}

ls -lh ${workdir}${heat/?.png/clusteredRegions.bed}
head ${workdir}${heat/?.png/clusteredRegions.bed}

##first 6 clusters are unspeciffic signal, so remove them using awk
gawk -F "\t" 'NR==1{print}; /cluster_7/' ${workdir}${heat/?.png/clusteredRegions.bed} > ${workdir}${heat/?.png/cluster_7.bed}

## Check
ls -lh ${workdir}${heat/?.png/cluster_7.bed}
head ${workdir}${heat/?.png/cluster_7.bed}

## ============ re-plot using this clean regions =========== ##
reg="/storage/scratch01/users/aquevedo/MCF10A/heatmap-siNIPBLsiSA_merged_replicates_chips@SuperEnhancersMCF10A-cluster_7.bed"
mat="/storage/scratch01/users/aquevedo/MCF10A/matrix-siNIPBLsiSA_merged_replicates_chips@SuperEnhancersMCF10A-cluster_7.gz"

conda activate deeptools

sbatch \
-J matrixSE \
-N 1 \
-c 12 \
--mem=12G \
-t 40 \
--wrap "computeMatrix scale-regions \
--regionBodyLength ${medianSE} \
--startLabel 5' --endLabel 3' \
--scoreFileName ${bw} \
--regionsFileName ${reg} \
--upstream 9500--downstream 9500 \
--binSize ${binSE} \
--missingDataAsZero \
--smartLabels \
--numberOfProcessors 12 \
--outFileName ${mat} |& tee ${mat}.log \
"

squeue -u aquevedo
less matrixSE.txt

## Heatmap
heat="/Users/aqo/Desktop/MCF10A/enhancers-MCF10A/heatmap-siNIPBLsiSA_merged_replicates_chips@SuperEnhancersMCF10A-1.png"

plotHeatmap \
--matrixFile ${mat} \
--outFileName ${heat} \
--startLabel "5'" --endLabel "3'" \
--colorMap Reds Reds Reds Reds Blues Blues Blues Purples \
--plotType lines \
--heatmapHeight 36 \
--heatmapWidth 6 \
