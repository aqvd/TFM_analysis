############################################################
## MATRIX & HEATMAP  of cohesin distribution arroun TSS of
## genes classified as Up, Down and 5000 random genes
##
## Genomic coordenates for those positions are generated 
## in "TSS_hg19_forHeatmapsCohesin.R"
############################################################

wd="/storage/scratch01/users/aquevedo/6Nov-siNipbl-newSnak/heatmaps_DEGs_TSS"
mkdir -p ${wd}/log; cd ${wd}

## ============================== MATRIX ================================ ##
bws="\
/storage/scratch01/users/aquevedo/dani_files/MCF10A_CTCF_treat_pileup.bw \
/storage/scratch01/users/aquevedo/6Nov-siNipbl-newSnak/res/bw/SMC1_siC_mean.bw \
/storage/scratch01/users/aquevedo/6Nov-siNipbl-newSnak/res/bw/SA2_siC_mean.bw \
/storage/scratch01/users/aquevedo/6Nov-siNipbl-newSnak/res/bw/SA1_siC_S1_RPKM_scaled.bw \
/storage/scratch01/users/aquevedo/6Nov-siNipbl-newSnak/res/bw/SMC1_siNipbl_S7_RPKM_scaled.bw \
/storage/scratch01/users/aquevedo/6Nov-siNipbl-newSnak/res/bw/SA2_siNipbl_S6_RPKM_scaled.bw \
/storage/scratch01/users/aquevedo/6Nov-siNipbl-newSnak/res/bw/SA1_siNipbl_S5_RPKM_scaled.bw \
/storage/scratch01/users/aquevedo/BRD4_PRJNA295370/res/bw/BRD4_noTreat_RPKM.bw \
/storage/scratch01/users/aquevedo/dani_files/GSE107176_H3K27ac_final.bw \
/storage/scratch01/users/aquevedo/dani_files/GSE107176_H3K27me3_final.bw \
"

ls -l $bws

bwsLab="\
CTCF \
SMC1-siC \
SA2-siC \
SA1-siC \
SMC1-siNipbl \
SA2-siNipbl \
SA1-siNipbl \
BRD4 \
H3K27ac \
H3K27me3 \
"

echo $bwsLab
regionsDir="${wd}/DEGs_TSS_bed"

ls ${regionsDir} | tr " " "\n"

reg="\
${regionsDir}/siNIPBL_up_TSS.bed \
${regionsDir}/siNIPBL_noDEG_TSS.bed \
${regionsDir}/siNIPBL_down_TSS.bed \
" 

wc -l ${reg}

regLab="\
\"UP siNIPBL (937)\" \
\"No DEGs (5000)\" \
\"DOWN siNIPBL (670)\" \
"

head ${reg}
echo $regLab

matFilename="matrix-CTCF_siNipblChIP_H3K27ac-me3@noDEG_UP_DOWN_siNIPBL.gz"
mat="${wd}/${matFilename}"

echo $mat

comm="\
computeMatrix reference-point \
--referencePoint center \
--scoreFileName ${bws} \
--samplesLabel ${bwsLab} \
--regionsFileName ${reg} \
--outFileName ${mat} \
--upstream 2500 --downstream 2500 \
--binSize 50 \
--missingDataAsZero \
--numberOfProcessors 10"

echo "${comm} |& tee log/${matFilename}.log"

conda activate deeptools

sbatch \
-J matrix \
-N 1 \
-c 10 \
--mem=12G \
-t 40 \
--wrap "${comm} |& tee log/${matFilename/.gz/.log}"

squeue -u aquevedo

less "log/${matFilename/.gz/.log}"

ls -lah ${mat}

## ================= HEATMAP ==================== ##
matFilename="matrix-CTCF_siNipblChIP_H3K27ac-me3@noDEG_UP_DOWN_siNIPBL.gz"
mat="${wd}/${matFilename}"

heatFilename="heatmap-CTCF_siNipblChIP_H3K27ac-me3@noDEG_UP_DOWN_siNIPB-1.png"
heat="${wd}/${heatFilename}"

bwsLab="\
CTCF \
SMC1-siC \
SA2-siC \
SA1-siC \
SMC1-siNipbl \
SA2-siNipbl \
SA1-siNipbl \
BRD4 \
H3K27ac \
H3K27me3 \
"

comm="\
plotHeatmap \
--matrixFile ${mat} \
--samplesLabel ${bwsLab} \
--regionsLabel ${regLab} \
--sortRegions descend \
--sortUsing mean \
--sortUsingSamples 9 \
--colorMap \
Oranges \
Reds Reds Reds Reds Reds Reds \
Greens RdPu Blues \
--zMax \
5 \
0.15 0.15 0.15 0.15 0.15 0.15 \
10 14 1.5 \
--zMin 0 \
--yMax \
6 \
0.2 0.2 0.2 0.2 0.2 0.2 \
20 22 2 \
--yMin 0 \
--outFileName ${heat} \
--plotType lines \
--heatmapHeight 20 \
--heatmapWidth 6 \
--xAxisLabel \"\" \
--refPointLabel \"\" |& tee log/${heatFilename/.png/.log} \
"

echo "${comm}"

sbatch \
-J heat \
-N 1 \
-c 1 \
--mem=4G \
-t 5 \
--wrap "${comm}"

squeue -u aquevedo

less log/${heatFilename/.png/.log}
ls ${heatFilename}

