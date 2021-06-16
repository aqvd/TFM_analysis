##########################################################################
## HEATMAPS log2FC cohesin and STAGs in positions_ctcf, psitions_no_CTCF
##
##########################################################################

basedir="/storage/scratch01/users/aquevedo/log2ratio_siSTAG_siNIPBL_ChIPs_NewPseudocount"
ls $basedir
wd="${basedir}/cohesinCTCF_noCTCF"
mkdir -p ${wd}/log;cd ${wd}

## ============================== MATRIX ================================ ##
bwdir="${basedir}"
ls ${bwdir}/*log2ratio.bw | xargs -n 1 basename

bws="\
${bwdir}/SMC1_siNipbl_log2ratio.bw \
${bwdir}/SA2_siNipbl_log2ratio.bw \
${bwdir}/SA1_siNipbl_log2ratio.bw \
"

ls -lh ${bws}

bwsLab="\
SMC1-siNipbl \
SA2-siNipbl \
SA1-siNipbl \
"

reg="\
/storage/scratch01/users/aquevedo/dani_files/positions_ctcf \
/storage/scratch01/users/aquevedo/dani_files/positions_no_ctcf \
"

wc -l ${reg}

regLab="\
\"Positions CTCF (24912)\" \
\"Positions No CTCF (14607)\" \
"

echo $regLab

matFilename="matrix-log2ratio_siNIPBL_siSTAGs-@positions_ctcf_no_ctcf.gz"
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
--numberOfProcessors 15"

echo "${comm} |& tee log/${matFilename}.log"

conda activate deeptools

sbatch \
-J matrix \
-N 1 \
-c 15 \
--mem=14G \
-t 120 \
--wrap "${comm} |& tee log/${matFilename/.gz/.log}"

squeue -u aquevedo

less "log/${matFilename/.gz/.log}"

ls -lah ${mat}

## ================= HEATMAP ==================== ##
matFilename="matrix-log2ratio_siNIPBL_siSTAGs-@positions_ctcf_no_ctcf.gz"
mat="${wd}/${matFilename}"

heatFilename="heatmap-log2ratio_siNIPBL_siSTAGs-@positions_ctcf_no_ctcf-1.png"
heat="${wd}/${heatFilename}"


comm="\
plotHeatmap \
--matrixFile ${mat} \
--samplesLabel ${bwsLab} \
--regionsLabel ${regLab} \
--colorMap bwr \
--yMin -2.5 --yMax 2.5 \
--zMin -2 --zMax 2 \
--outFileName ${heat} \
--plotType lines \
--heatmapHeight 25 \
--heatmapWidth 6 \
--xAxisLabel \"\" \
--refPointLabel \"\" |& tee log/${heatFilename/.png/.log} \
"

echo "${comm}"

sbatch \
-J heat \
-N 1 \
-c 1 \
--mem=6G \
-t 5 \
--wrap "${comm}"

squeue -u aquevedo

less log/${heatFilename/.png/.log}
ls ${heatFilename}


