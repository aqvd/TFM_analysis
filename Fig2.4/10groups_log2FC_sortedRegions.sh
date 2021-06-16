## Just run script to divide descending log2FC sorted regins in 10 groups.
## the script is available at github aqvd/snakemake_workflows/Chip-Seq/scripts

basedir="/Users/aqo/Desktop/cluster/MCF10A/siNipbl-siSA_combined/log2ratio_siNipbl_siSA_ChIPs_newPseudocount/sorted_regions_log2desc"
cd ${basedir}

## Divide sorted regions in groups
~/Desktop/snakemake_workflows/Chip-Seq/scripts/deepToolsRegions_Ngroups.R \
--outDir 10groups --dir $(pwd) --nGroups 10 --regionsRegex ".+SortedRegions.?[.]bed"






