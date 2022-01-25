workDir=Data/

COOL=GSE82144_Kieffer-Kwon-2017-activated_B_cells_72_hours_WT_30_10000.cool

CRS_BED=${workDir}chromosight_all_merged.bed
ZEB_BED=${workDir}zebra_all.bed
STR_BED=${workDir}stripenn_all.bed
DCR_BED=${workDir}dcr_all.bed
CRS_OUT=chromosight_all_merged.out
ZEB_OUT=zebra_all.out
STR_OUT=stripenn_all.out
DCR_OUT=dcr_all.out


coolpup.py $COOL $CRS_BED --outname $CRS_OUT --local --rescale --unbalanced
coolpup.py $COOL $ZEB_BED --outname $ZEB_OUT --local --rescale --unbalanced
coolpup.py $COOL $DCR_BED --outname $DCR_OUT --local --rescale --unbalanced
coolpup.py $COOL $STR_BED --outname $STR_OUT --local --rescale --unbalanced

plotpup.py $STR_OUT $CRS_OUT $ZEB_OUT $DCR_OUT --output pileup.pdf --n_cols 4

