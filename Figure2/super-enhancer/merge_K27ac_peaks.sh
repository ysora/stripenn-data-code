#!/bin/bash
# Sora Yoon (Modified Golnaz Vahedi's code)
# Date: September 14th, 2021
# Figure 3: the goal here is to find super-enhancers from K27ac peaks merged by 12.5kb following Rick Young's approach
# merge_K27ac_peaks.sh combined two technical replicates and merged the final peaks- I then counted tags at merged peaks and used MapSEs_YoungMethod.R to calculate super and typical enhancers
#
mydir=/mnt/data0/sora/stripe/human_mouse_conserve/Super-enhancers
Tempdir=/mnt/data0/sora/stripe/human_mouse_conserve/Super-enhancers/TempFiles
SEdir=/mnt/data0/sora/stripe/human_mouse_conserve/Super-enhancers/SE

BEDDIR=/mnt/alvand/VahediLab/GenomicsData/Mouse/T_lymphocytes/ChIPseq/Thymus/04.macs2/
bamDir=/mnt/alvand/VahediLab/GenomicsData/Mouse/T_lymphocytes/ChIPseq/Thymus/02.bam/
#
windowsize=50
mergesize=12500

cd ${mydir}
mkdir ${Tempdir}
# Merge peaks of two repeats and those within 12.5kbp
filename=ChIP-seq_K27ac_DP_BL6
cat $BEDDIR/ChIP-seq_K27ac_DP_BL6_Rep1_36160142_S1_peaks.bed $BEDDIR/ChIP-seq_K27ac_DP_BL6_Rep2_36160142_S2_peaks.bed | sortBed -i - | mergeBed -i - -d $mergesize > ${mydir}/${filename}_merged.bed

filename=ChIP-seq_K27ac_DP_NOD
cat $BEDDIR/ChIP-seq_K27ac_DP_NOD_Rep1_36160142_S3_peaks.bed $BEDDIR/ChIP-seq_K27ac_DP_NOD_Rep2_36160142_S4_peaks.bed | sortBed -i - | mergeBed -i - -d $mergesize > ${mydir}/${filename}_merged.bed
#-----------------------------------------------------------------------------
# Count tags for height analysis in SE detection
REFFILE=${mydir}/ChIP-seq_K27ac_DP_BL6_merged.bed
cd $bamDir
for i in `ls *K27ac*BL6*.bam`;
do
 filename="${i%.*}"
 echo $filename
 bedtools coverage -counts -a $REFFILE -b ${filename}.bam  > ${Tempdir}/${filename}_merged.bedgraph
done

REFFILE=${mydir}/ChIP-seq_K27ac_DP_NOD_merged.bed
cd $bamDir
for i in `ls *K27ac*NOD*.bam`;
do
 filename="${i%.*}"
 echo $filename
 bedtools coverage -counts -a $REFFILE -b ${filename}.bam  > ${Tempdir}/${filename}_merged.bedgraph
done
