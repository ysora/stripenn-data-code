# mm10 stripe

workdir=/mnt/data0/sora/stripe/human_mouse_conserve
cd $workdir

mm10=/mnt/data0/sora/stripe/stripenn_out/human_mouse_conserve/mouse_CD4/Th1/result_filtered.tsv
mapchain=/mnt/data0/sora/reference/liftOver/mm10ToHg38.over.chain.gz

awk ' NR > 1 && $11 < 0.05 && $12 > 0 { if ( $2 == $5 ) print "chr"$1"\t"$2"\t"$2+49999 ; else print "chr"$1"\t"$3-49999"\t"$3}' $mm10 > mm10_anchor_tsv
awk ' NR > 1 && $11 < 0.05 && $12 > 0 { print "chr"$4"\t"$5"\t"$6 }' $mm10 > mm10_domain.tsv

liftOver mm10_anchor.tsv $mapchain liftOver_mm10_to_hg38_anchor.tsv liftOver_mm10_to_hg38_anchor_UnMapped.tsv -minMatch=0.1
liftOver mm10_domain.tsv $mapchain liftOver_mm10_to_hg38_domain.tsv liftOver_mm10_to_hg38_domain_UnMapped.tsv -minMatch=0.1

# PairToPair

hg38=/mnt/data0/sora/stripe/stripenn_out/human_mouse_conserve/human_CD4/result_filtered.tsv

awk ' NR > 1 && $11 < 0.05 && $12 > 0 { if ( $2 == $5 ) print "chr"$1"\t"$2"\t"$2+49999 ; else print "chr"$1"\t"$3-49999"\t"$3}' $hg38 > hg38_anchor_tsv
awk ' NR > 1 && $11 < 0.05 && $12 > 0 { print "chr"$4"\t"$5"\t"$6 }' $hg38 > hg38_domain.tsv
bedtools intersect -a hg38_anchor.tsv -b liftOver_mm10_to_hg38_anchor.tsv -wb > anchor_overlap.tsv
bedtools intersect -a hg38_domain.tsv -b liftOver_mm10_to_hg38_anchor.tsv -wb > mm_anchor_to_hg_domain_overlap.tsv
bedtools intersect -a hg38_domain.tsv -b liftOver_mm10_to_hg38_anchor.tsv -wa > mm_anchor_to_hg_domain_overlap_domain.tsv

mapchain=/mnt/data0/sora/reference/liftOver/hg38ToMm10.over.chain.gz
liftOver mm_anchor_to_hg_domain_overlap.tsv $mapchain mm_anchor_to_hg_domain_overlap_to_mm10.tsv mm_anchor_to_hg_domain_overlap_to_mm10_UnMapped.tsv -minMatch=0.1

#bedtools intersect -a mm_anchor_to_hg_domain_overlap_to_mm10.tsv -b 
