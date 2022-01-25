gene_exp = read.delim('/mnt/data0/sora/stripe/GeneExpression/DPT/T_cell_NOD_BL_RNAseq_TPM.txt')
BL = apply(gene_exp[,3:4],1,mean)
NOD= apply(gene_exp[,1:2],1,mean)
DESeq = read.delim('/mnt/data0/sora/stripe/GeneExpression/DPT/DESeq/DESeq_results_with_IDs.txt')

stripe_gene = read.delim('/mnt/data0/sora/stripe/human_mouse_conserve/t1d_stripe_genes.tsv', header = F)
stripe_genes = unique(stripe_gene$V7)

h_stripe_genes = intersect(names(BL[which(BL>10)]),stripe_genes)

write(h_stripe_genes, '/mnt/data0/sora/stripe/human_mouse_conserve/t1d_stripe_genes_TPM_gt_10.tsv')

exp = gene_exp[which(rownames(gene_exp)%in%h_stripe_genes),]

# 2021 08 04
# sushi 
library(Sushi)

bed = read.delim('/mnt/data0/sora/stripe/GeneExpression/210720_stripe_TPM_gt_10_genes_coordinate.txt')
BL6 = read.delim('/mnt/data0/sora/stripe/testdata/DP_T_10000/chr5.matrix', header=T,row.names=1)
NOD = read.delim('/mnt/data0/sora/stripe/testdata/NOD_DP_T/chr5.matrix', header=T,row.names=1)
geneinfo = read.delim('/mnt/data0/golnaz/references/mouse/mm10_RefSeq.bed',header=F)
CHR = '5'
START = 53345001-50000
END = 	53680000+50000

system(
paste0("
BB=/mnt/data0/wenliang/software/UCSC/linux.x86_64/bigWigToBedGraph;
name=T1D_2;
gene=/mnt/data0/golnaz/references/mouse/mm10_RefSeq.bed;
folder=/mnt/alvand/VahediLab/GenomicsData/Mouse/T_lymphocytes/ChIPseq/Thymus/03.bigwig/;
CTCF_BL6=ChIP-seq_CTCF_DP_BL6_Rep2B_189927787_S1;
CTCF_NOD=ChIP-seq_CTCF_DP_NOD_Rep2B_189927787_S2;
SMC1_BL6=ChIP_seq_smc1_BL6_DP_rep2B_190921758_S2;
SMC1_NOD=ChIP_seq_smc1_NOD_DP_rep2B_190921758_S4;
H3K27ac_BL6=ChIP-seq_K27ac_DP_BL6_Rep2_36160142_S2;
H3K27ac_NOD=ChIP-seq_K27ac_DP_NOD_Rep2_36160142_S4;
ATAC_BL6=/mnt/alvand/VahediLab/GenomicsData/Mouse/T_lymphocytes/ATACseq/Thymus/05.bigwig/ATAC-seq_DP_170822_C57BL6_Rep3_33396376_MF_S5.bw;
ATAC_NOD=/mnt/alvand/VahediLab/GenomicsData/Mouse/T_lymphocytes/ATACseq/Thymus/05.bigwig/ATAC-seq_DP_170822_NOD_Rep3_33396376_MF_S9.bw;

chr=",CHR,";
start=",START,";
end=",END,";
cd /mnt/data0/sora/stripe/hicPlotter/NOD_vs_BL6/T1D_2/
$BB ${folder}${CTCF_BL6}.bw ${CTCF_BL6}.${name}.bedgraph -chrom=$chr -start=$start -end=$end;
$BB ${folder}${CTCF_NOD}.bw ${CTCF_NOD}.${name}.bedgraph -chrom=$chr -start=$start -end=$end;
$BB ${folder}${SMC1_BL6}.bw ${SMC1_BL6}.${name}.bedgraph -chrom=$chr -start=$start -end=$end;
$BB ${folder}${SMC1_NOD}.bw ${SMC1_NOD}.${name}.bedgraph -chrom=$chr -start=$start -end=$end;
$BB ${folder}${H3K27ac_BL6}.bw ${H3K27ac_BL6}.${name}.bedgraph -chrom=$chr -start=$start -end=$end;
$BB ${folder}${H3K27ac_NOD}.bw ${H3K27ac_NOD}.${name}.bedgraph -chrom=$chr -start=$start -end=$end;
$BB ${ATAC_BL6} ATAC-seq_BL6.${name}.bedgraph -chrom=$chr -start=$start -end=$end;
$BB ${ATAC_NOD} ATAC-seq_NOD.${name}.bedgraph -chrom=$chr -start=$start -end=$end")
)

snpbed=read.delim('/mnt/data0/sora/stripe/human_mouse_conserve/T1D_SNPs_in_BL6_stripes.bed',header=F)

bl6_smc1= read.delim('/mnt/data0/sora/stripe/hicPlotter/NOD_vs_BL6/T1D_2/ChIP_seq_smc1_BL6_DP_rep2B_190921758_S2.T1D_2.bedgraph', header=F)
bl6_ctcf= read.delim('/mnt/data0/sora/stripe/hicPlotter/NOD_vs_BL6/T1D_2/ChIP-seq_CTCF_DP_BL6_Rep2B_189927787_S1.T1D_2.bedgraph', header=F)
bl6_k27= read.delim('/mnt/data0/sora/stripe/hicPlotter/NOD_vs_BL6/T1D_2/ChIP-seq_K27ac_DP_BL6_Rep2_36160142_S2.T1D_2.bedgraph', header=F)
bl6_atac= read.delim('/mnt/data0/sora/stripe/hicPlotter/NOD_vs_BL6/T1D_2/ATAC-seq_BL6.T1D_2.bedgraph', header=F)

nod_smc1 = read.delim('/mnt/data0/sora/stripe/hicPlotter/NOD_vs_BL6/T1D_2/ChIP_seq_smc1_NOD_DP_rep2B_190921758_S4.T1D_2.bedgraph', header=F)
nod_ctcf = read.delim('/mnt/data0/sora/stripe/hicPlotter/NOD_vs_BL6/T1D_2/ChIP-seq_CTCF_DP_NOD_Rep2B_189927787_S2.T1D_2.bedgraph', header=F)
nod_k27 = read.delim('/mnt/data0/sora/stripe/hicPlotter/NOD_vs_BL6/T1D_2/ChIP-seq_K27ac_DP_NOD_Rep2_36160142_S4.T1D_2.bedgraph', header=F)
nod_atac = read.delim('/mnt/data0/sora/stripe/hicPlotter/NOD_vs_BL6/T1D_2/ATAC-seq_NOD.T1D_2.bedgraph', header=F)

rownames(BL6)  = rownames(NOD) = colnames(BL6) = colnames(NOD) = c(1:nrow(BL6))*10000-9999

#CHR = paste0('chr',CHR)
pdf('/mnt/data0/sora/stripe/hicPlotter/NOD_vs_BL6/T1D_2/R_figure.pdf', width = 10, height = 8)
plot.new()
par(cex=0.7, mai=c(.1,0.1,0.2,0.1))
par(fig=c(0.1,0.5,0.6,0.9))
phic_BL6 = plotHic(BL6, CHR,START,END, max_y = 50,
                   zrange=c(0,25), palette=colorRampPalette(c("white",'red'))
)
par(fig=c(0.1,0.5,0.51,0.65), new=T)
colnames(geneinfo)[4] = 'gene'
plotGenes(geneinfo[which(geneinfo$V1==paste0('chr',CHR) &geneinfo$V2>=START & geneinfo$V3<=END),],CHR,START,END ,
          maxrows=1,height=0.5,plotgenetype="arrow",bentline=TRUE,col="#465075", labeltext = TRUE,
          labeloffset= 0.4,fontsize=1.2, arrowlength = 0.01, bheight = 0.2)

par(fig=c(0.1,0.5,0.4,0.5), new=T)
plotBedgraph(bl6_ctcf ,CHR,START,END,color='#22854a',range=c(0,2),addscale = T)

par(fig=c(0.1,0.5,0.3,0.4), new=T)
plotBedgraph(bl6_smc1, CHR,START,END,color='#ff2185',range=c(0,3),addscale = T)

par(fig=c(0.1,0.5,0.2,0.3), new=T)
plotBedgraph(bl6_k27, CHR,START,END,color='#0600b5',range=c(0,5),addscale = T)

par(fig=c(0.1,0.5,0.1,0.2), new=T)
plotBedgraph(bl6_atac, CHR,START,END,color='black',range=c(0,1.8),addscale = T)

par(fig=c(0.1,0.5,0.05,0.1), new=T)
plotBed(beddata = snpbed, chrom = paste0('chr',CHR), chromstart = START, chromend = END)
labelgenome(CHR,START,END,n=5,scale="Mb")

###

par(fig=c(0.5,0.9,0.6,0.9), new=T)
phic_NOD = plotHic(NOD, CHR, START, END, max_y = 50,
                   zrange=c(0,25), palette=colorRampPalette(c("white",'red'))
)
par(fig=c(0.5,0.9,0.51,0.65), new=T)
colnames(geneinfo)[4] = 'gene'
plotGenes(geneinfo[which(geneinfo$V1==paste0('chr',CHR) &geneinfo$V2>=START & geneinfo$V3<=END),],CHR,START,END ,
          maxrows=1,height=0.5,plotgenetype="arrow",bentline=TRUE,col="#465075", labeltext = TRUE,
          labeloffset= 0.4,fontsize=1.2, arrowlength = 0.01, bheight = 0.2)

par(fig=c(0.5,0.9,0.4,0.5), new=T)
plotBedgraph(nod_ctcf ,CHR,START,END,color='#22854a',range=c(0,2),addscale = T)

par(fig=c(0.5,0.9,0.3,0.4), new=T)
plotBedgraph(nod_smc1, CHR,START,END,color='#ff2185',range=c(0,3),addscale = T)

par(fig=c(0.5,0.9,0.2,0.3), new=T)
plotBedgraph(nod_k27, CHR,START,END,color='#0600b5',range=c(0,5),addscale = T)

par(fig=c(0.5,0.9,0.1,0.2), new=T)
plotBedgraph(nod_atac, CHR,START,END,color='black',range=c(0,1.8),addscale = T)

par(fig=c(0.5,0.9,0.05,0.1), new=T)
plotBed(beddata = snpbed, chrom = paste0('chr',CHR), chromstart = START, chromend = END)
labelgenome(CHR,START,END,n=5,scale="Mb")

dev.off()

