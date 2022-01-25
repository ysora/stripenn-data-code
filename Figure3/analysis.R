library(GenomicRanges)

ch = read.delim('/mnt/data0/sora/stripe/CH12_DN3_compare/union_stripes_score_CH12.tsv')
dn = read.delim('/mnt/data0/sora/stripe/CH12_DN3_compare/union_stripes_score_DN3EV.tsv')

stripiness_ch = ch$Stripiness_added
stripiness_dn = dn$Stripiness_added

#tsum_ch = ch$O_Sum_added
#tsum_dn = dn$O_Sum_added

tsum_ch = ch$O.E_Total_added
tsum_dn = dn$O.E_Total_added

plot(stripiness_ch, stripiness_dn, xlab="CH12", ylab="DN3", xlim=c(-25,35),ylim=c(-25,35), col=rgb(0,0,0,0.5), pch=20,cex=0.7,las=1)
cor.test(stripiness_ch, stripiness_dn)
x=y=c(0,1)
abline(lm(y~x), lwd=2,col='gray80',lty=2)
abline(v=c(0,2),h=c(0,2), lwd=2,col=c('gray80','red',"gray80","red"),lty=2)

plot(tsum_ch, tsum_dn, xlab="CH12", ylab="DN3", col=rgb(0,0,0,0.5), pch=20)
cor.test(tsum_ch, tsum_dn)

# sort by difference
union = data.frame(ch, dn[,13:18])
union = union[order(abs(stripiness_ch - stripiness_dn), decreasing = T),]
head(union)


# Gene expression

# Biomart
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version = 100)
datasets = listDatasets(ensembl)
ensembl = useDataset(dataset = 'mmusculus_gene_ensembl', mart = ensembl)
bm = getBM(attributes = c('ensembl_gene_id', 'external_gene_name' ,'entrezgene_id','transcript_length'), mart = ensembl)
bm2 = bm %>%
  group_by(ensembl_gene_id) %>%
  summarise_each(max)

# B cell rna-seq 
rna_ch1 = read.delim('/mnt/alvand/sora/public_data/CH12/RNAseq/04.raw_count/4DNFIUULMHVG.tsv', header=F, skip=4)
rna_ch2 = read.delim('/mnt/alvand/sora/public_data/CH12/RNAseq/04.raw_count/4DNFIAQ3ZY3Z.tsv', header=F, skip=4)
rna_ch3 = read.delim('/mnt/alvand/sora/public_data/CH12/RNAseq/04.raw_count/4DNFIUL67HUG.tsv', header=F, skip=4)
rna_ch4 = read.delim('/mnt/alvand/sora/public_data/CH12/RNAseq/04.raw_count/4DNFIQTM1BBV.tsv', header=F, skip=4)

raw_ch = data.frame(ch1=rna_ch1$V2,ch2=rna_ch2$V2,ch3=rna_ch3$V2,ch4=rna_ch4$V2)
head(raw_ch)
rownames(raw_ch) = rna_ch1$V1
raw_ch$ensembl_gene_id = rownames(raw_ch)
raw_ch = merge(raw_ch, bm2, by="ensembl_gene_id")

rna_dn1 = read.delim('/mnt/alvand/VahediLab/GenomicsData/Mouse/T_lymphocytes/RNAseq/DN3_Scid.adh/04.raw_count/RNA_DN3_EmptyVector_2_FC75050_S6.tsv', header=F, skip=4)
rna_dn2 = read.delim('/mnt/alvand/VahediLab/GenomicsData/Mouse/T_lymphocytes/RNAseq/DN3_Scid.adh/04.raw_count/RNA_DN3_EmptyVector_1_FC75050_S5.tsv', header=F, skip=4)
rna_dn = data.frame(dn1=rna_dn1$V2, dn2=rna_dn2$V2)
rna_dn$ensembl_gene_id = rna_dn1$V1

ch_dn_rna = merge(raw_ch, rna_dn, by='ensembl_gene_id')

ch_dn_rna[which(ch_dn_rna$external_gene_name == "Ets1"),]


# TPM 
tpm = function(rawcount, transcript_length)
{
  tl = transcript_length/1000
  tc = rawcount / tl
  sumv = sum(tc) / 1000000
  tc = tc / sumv
  return(tc)
}

ch1 = tpm(ch_dn_rna$ch1, ch_dn_rna$transcript_length)
ch2 = tpm(ch_dn_rna$ch2, ch_dn_rna$transcript_length)
ch3 = tpm(ch_dn_rna$ch3, ch_dn_rna$transcript_length)
ch4 = tpm(ch_dn_rna$ch4, ch_dn_rna$transcript_length)

dn1 = tpm(ch_dn_rna$dn1, ch_dn_rna$transcript_length)
dn2 = tpm(ch_dn_rna$dn2, ch_dn_rna$transcript_length)

ch = data.frame(ch1,ch2,ch3,ch4)
dn = data.frame(dn1,dn2)
rownames(ch) = rownames(dn) = ch_dn_rna$ensembl_gene_id
ch$gene = dn$gene = ch_dn_rna$external_gene_name
ch_avg = apply(ch[,1:4],1,mean)
dn_avg = apply(dn[,1:2],1,mean)
names(ch_avg) = names(dn_avg) = ch_dn_rna$external_gene_name

# Test
ch[which(ch_dn_rna$external_gene_name == 'Ebf1'),]
dn[which(ch_dn_rna$external_gene_name == 'Ebf1'),]

# stripe and gene expression 
##  (1): The most highly expressed gene in each stripe

ref = read.delim('/mnt/data0/golnaz/references/mouse/mm10_RefSeq.bed', header=F)
colnames(ref)[1:4] = c('chr','start','end','gene')
union_df = union[,4:6]
colnames(union_df) = c('chr','start','end')
union_df$chr = paste0("chr",union_df$chr)

ov = findOverlaps(makeGRangesFromDataFrame(union_df), makeGRangesFromDataFrame(ref, keep.extra.columns = T))

HEG_B = c()
HEG_T = c()
HEG_all=c()
HEG_B_cnt = c()
HEG_T_cnt = c()
HEG_all_cnt=c()
MEAN_B = c()
MED_B = c()
MEAN_T = c()
MED_T = c()

for(i in 1:nrow(union_df))
{
  idx = ov@to[which(ov@from == i)]
  idx = unique(idx)
  gene = ref$gene[idx]
  gene = intersect(gene, names(ch_avg))
  if(length(gene)>0){
    b = ch_avg[gene]
    t = dn_avg[gene]
    HEG_B = append(HEG_B, names(b)[which.max(b)])
    HEG_T = append(HEG_T, names(t)[which.max(t)])
    HEG_all = append(HEG_all, names(c(b,t))[which.max(c(b,t))])
    HEG_B_cnt = append(HEG_B_cnt, b[which.max(b)])
    HEG_T_cnt = append(HEG_T_cnt, t[which.max(t)])  
    HEG_all_cnt = append(HEG_all_cnt, c(b,t)[which.max(c(b,t))])
    MEAN_T = append(MEAN_T, mean(t))
    MEAN_B = append(MEAN_B, mean(b))
    MED_T = append(MED_T, median(t))
    MED_B = append(MED_B, median(b))
  }else{
    HEG_B = append(HEG_B, NA)
    HEG_T = append(HEG_T, NA)
    HEG_all = append(HEG_all, NA)
    HEG_B_cnt = append(HEG_B_cnt, NA)
    HEG_T_cnt = append(HEG_T_cnt, NA)  
    HEG_all_cnt = append(HEG_all_cnt, NA)
    MEAN_T = append(MEAN_T, NA)
    MEAN_B = append(MEAN_B, NA)
    MED_T = append(MED_T, NA)
    MED_B = append(MED_B, NA)
    
  }
}

union$HEG_B = HEG_B
union$HEG_T = HEG_T
union$HEG_all = HEG_all
union$HEG_B_cnt = HEG_B_cnt
union$HEG_T_cnt = HEG_T_cnt
union$HEG_all_cnt = HEG_all_cnt

head(union)
union$HEG_all_B_cnt = ch_avg[union$HEG_all]
union$HEG_all_T_cnt = dn_avg[union$HEG_all]

idx1 = which(union$Stripiness_added>2 & union$Stripiness_added.1 < 0)
idx2 = which(union$Stripiness_added<0 & union$Stripiness_added.1 > 2)
idx = sort(c(idx1,idx2))
FC_Total = log2((union$O.E_Total_added[idx]+1) / (union$O.E_Total_added.1[idx]+1)) # Bcell/Tcell
SUB_stri = union$Stripiness_added[idx] - union$Stripiness_added.1[idx]
FC_GE    = (union$HEG_all_B_cnt[idx]+1) / (union$HEG_all_T_cnt[idx]+1)
FC_GE_MEAN=(MEAN_B[idx]+1) / (MEAN_T[idx]+1)
FC_GE_MED =(MED_B[idx]+1) / (MED_T[idx]+1)

cor.test(FC_Total, FC_GE, na.rm=T)
cor.test(FC_Total, FC_GE_MEAN, na.rm=T)
cor.test(FC_Total, FC_GE_MED, na.rm=T)

cor.test(SUB_stri, FC_GE, na.rm=T)
cor.test(SUB_stri, FC_GE_MEAN, na.rm=T)
cor.test(SUB_stri, FC_GE_MED, na.rm=T)

plot(SUB_stri, log2(FC_GE), xlab="Stripiness change (B cell - T cell)", ylab="Log2 of gene expression fold change (B cell / T cell)", cex=0.8, pch=20)
abline(lm(log2(FC_GE) ~ SUB_stri), col='red', lty=2,lwd=2)

pdf('/mnt/data0/sora/stripe/CH12_DN3_compare/HEG_FC_GE_FC_TSUM.pdf', width = 7, height = 7)
plot(FC_Total, log2(FC_GE), xlab="Fold change of total sum (B cell / T cell)", ylab="Log2 of gene expression fold change (B cell / T cell)", cex=0.8, pch=20, main="The most highly expressed genes")
abline(lm(log2(FC_GE) ~ FC_Total), col='red', lty=2,lwd=2)
text(x=1, y=-5, labels = paste0("R=",round(cor.test(FC_GE,FC_Total)$estimate,3), "\nP=",signif(cor.test(FC_GE,FC_Total)$p.value,3)))
dev.off()

pdf('/mnt/data0/sora/stripe/CH12_DN3_compare/MED_FC_GE_FC_TSUM.pdf', width = 7, height = 7)
plot(FC_Total, log2(FC_GE_MED),xlab="Fold change of total sum (B cell / T cell)", ylab="Log2 of median gene expression fold change (B cell / T cell)", cex=0.8, pch=20, main="Median gene expression")
abline(lm(log2(FC_GE_MED) ~ FC_Total), col='red', lty=2,lwd=2)
text(x=1, y=-5, labels = paste0("R=",round(cor.test(log2(FC_GE_MED),FC_Total)$estimate,3), "\nP=",signif(cor.test(log2(FC_GE_MED),FC_Total)$p.value,3)))
dev.off()

pdf('/mnt/data0/sora/stripe/CH12_DN3_compare/MEAN_FC_GE_FC_TSUM.pdf', width = 7, height = 7)
plot(FC_Total, log2(FC_GE_MEAN),xlab="Fold change of total sum (B cell / T cell)", ylab="Log2 of mean gene expression fold change (B cell / T cell)", cex=0.8, pch=20, main="Mean gene expression")
abline(lm(log2(FC_GE_MEAN) ~ FC_Total), col='red', lty=2,lwd=2)
text(x=1, y=-5, labels = paste0("R=",round(cor.test(log2(FC_GE_MEAN),FC_Total)$estimate,3), "\nP=",signif(cor.test(log2(FC_GE_MEAN),FC_Total)$p.value,3)))
dev.off()

pdf('/mnt/data0/sora/stripe/CH12_DN3_compare/MED_FC_GE_SUB_STRIPINESS.pdf', width = 7, height = 7)
plot(SUB_stri, log2(FC_GE_MED), xlab="Stripiness change (B cell - T cell)", ylab="Log2 of gene expression fold change (B cell / T cell)", main="Median gene expression", cex=0.8, pch=20)
abline(lm(log2(FC_GE_MED) ~ SUB_stri), col='red', lty=2,lwd=2)
text(x=1, y=-5, labels = paste0("R=",round(cor.test(log2(FC_GE_MED),SUB_stri)$estimate,3), "\nP=",signif(cor.test(log2(FC_GE_MED),SUB_stri)$p.value,3)))
dev.off()

pdf('/mnt/data0/sora/stripe/CH12_DN3_compare/MEAN_FC_GE_SUB_STRIPINESS.pdf', width = 7, height = 7)
plot(SUB_stri, log2(FC_GE_MEAN), xlab="Stripiness change (B cell - T cell)", ylab="Log2 of gene expression fold change (B cell / T cell)", main="Mean gene expression", cex=0.8, pch=20)
abline(lm(log2(FC_GE_MEAN) ~ SUB_stri), col='red', lty=2,lwd=2)
text(x=1, y=-5, labels = paste0("R=",round(cor.test(log2(FC_GE_MEAN),SUB_stri)$estimate,3), "\nP=",signif(cor.test(log2(FC_GE_MEAN),SUB_stri)$p.value,3)))
dev.off()

plot(FC_Total, log2(FC_GE_MEAN))
abline(lm(log2(FC_GE_MEAN) ~ FC_Total), col='red', lty=2,lwd=2)

cor.test(FC_GE, FC_GE_MED)
plot(FC_GE, FC_GE_MEAN,ylim=c(0,100),xlim=c(0,100))
plot(FC_GE, FC_GE_MEAN,ylim=c(0,100),xlim=c(0,100))


## Candidates to draw
idx3 =  which(SUB_stri > 10 & log2(FC_GE)>4)
idx3 = idx3[1]

union[idx[idx3],]


idx4 =  which(SUB_stri < -15 & log2(FC_GE)< -5)
idx4 = idx4[1]
union[idx[idx4],]

## Sushi image --- B cell specific
library(Sushi)

CHR = 1  
START = 106515001-50000
END = 106715000 +50000
coord = paste0(CHR,":",START,"-",END)
cool2="/mnt/data0/sora/stripe/hic/DN3EV.allValidPairs.mcool"
cool1="/mnt/data0/sora/stripe/hic/CH12_4DNFI8KBXYNL.mcool"
cool = paste(cool1, cool2, sep=",")

BB='/mnt/data0/wenliang/software/UCSC/linux.x86_64/bigWigToBedGraph'
ch_bw = '/mnt/alvand/sora/public_data/CH12/RNAseq/03.bigwig/4DNFIAQ3ZY3Z.str1.bw'
dn_bw = '/mnt/alvand/VahediLab/GenomicsData/Mouse/T_lymphocytes/RNAseq/DN3_Scid.adh/03.bigwig/RNA_DN3_EmptyVector_1_FC75050_S5.str1.bw'
out = "/mnt/data0/sora/stripe/hicPlotter/Bcl2_"

system(paste("mkdir ",out))

system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq.py --cool ",cool," --coord ",coord," --norm VC_SQRT --out ",out," --res 5000"))
system(paste0(BB ," ",ch_bw, " ",out,"CH12.bedgraph " ,"-chrom=",CHR," -start=",START, " -end=",END))
system(paste0(BB ," ",dn_bw, " ",out,"DN3.bedgraph " ,"-chrom=",CHR," -start=",START, " -end=",END))

ch_bg = read.delim(paste0(out,"CH12.bedgraph"),header=F)
dn_bg = read.delim(paste0(out,"DN3.bedgraph"),header=F)

EV = read.delim(paste0(out,'0.txt'),header=T,row.names=1)
KO = read.delim(paste0(out,'1.txt'),header=T,row.names=1)
colnames(EV) = colnames(KO) = row.names(EV)

pdf(paste0(out,'VC_SQRT.pdf'), width = 10, height = 3)
par(cex=0.7, mai=c(.1,0.1,0.2,0.1))
par(fig=c(0.1,0.35,0.3,0.9))
phic1 = plotHic(EV, CHR,START,END, max_y = 30,
                zrange=c(0,13), palette=colorRampPalette(c("white",'red'))
)

par(fig=c(0.36,0.6,0.3,0.9),new=T)
phic1 = plotHic(KO, CHR,START,END, max_y = 30,
                zrange=c(0,13), palette=colorRampPalette(c("white",'red'))
)

par(fig=c(0.1,0.35,0.1,0.3), new=T)
colnames(ref)[4] = 'gene'
plotGenes(ref[which(ref$chr==paste0('chr',CHR) &ref$start>=START & ref$end<=END),],CHR,START,END ,
          height=0.5,plotgenetype="arrow",bentline=F,col="#465075", labeltext = TRUE, type='box',
          labeloffset= 0.4,fontsize=1.5, arrowlength = 0.01, bheight = 0.1, packrow = TRUE)

par(fig=c(0.36,0.6,0.1,0.3), new=T)
colnames(ref)[4] = 'gene'
plotGenes(ref[which(ref$chr==paste0('chr',CHR) &ref$start>=START & ref$end<=END),],CHR,START,END ,
          height=0.5,plotgenetype="arrow",col="#465075", labeltext = TRUE, type='box',
          labeloffset= 0.4,fontsize=1.5, arrowlength = 0.01, bheight = 0.1, packrow = TRUE)

par(fig=c(0.61,0.9,0.1,0.9), new=T)
boxplot(list(CH12=unlist(ch[which(ch$gene=='Bcl2'),1:4]), DN3=unlist(dn[which(dn$gene=='Bcl2'),1:2])))


dev.off()



## Sushi image --- T cell specific
library(Sushi)

CHR = 12
START = 106585001-100000
END = 108005000 +100000
coord = paste0(CHR,":",START,"-",END)
cool2="/mnt/data0/sora/stripe/hic/DN3EV.allValidPairs.mcool"
cool1="/mnt/data0/sora/stripe/hic/CH12_4DNFI8KBXYNL.mcool"
cool = paste(cool1, cool2, sep=",")


BB='/mnt/data0/wenliang/software/UCSC/linux.x86_64/bigWigToBedGraph'
ch_bw = '/mnt/alvand/sora/public_data/CH12/RNAseq/03.bigwig/4DNFIAQ3ZY3Z.str1.bw'
dn_bw = '/mnt/alvand/VahediLab/GenomicsData/Mouse/T_lymphocytes/RNAseq/DN3_Scid.adh/03.bigwig/RNA_DN3_EmptyVector_1_FC75050_S5.str1.bw'
out = "/mnt/data0/sora/stripe/hicPlotter/Bcl11b_"

system(paste("mkdir ",out))

system(paste0("/mnt/data0/apps/anaconda/Anaconda2-5.2/envs/py38/bin/python3 /mnt/data0/sora/bin/py/extract_contact_freq.py --cool ",cool," --coord ",coord," --norm VC_SQRT --out ",out," --res 5000"))
system(paste0(BB ," ",ch_bw, " ",out,"CH12.bedgraph " ,"-chrom=",CHR," -start=",START, " -end=",END))
system(paste0(BB ," ",dn_bw, " ",out,"DN3.bedgraph " ,"-chrom=",CHR," -start=",START, " -end=",END))

ch_bg = read.delim(paste0(out,"CH12.bedgraph"),header=F)
dn_bg = read.delim(paste0(out,"DN3.bedgraph"),header=F)

EV = read.delim(paste0(out,'0.txt'),header=T,row.names=1)
KO = read.delim(paste0(out,'1.txt'),header=T,row.names=1)
colnames(EV) = colnames(KO) = row.names(EV)

pdf(paste0(out,'VC_SQRT.pdf'), width = 10, height = 3)
par(cex=0.7, mai=c(.1,0.1,0.2,0.1))
par(fig=c(0.1,0.35,0.3,0.9))
phic1 = plotHic(EV, CHR,START,END, max_y = 165,
                zrange=c(0,5), palette=colorRampPalette(c("white",'red'))
)

par(fig=c(0.36,0.6,0.3,0.9),new=T)
phic1 = plotHic(KO, CHR,START,END, max_y = 165,
                zrange=c(0,5), palette=colorRampPalette(c("white",'red'))
)

par(fig=c(0.1,0.35,0.1,0.3), new=T)
colnames(ref)[4] = 'gene'
plotGenes(ref[which(ref$chr==paste0('chr',CHR) &ref$start>=START & ref$end<=END),],CHR,START,END ,
          height=0.5,plotgenetype="arrow",bentline=F,col="#465075", labeltext = TRUE, type='box',
          labeloffset= 0.4,fontsize=1.5, arrowlength = 0.01, bheight = 0.1, packrow = TRUE)

par(fig=c(0.36,0.6,0.1,0.3), new=T)
colnames(ref)[4] = 'gene'
plotGenes(ref[which(ref$chr==paste0('chr',CHR) &ref$start>=START & ref$end<=END),],CHR,START,END ,
          height=0.5,plotgenetype="arrow",col="#465075", labeltext = TRUE, type='box',
          labeloffset= 0.4,fontsize=1.5, arrowlength = 0.01, bheight = 0.1, packrow = TRUE)

par(fig=c(0.61,0.9,0.1,0.9), new=T)
boxplot(list(CH12=unlist(ch[which(ch$gene=='Bcl11b'),1:4]), DN3=unlist(dn[which(dn$gene=='Bcl11b'),1:2])))


dev.off()



