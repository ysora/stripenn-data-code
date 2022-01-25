
### ENHANCER MARKS TEST ###
library(GenomicRanges)
library(rtracklayer)

bw=import('/mnt/data0/sora/stripe/ChIP_hiC/DPT/H3K27ac/05.bigwig/merged.raw.bw',format='bigwig')

stripy_top1=read.delim('/mnt/data0/sora/stripe/GeneExpression/210702_top1_genes_coordinate_interaction.txt')
stripy_oth=read.delim('/mnt/data0/sora/stripe/GeneExpression/210702_other_genes_coordinate_interaction.txt')

str_top_genes = stripy_top1$name
str_oth_genes = stripy_oth$name

nostr_top_genes = top1$name
nostr_oth_genes = oth$name

getEnhancerScore = function(bw, gene)
{
  idx = which(promoter$V4 %in% gene)
  genepromchr = promoter$V1[idx]
  genepromsta = promoter$V2[idx]
  genepromend = promoter$V3[idx]
  df = data.frame(chr=genepromchr, start=genepromsta-1000, end = genepromend+1000)
  ov = findOverlaps(makeGRangesFromDataFrame(df), bw)
  
  mp= c()
  for(g in gene)
  {
    idx2 = which(promoter$V4 == g)
    idx3 = which(idx %in% idx2)
    maxPeak = max(bw$score[ov@to[which(ov@from %in% idx3)]])
    mp = append(mp,maxPeak)
  }
  
  return(mp)
}

str_top_genes_peak =getEnhancerScore(bw=bw, str_top_genes)
str_oth_genes_peak = getEnhancerScore(bw, str_oth_genes)
nostr_top_genes_peak = getEnhancerScore(bw, nostr_top_genes)
nostr_oth_genes_peak = getEnhancerScore(bw, nostr_oth_genes)
boxplot(list(stripy_top1 = str_top_genes_peak, nostripy_top1 = nostr_top_genes_peak, stripy_others=str_oth_genes_peak,nostripy_others=nostr_oth_genes_peak), main = 'H3K27ac bias to Top1 gene', ylab="Maximum H3K27ac peak")

wilcox.test(str_top_genes_peak, nostr_top_genes_peak)
