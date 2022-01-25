options = commandArgs(trailingOnly = TRUE)
# 
# options = c('/home/sora/Desktop/sora_alborz/sora/stripe/homer/PCA/Vian_aB_30hrs/Vian_30hrs_PCA.PC1.txt',
#             '/home/sora/Desktop/sora_alborz/sora/stripe/stripenn_out/Vian_aB_30hrs_10kb_VCSQRT/result_unfiltered.txt',
#             '0.05,0.1',
#             '-1000,0', 
#             '/home/sora/Desktop/sora_alborz/sora/stripe/homer/PCA/Vian_aB_30hrs/','Vian_aB_30hrs_10kb_VCSQRT')
pca = options[1]
stripe = options[2]
pval=options[3]
scut=options[4]
outdir=options[5]
outprefix=options[6]
tadfile=options[7]


library(GenomicRanges)
library(ggplot2)
library(scales)
# 
# stripe = read.delim('/home/sora/Desktop/sora_alborz/sora/stripe/result/DPT_20200714_2/stripes.txt',header=T)
# pval = '0.2,0.1,0.05'
# scut = '-1000,0'
# outdir = '/home/sora/Desktop/sora_alborz/sora/stripe/homer/PCA/BL6DPT/'
# outprefix = 'BL6_10kb_previous'
# pca = read.delim('/home/sora/Desktop/sora_alborz/sora/stripe/homer/PCA/BL6DPT/DP_T_50K.PC1.txt')
# 

pca = read.delim(pca)
pca_df = makeGRangesFromDataFrame(pca)
stripe = read.delim(stripe,header=T)
pval = as.numeric(unlist(strsplit(pval,split = ',',fixed=T)))
scut = as.numeric(unlist(strsplit(scut,split = ',',fixed=T)))

test = substr(stripe$chr[1],1,3) == 'chr'
if(!test)
{
  stripe$chr = paste0('chr',stripe$chr)
  stripe$chr2 = stripe$chr
}


for(p in pval)
{
  for(s in scut)
  {
    dat = stripe[which(stripe$pvalue<p & stripe$Stripiness>s),]
    dat_df = makeGRangesFromDataFrame(data.frame(chr=dat$chr,start=dat$pos3,end=dat$pos4))
    ud = c()
    for(i in 1:nrow(dat))
    {
      if(dat$pos1[i] == dat$pos3[i]){ud = append(ud,2)}
      if(dat$pos2[i] == dat$pos4[i]){ud = append(ud,1)}
    }
    
    ov = findOverlaps(dat_df, pca_df)
    n_stripe = unique(ov@from)
    PROP = array(NA,nrow(dat))
    for(i in n_stripe)
    {
      pcaidx = which(ov@from == i)
      pc1 = pca$PC1[ov@to[pcaidx]]
      prop = length(which(pc1>0))/length(pc1)
      PROP[i] = prop
    }
    dat$PROP = PROP
    gg_b = ggplot_build(
      ggplot() + geom_histogram(aes(x = PROP), binwidth=0.1)
    )
    nu_bins = dim(gg_b$data[[1]])[1]
    
    pdf(paste0(outdir,'/',outprefix,'_p_',p,'_s_',s,'.pdf'))
    P = ggplot(dat, aes(PROP,fill=..x..)) + geom_histogram(aes(x=PROP,y = ..count../length(PROP)),group=1, binwidth=0.1)
    mid = 0.5
    P2=P + scale_fill_gradient2(low='blue',high='red',midpoint=mid,mid='white')+
      annotate(geom="text",x=0.25,y=0.5,label=paste0('Only A=',round(length(which(dat$PROP==1.0))/length(dat$PROP),3)))+
      annotate(geom="text",x=0.25,y=0.3,label=paste0('Only B=',round(length(which(dat$PROP==0.0))/length(dat$PROP),3)))+
      xlab('Proportion of A compartment') + 
      ylab('Frequency')
    print(P2)
    dev.off()
  }
}


# TAD
dat = read.delim(tadfile,header=F)
dat_df = makeGRangesFromDataFrame(data.frame(chr=dat$V1,start=dat$V2,end=dat$V3))

ov = findOverlaps(dat_df, pca_df)
n_stripe = unique(ov@from)
PROP = array(NA,nrow(dat))
for(i in n_stripe)
{
  pcaidx = which(ov@from == i)
  pc1 = pca$PC1[ov@to[pcaidx]]
  prop = length(which(pc1>0))/length(pc1)
  PROP[i] = prop
}
dat$PROP = PROP
gg_b = ggplot_build(
  ggplot() + geom_histogram(aes(x = PROP), binwidth=0.1)
)
nu_bins = dim(gg_b$data[[1]])[1]

pdf(paste0(outdir,'/TAD_ABcompartment.pdf'))
P = ggplot(dat, aes(PROP,fill=..x..)) + geom_histogram(aes(x=PROP,y = ..count../length(PROP)),group=1, binwidth=0.1)
mid = 0.5
P2=P + scale_fill_gradient2(low='blue',high='red',midpoint=mid,mid='white')+
  annotate(geom="text",x=0.25,y=0.5,label=paste0('Only A=',round(length(which(dat$PROP==1.0))/length(dat$PROP),3)))+
  annotate(geom="text",x=0.25,y=0.3,label=paste0('Only B=',round(length(which(dat$PROP==0.0))/length(dat$PROP),3)))+
  xlab('Proportion of A compartment') + 
  ylab('Frequency')+ylim(0,0.95)
print(P2)
dev.off()