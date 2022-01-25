# By: Golnaz Vahedi, PhD
# Date: April 25th, 2019
# Figure 3: the goal here is to find super-enhancers from K27ac peaks merged by 12.5kb following Rick Young's approach
# merge_K27ac_peaks.sh combined two technical replicates and merged the final peaks- I then counted tags at merged peaks and used MapSEs_YoungMethod.R to calculate super and typical enhancers

rm(list=ls())
library('GenomicRanges')
library('gplots')
library('pheatmap')
library("reshape")
library("plyr")
library('edgeR')
library('ggplot2')
#*****************************************************************
working_folder='/mnt/data0/sora/stripe/human_mouse_conserve/Super-enhancers/'
setwd(working_folder)
#*****************************************************************
# To define genes in Granges
#*****************************************************************
TSS_minus = 1000
TSS_plus = 1000
# Promoter Coordinates in RefSeq
RefSeq = read.table(file.path('/mnt/data0/golnaz/references/mouse/mm10_RefSeq.bed'),
                    header=F,stringsAsFactors=F,sep='\t')
RefSeq$V1 = gsub("chr","",RefSeq$V1)
RefSeq_gr <- GRanges(seqnames= Rle(RefSeq[,1]),ranges=IRanges(RefSeq[,2],RefSeq[,3]),strand=RefSeq[,6],intens=RefSeq[,4])

RefSeqPos1_gr <- RefSeq_gr[strand(RefSeq_gr)=='+',]
RefSeqPos2_gr <- GRanges(seqnames= seqnames(RefSeqPos1_gr),
                         ranges=IRanges(start(RefSeqPos1_gr)-TSS_minus,
                                        start(RefSeqPos1_gr)+TSS_plus),strand=strand(RefSeqPos1_gr),
                         intens=elementMetadata(RefSeqPos1_gr)[,"intens"])
# 
RefSeqNeg1_gr <- RefSeq_gr[strand(RefSeq_gr)=='-',]
RefSeqNeg2_gr <- GRanges(seqnames= seqnames(RefSeqNeg1_gr),
                         ranges=IRanges(end(RefSeqNeg1_gr)-TSS_plus,end(RefSeqNeg1_gr)+TSS_minus),
                         strand=strand(RefSeqNeg1_gr),
                         intens=elementMetadata(RefSeqNeg1_gr)[,"intens"])

RefSeq_excludePromoter <- c(RefSeqPos2_gr,RefSeqNeg2_gr)
#******************************************************************************************************************
#This function calculates the cutoff by sliding a diagonal line and finding where it is tangential (or as close as possible)
#******************************************************************************************************************
calculate_cutoff <- function(inputVector, drawPlot,save_file,...){
  inputVector <- sort(inputVector)
  inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
  slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
  xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum)
  #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
  y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
  
  if(drawPlot){  #if TRUE, draw the plot
    pdf(save_file)
    plot(1:length(inputVector), inputVector,type="p",...)
    b <- y_cutoff-(slope* xPt)
    abline(v= xPt,h= y_cutoff,lty=2,col=8)
    #lines(xPt,y_cutoff,col=1)
    points(xPt,y_cutoff,pch=16,cex=0.1,col=2)
    abline(coef=c(b,slope),col=2)
    title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
    axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
    dev.off()
  }
  return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector)))
}

#******************************************************************************************************************
#this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
#******************************************************************************************************************
numPts_below_line <- function(myVector,slope,x){
  yPt <- myVector[x]
  b <- yPt-(slope*x)
  xPts <- 1:length(myVector)
  return(sum(myVector<=(xPts*slope+b)))
}

#******************************************************************************************************************
# major function
#******************************************************************************************************************
find_SEs_fromMerged <- function(merged_file)
#function(working_folder, filedir, merged_file)
{#for outputs
  output_filename=strsplit(merged_file,"_merged.bedgraph")
  SEdir=file.path(working_folder,"SE")
  TEdir=file.path(working_folder,"TE")
  
  merged_data = read.table(file.path(filedir,merged_file),
                           header=F,stringsAsFactors=F,sep='\t')
  
  All_Enhancers_gr <- sort(GRanges(seqnames= Rle(merged_data[,1]),
                                    ranges=IRanges(merged_data[,2],merged_data[,3])))
  
# merged_data[,4]=log10(merged_data[,4])
  merged_data[,4] = (merged_data[,4])/max(merged_data[,4])
  All_Enhancers_gr$metadata = merged_data[,4]
  sorted_merged_data = merged_data[(sort(merged_data[,4],index.return = T))$ix,]
  
  inputVector = sorted_merged_data[,4]
  cutoff_options <- calculate_cutoff(inputVector, 
                                     drawPlot='TRUE',
                                     paste(output_filename,'_SEplot.pdf',sep=''),
                                    # xlab=paste(rankBy_factor,'_enhancers'),
                                     #ylab=paste(rankByFactor,' Signal','- ',wceName),
                                     lwd=2,col=4)
  superEnhancerRows <- which(inputVector > cutoff_options$absolute)
  # find a1 is the question for super enhancers
  
  OBS <- sorted_merged_data
  SE_merged_data = OBS[OBS[,4]> cutoff_options$absolute,]
  SE_merged_data_gr <- sort(GRanges(seqnames= Rle(SE_merged_data[,1]),
                                    ranges=IRanges(SE_merged_data[,2],SE_merged_data[,3]),
                                    intens=SE_merged_data[,4]))
  
  temp = countOverlaps(All_Enhancers_gr,SE_merged_data_gr)
  TE_gr = All_Enhancers_gr[temp==0,]
  
  ListofPeaks <- TE_gr
  Intermediate <- countOverlaps(ListofPeaks,RefSeq_excludePromoter)
  ListofPeaks_nopromoters <- ListofPeaks[Intermediate >= 0, ]
  mytable = data.frame(as.vector(seqnames(ListofPeaks_nopromoters)),
                       start(ListofPeaks_nopromoters), end(ListofPeaks_nopromoters), elementMetadata(ListofPeaks_nopromoters))
                           
  write.table(mytable, file.path(TEdir, paste(output_filename,"_TE.bedgraph",sep="")),
                                       row.names=F,quote=F,sep='\t',col.names=F)
                           
  mytable = data.frame(as.vector(seqnames(SE_merged_data_gr)),
                       start(SE_merged_data_gr), end(SE_merged_data_gr),
                       elementMetadata(SE_merged_data_gr)[,"intens"])
  
  write.table(mytable, file.path(SEdir, paste(output_filename,"_SE.bedgraph",sep="")),
              row.names=F,quote=F,sep='\t',col.names=F)
  
  
  return((SE_merged_data_gr))
}

#-----------------------------------------------------------------------------
# find all TEs and SEs
#-----------------------------------------------------------------------------
#system('sh merge_K27ac_peaks.sh ')
mydir=file.path(working_folder,"SE")
filedir=file.path(working_folder,"TempFiles")
setwd(mydir)

myfiles=dir(filedir)
merged_files=myfiles[grep("_merged.bedgraph",myfiles)]
out_GR = lapply(merged_files,find_SEs_fromMerged)

#-----------------------------------------------------------------------------
# Intersection with stripes
#-----------------------------------------------------------------------------
stripe='/mnt/data0/sora/stripe/stripenn_out/DP_T_5kb_VCSQRT3/result_filtered.txt'
#loop=/
se='/mnt/data0/sora/stripe/human_mouse_conserve/Super-enhancers/SE/ChIP-seq_K27ac_DP_BL6_Rep1_36160142_S1.sorted_SE.bedgraph'
te='/mnt/data0/sora/stripe/human_mouse_conserve/Super-enhancers/TE/ChIP-seq_K27ac_DP_BL6_Rep1_36160142_S1.sorted_TE.bedgraph'

stripe = read.delim(stripe,header=T,sep='\t')
stripe = stripe[which(stripe$Stripiness>0 & stripe$pvalue<0.05),]
se = read.delim(se, header=F)
te = read.delim(te,header=F)

colnames(stripe)[4:6] = colnames(se) = colnames(te) = c('chr','start','end')
colnames(se)[4] = colnames(te)[4] = 'elementMetadata'
A = findOverlaps(makeGRangesFromDataFrame(se,keep.extra.columns = T), makeGRangesFromDataFrame(stripe[,4:6]))
idx_se= unique(A@from)
B = findOverlaps(makeGRangesFromDataFrame(te,keep.extra.columns = T), makeGRangesFromDataFrame(stripe[,4:6]))
idx_te = unique(B@from)

intens = data.frame(intens=c(se$elementMetadata, te$elementMetadata), stripe=1)
intens$stripe[c(idx_se, (idx_te+nrow(se)))] = 2
intens = intens[order(intens$intens,decreasing = T),]

plot.new()
par(fig=c(0.1,0.9,0.1,0.9))
plot(x=1:nrow(intens),y=intens$intens, col=intens$stripe, pch=20,cex=2, ylab='Super-enhancer intensity',xlab='Index',las=1)
