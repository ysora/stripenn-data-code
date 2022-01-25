library(rtracklayer)
library(AnnotationHub)
library(Rsamtools)
library(gplots)
library(RColorBrewer)
library(ggpubr)
library(tibble)
library(viridis)
library(DataCombine)
library(tidyr)

# interpolation = function(scores, binsize)
# {
#   L1 = length(scores)
#   L2 = binsize
#   result = array(data = NA, L2)
#   for(i in 1:L2)
#   {
#     start = round((i-1)*L1/L2)+1
#     end  = round((i)*L1/L2)
#     result[i] = mean(scores[start:end])
#   }
#   return(result)
# }

print("AAA")

interpolation = function(df, binsize)
{
  # 2020.12.16 modified
  dif = df$ov_range_start[2:nrow(df)] - df$ov_range_end[1:(nrow(df)-1)]
  nonzero_idx = which(dif!=0)
  if(length(nonzero_idx)>0)
  {
    names_df = 1:nrow(df)
    names_df_update = names_df
    for(i in nonzero_idx)
    {
      names_df_update[(i+1):length(names_df_update)] =names_df_update[(i+1):length(names_df_update)]+1 
    }
    df2=data.frame(idx=names_df_update, df)
    adding = data.frame(idx=setdiff(1:names_df_update[length(names_df_update)], names_df_update), ov_range_start=df$ov_range_end[nonzero_idx], ov_range_end=df$ov_range_start[nonzero_idx + 1],ov_score=0)
    df2 = rbind(df2, adding)
    df2 = df2[order(df2$idx),]
    
    df = df2
  }
  idxNA = which(is.na(df$ov_range_start))
  if(length(idxNA)>0){df=df[-idxNA,]}
  
  width = df$ov_range_end - df$ov_range_start
  score = df$ov_score
  sum_width = sum(width)
  accum_width = cumsum(width)
  unit_size = floor(sum_width/binsize)
  
  result = array(data=NA, binsize)
  for(i in 1:binsize)
  {
    start = unit_size * (i-1) + 1
    end = unit_size * i
    start_idx = which(accum_width<=start)
    start_idx = start_idx[length(start_idx)]+1
    end_idx = which(accum_width<=end)
    if(length(start_idx)==0){start_idx = 1}
    if(length(end_idx) == 0){end_idx = 1}
    end_idx = end_idx[length(end_idx)] + 1
    if(end_idx > nrow(df)){end_idx = nrow(df)}

    if(start_idx == end_idx){
      result[i] = score[start_idx]
    }else{
      width_temp = width[start_idx:end_idx]
      score2 = score[start_idx:end_idx]
      if(start_idx==1){width_temp[1] = width_temp[1] - start + 1}
      if(start_idx!=1){width_temp[1] = length(c(start:accum_width[start_idx+1]))}
      width_temp[length(width_temp)] = length(accum_width[end_idx-1] : end )
      result[i] = sum(score2 * width_temp)/sum(width_temp)
    }
  }
  return(result)
}

options = commandArgs(trailingOnly = TRUE)
# options=c("/home/sora/Desktop/sora_alborz/sora/stripe/ChIP_hiC/DPT/H3K27ac/ChIP-seq_K27ac_DP_BL6_Rep2_36160142_S2.bw", '/home/sora/Desktop/sora_alborz/sora/stripe/stripenn_out/DP_T_5kb_VCSQRT/result_filtered_p0.1_S0_up.bed',
#         '/home/sora/Desktop/sora_alborz/sora/stripe/deeptools/2020_12_01_BL6_5kb_p_0.1_s0/H3K27ac_down.pdf','/home/sora/Documents/stripe/CTCF_DPT.png')
bw1 = import(options[1], format = 'BigWig')
bed = read.delim(options[2],header=F)
test = substr(bed$V1[1],1,3) == 'chr'
if(!test)
{
  bed$V1=paste0('chr',bed$V1)
}
a=10000
b=10000
m=3000

r1 = makeGRangesFromDataFrame(data.frame(chr=bed$V1, start = bed$V2, end = bed$V3))
ov = findOverlaps(bw1, r1)
ov_range = sort(ov@from)
ov_score = bw1$score[ov_range]

r1 = makeGRangesFromDataFrame(data.frame(chr=bed$V1, start = bed$V2, end = bed$V3))
ov = findOverlaps(bw1, r1)
ov_range = sort(ov@from)
ov_score = bw1$score[ov_range]

res = matrix(0,nrow=nrow(bed),ncol=300)
for(i in 1:nrow(bed))
{
  print(i)
  ov1 = ov[which(ov@to == i)]
  ov_range = sort(ov1@from)
  ov_range_start = bw1@ranges@start[ov_range]
  ov_range_end = ov_range_start + bw1@ranges@width[ov_range]
  ov_score = bw1$score[ov_range]
  ov_dat = data.frame(ov_range_start,ov_range_end,ov_score)
  start_site = bed$V2[i]
  if(ov_dat$ov_range_start[1] > start_site)
  {
    first_start_site = start_site
    first_end_site = ov_dat$ov_range_start[1]
    first_score = 0
    ov_dat=InsertRow(ov_dat,NewRow = c(first_start_site,first_end_site,first_score),RowNum = 1)
  }else{
    ov_dat$ov_range_start[1] = start_site
  }
  end_site = bed$V3[i]
  if(ov_dat$ov_range_end[nrow(ov_dat)] < end_site)
  {
    final_start_site = ov_dat$ov_range_end[nrow(ov_dat)]
    final_end_site = end_site
    final_score = 0
    ov_dat=InsertRow(ov_dat,NewRow = c(final_start_site,final_end_site,final_score),RowNum = nrow(ov_dat)+1)
  }else{
    ov_dat$ov_range_end[nrow(ov_dat)] = end_site
  }
  res[i,] = interpolation(ov_dat, 300)
}

res_b = matrix(0,nrow=nrow(bed),ncol=100)
r1 = makeGRangesFromDataFrame(data.frame(chr=bed$V1, start = bed$V2-b+1, end = bed$V2-1))
ov = findOverlaps(bw1, r1)
for(i in 1:nrow(bed))
{
  print(i)
  ov1 = ov[which(ov@to == i)]
  
  if(length(ov1) == 0)
  {
    res_b[i,] = 0
  }else{
    ov_range = sort(ov1@from)
    ov_range_start = bw1@ranges@start[ov_range]
    ov_range_end = ov_range_start + bw1@ranges@width[ov_range]
    ov_score = bw1$score[ov_range]
    
    
    ov_dat = data.frame(ov_range_start,ov_range_end,ov_score)
    start_site = bed$V2[i] - b + 1
    if(ov_dat$ov_range_start[1] > start_site)
    {
      first_start_site = start_site
      first_end_site = ov_dat$ov_range_start[1]
      first_score = 0
      ov_dat=InsertRow(ov_dat,NewRow = c(first_start_site,first_end_site,first_score),RowNum = 1)
    }else{
      ov_dat$ov_range_start[1] = start_site
    }
    end_site = bed$V2[i] -1
    if(ov_dat$ov_range_end[nrow(ov_dat)] < end_site)
    {
      final_start_site = ov_dat$ov_range_end[nrow(ov_dat)]
      final_end_site = end_site
      final_score = 0
      ov_dat=InsertRow(ov_dat,NewRow = c(final_start_site,final_end_site,final_score),RowNum = nrow(ov_dat)+1)
    }else{
      ov_dat$ov_range_end[nrow(ov_dat)] = end_site
    }
    res_b[i,] = interpolation(ov_dat, 100)
  }
}
res_a = matrix(0,nrow=nrow(bed),ncol=100)
r1 = makeGRangesFromDataFrame(data.frame(chr=bed$V1, start = bed$V3+1, end = bed$V3+a-1))
ov = findOverlaps(bw1, r1)
for(i in 1:nrow(bed))
{
  ov1 = ov[which(ov@to == i)]
  if(length(ov1) == 0)
  {
    res_a[i,] = 0
  }else{
    ov_range = sort(ov1@from)
    ov_range_start = bw1@ranges@start[ov_range]
    ov_range_end = ov_range_start + bw1@ranges@width[ov_range]
    ov_score = bw1$score[ov_range]
    ov_dat = data.frame(ov_range_start,ov_range_end,ov_score)
    start_site = bed$V3[i] +1
    if(ov_dat$ov_range_start[1] > start_site)
    {
      first_start_site = start_site
      first_end_site = ov_dat$ov_range_start[1]
      first_score = 0
      ov_dat=InsertRow(ov_dat,NewRow = c(first_start_site,first_end_site,first_score),RowNum = 1)
    }else{
      ov_dat$ov_range_start[1] = start_site
    }
    end_site = bed$V3[i] + a
    if(ov_dat$ov_range_end[nrow(ov_dat)] < end_site)
    {
      final_start_site = ov_dat$ov_range_end[nrow(ov_dat)]
      final_end_site = end_site
      final_score = 0
      ov_dat=InsertRow(ov_dat,NewRow = c(final_start_site,final_end_site,final_score),RowNum = nrow(ov_dat)+1)
    }else{
      ov_dat$ov_range_end[nrow(ov_dat)] = end_site
    }
    res_a[i,] = interpolation(ov_dat, 100)
    
  }
}
res_combined = cbind(res_b, res, res_a)

# res = matrix(0,nrow=nrow(bed),ncol=300)
# for(i in 1:nrow(bed))
# {
#   ov1 = ov[which(ov@to == i)]
#   if(length(ov1) == 0)
#   {
#     res[i,] = rep(0,300)
#   }else{
#     ov_range = sort(ov1@from)
#     ov_score = bw1$score[ov_range]  
#     res[i,] = interpolation(ov_score, 300)
#   }
# }
# 
# res_b = matrix(0,nrow=nrow(bed),ncol=100)
# r1 = makeGRangesFromDataFrame(data.frame(chr=bed$V1, start = bed$V2-b+1, end = bed$V2-1))
# ov = findOverlaps(bw1, r1)
# for(i in 1:nrow(bed))
# {
#   ov1 = ov[which(ov@to == i)]
#   if(length(ov1) == 0)
#   {
#     res_b[i,] = rep(0,100)
#   }else{
#     ov_range = sort(ov1@from)
#     ov_score = bw1$score[ov_range]  
#     res_b[i,] = interpolation(ov_score, 100)
#   }
# }
# res_a = matrix(0,nrow=nrow(bed),ncol=100)
# r1 = makeGRangesFromDataFrame(data.frame(chr=bed$V1, start = bed$V3+1, end = bed$V3+a-1))
# ov = findOverlaps(bw1, r1)
# for(i in 1:nrow(bed))
# {
#   ov1 = ov[which(ov@to == i)]
#   if(length(ov1) == 0){
#     res_a[i,] = rep(0,100)
#   }else{
#     ov_range = sort(ov1@from)
#     ov_score = bw1$score[ov_range]
#     res_a[i,] = interpolation(ov_score, 100)
#   }
# }
# res_combined = cbind(res_b, res, res_a)

qv = quantile(res_combined, 0.9, na.rm=T)
res_combined[which(res_combined>qv)] = qv
idxNA = which(is.na(apply(res_combined,1,sum)))
if(length(idxNA)>0){res_combined = res_combined[-idxNA,]}

class(res_combined)
colnames(res_combined) = rep("", ncol(res_combined))
rownames(res_combined) = rep("", nrow(res_combined))

rowRank = order(apply(res_combined,1,sum))
res_combined = res_combined[rowRank,]

colnames(res_combined) = c(c(-100:-1),c(1:400))
#colnames(res_combined) = c("-1kb", rep('',99),'Anchor',rep('',298),'End',rep('',99),'1kb')


df = as_tibble(res_combined, .name_repair = 'check_unique')
df <- df %>%
  rowid_to_column("row") %>%
  gather(col, Signals, -row) %>%
  mutate(col = as.integer(col))

my.lines<-data.frame(x=c(-1,301), y=c(0,0), xend=c(-1,301), yend=c(nrow(res_combined),nrow(res_combined)))


HM = ggplot(df, aes(col, row, fill = Signals)) +
  geom_tile() +ylab("")+
  scale_x_continuous(expand=c(0,0),breaks=c(-100,1,300,400),labels = c('-10kb','Anchor','End','+10kb'))+
  scale_fill_gradientn(colours = viridis(256)) +     # Use your custom colour palette
  #theme_void() +                                   # Minimal theme
  theme(
    axis.title.x=element_blank(), #axis.text.x=element_blank()
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(),
    axis.text.x = element_text(vjust=6, size = 12),
    panel.background = element_blank(),
    #panel.border=element_blank(),
    plot.margin=unit(c(0,0.9,0,2.8),'lines'),
    plot.title = element_text(hjust = 1),        # Right-aligned text
    panel.border=element_blank(),
    legend.position="bottom") +                  # Legend at the bottom
  guides(fill = guide_colourbar(
    title.position = "bottom",                   # Legend title below bar
    barwidth = 15,                               # Extend bar length
    title.hjust = 0.5))+
  geom_segment(data=my.lines, aes(x,y,xend=xend, yend=yend), size=1, inherit.aes=F,linetype='solid')

res_combined = cbind(res_b, res, res_a)
#D = apply(res_combined, 2, function(x){a=quantile(res_combined,0.5,na.rm=T); length(which(x>a))/length(x)})
D = apply(res_combined, 2, mean)
my.lines2<-data.frame(x=c(100,400), y=c(min(D),min(D)), xend=c(100,400), yend=c(max(D)+0.01,max(D)+0.01))
DF  = data.frame(Position=1:500, Prop = D)
L = ggplot(data = DF, aes(x=Position, y=Prop, group=1))+geom_line()+ylab("Average signal")+
  theme(
    panel.background = element_blank(), 
    #panel.border = element_rect(),
    axis.title.x=element_blank(), 
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.margin = unit(c(0,0,0,0.0),'lines')
    )+
  theme(axis.line = element_line(color='black'))
  #geom_segment(data=my.lines2, aes(x,y,xend=xend, yend=yend), size=1, inherit.aes=F,linetype='dashed',col='blue')

D = apply(res_combined, 2, function(x){a=quantile(res_combined,0.75,na.rm=T); length(which(x>=a))/length(x)})
#D = apply(res_combined, 2, mean)
my.lines2<-data.frame(x=c(100,400), y=c(min(D),min(D)), xend=c(100,400), yend=c(max(D)+0.01,max(D)+0.01))
DF  = data.frame(Position=1:500, Prop = D)
L2 = ggplot(data = DF, aes(x=Position, y=Prop, group=1))+geom_line()+ylab("Average signal")+
  theme(
    panel.background = element_blank(),
    #panel.border = element_rect(),
    axis.title.x=element_blank(), 
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.margin = unit(c(0,0,0,0.0),'lines')
  )+
  theme(axis.line = element_line(color='black'))

D = apply(res_combined, 2, function(x){a=quantile(res_combined,0.5,na.rm=T); length(which(x>=a))/length(x)})
#D = apply(res_combined, 2, mean)
my.lines2<-data.frame(x=c(100,400), y=c(min(D),min(D)), xend=c(100,400), yend=c(max(D)+0.01,max(D)+0.01))
DF  = data.frame(Position=1:500, Prop = D)
L3 = ggplot(data = DF, aes(x=Position, y=Prop, group=1))+geom_line()+ylab("Average signal")+
  theme(
    panel.background = element_blank(),
    #panel.border = element_rect(),
    axis.title.x=element_blank(), 
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.margin = unit(c(0,0,0,0.0),'lines')
  )+
  theme(axis.line = element_line(color='black'))

D = apply(res_combined, 2, median)
my.lines2<-data.frame(x=c(100,400), y=c(min(D),min(D)), xend=c(100,400), yend=c(max(D)+0.01,max(D)+0.01))
DF  = data.frame(Position=1:500, Prop = D)
L4 = ggplot(data = DF, aes(x=Position, y=Prop, group=1))+geom_line()+ylab("Average signal")+
  theme(
    panel.background = element_blank(), 
    #panel.border = element_rect(),
    axis.title.x=element_blank(), 
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.margin = unit(c(0,0,0,0.0),'lines')
  )+
  theme(axis.line = element_line(color='black'))


pdf(options[3], width = 5, height=8)
R = ggarrange(L, HM, ncol=1, nrow=2, heights=c(0.8,2.5))
print(R)
dev.off()


pdf(gsub('.pdf','2.pdf',options[3]), width = 5, height=8)
R = ggarrange(L2, HM, ncol=1, nrow=2, heights=c(0.8,2.5))
print(R)
dev.off()

pdf(gsub('.pdf','3.pdf',options[3]), width = 5, height=8)
R = ggarrange(L3, HM, ncol=1, nrow=2, heights=c(0.8,2.5))
print(R)
dev.off()

pdf(gsub('.pdf','4.pdf',options[3]), width = 5, height=8)
R = ggarrange(L4, HM, ncol=1, nrow=2, heights=c(0.8,2.5))
print(R)
dev.off()
