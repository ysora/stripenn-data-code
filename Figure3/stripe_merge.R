library(GenomicRanges)

setwd('/mnt/data0/sora/stripe')
#uni = read.delim('stripenn_out/NOD_vs_BL6/union_stripiness.txt',header=T)
#plot(uni$BL6_stripiness, uni$NOD_stripiness)
#uni[which(uni$NOD_stripiness>20),]

ch12 = read.delim('stripenn_out/CH12/CH12/result_filtered.tsv',header=T)
dn3 = read.delim('stripenn_out/DN3EV2/DN3EV_VCSQRT/result_filtered.tsv', header=T)

ch12 = ch12[which(ch12$pvalue<0.05),]
dn3 = dn3[which(dn3$pvalue<0.05),]

getOverlaps = function(A,B)
{
  A_x = A[,1:3]
  A_y = A[,4:6]
  B_x = B[,1:3]
  B_y = B[,4:6]
  colnames(A_x) = colnames(A_y) = colnames(B_x) = colnames(B_y) = c("chr",'start','end')
  
  ov_x = findOverlaps(makeGRangesFromDataFrame(A_x), makeGRangesFromDataFrame(B_x))
  ov_y = findOverlaps(makeGRangesFromDataFrame(A_y), makeGRangesFromDataFrame(B_y))
  
  ovx = paste(ov_x@from , ov_x@to, sep="_")
  ovy = paste(ov_y@from, ov_y@to, sep='_')
  ints = intersect(ovx, ovy)
  
  SPLIT = strsplit(ints, split="_")
  A_ind = as.numeric(unlist(lapply(SPLIT, function(x) x[1])))
  B_ind = as.numeric(unlist(lapply(SPLIT, function(x) x[2])))
  
  A_unique = setdiff(1:nrow(A), A_ind)
  B_unique = setdiff(1:nrow(B), B_ind)
  
  chooseA = c()
  chooseB = c()
  
  for(i in 1:length(A_ind))
  {
    stripinessA = A$pvalue[A_ind[i]]
    stripinessB = B$pvalue[B_ind[i]]
    if(stripinessA <= stripinessB)
    {
      chooseA = append(chooseA, A_ind[i])
    }else{
      chooseB = append(chooseB, B_ind[i])
    }
  }
  
  A_ind = c(A_unique, chooseA)
  B_ind = c(B_unique, chooseB)
  
  Asub = A[A_ind,]
  Bsub = B[B_ind,]
  result = rbind(Asub, Bsub)
  return(result)
}  

unions = getOverlaps(ch12, dn3)

dim(unions)
head(unions)

del = grep('.',rownames(unions),fixed=T)
if(length(del)>0){unions = unions[-del,]}
write.table(unions, 'CH12_DN3_compare/union_stripes.tsv' , sep = '\t', row.names = F, col.names = T, quote = F)
