library(RobustRankAggreg)
# library(clusterProfiler)
# set.seed(123456789)
# data(geneList, package="DOSE")
# head(geneList)
# deg=data.frame(gene=names(geneList),
#               logFC=as.numeric(geneList),
#               P.Value=0)
# head(deg)
# # 大家的deg 来源于真实的差异分析
# deg2=deg;deg2$gene=sample(deg$gene,length(deg2$gene))
# deg3=deg;deg3$gene=sample(deg$gene,length(deg2$gene))
# deg4=deg;deg4$gene=sample(deg$gene,length(deg2$gene))
# get_up <- function(df){
#   df$g=ifelse(df$P.Value>0.01,'stable', #if 判断：如果这一基因的P.Value>0.01，则为stable基因
#               ifelse( df$logFC >2,'up', #接上句else 否则：接下来开始判断那些P.Value<0.01的基因，再if 判断：如果logFC >1.5,则为up（上调）基因
#                       ifelse( df$logFC < -2,'down','stable') )#接上句else 否则：接下来开始判断那些logFC <1.5 的基因，再if 判断：如果logFC <1.5，则为down（下调）基因，否则为stable基因
#   )
#  print( table(df$g))
#   df=df[order(df$logFC,decreasing = T),]
#  # rownames(df[df$g=='up',])
#   df[df$g=='up','gene']
# }
# glist=list(get_up(deg2)
#            ,get_up(deg3)
#            , get_up(deg4))
# ups=aggregateRanks(glist = glist, N = length(unique(unlist(glist))))
# tmp=as.data.frame(table(unlist(glist)))
# ups$Freq=tmp[match(ups$Name,tmp[,1]),2]
# head(ups)

library(dplyr)

# df<-read.table('/home/wus/2023_3_19-009-cxl_RNA/7-25/condition_KR76_vs_WT.txt')
# T.genes<-df %>% dplyr::filter(pvalue < 0.05)%>%arrange(desc(log2FoldChange))%>%dplyr::select(Row.names,log2FoldChange)
# head(T.genes)




flist<-c('condition_KR234_vs_WT.txt','condition_KR26_vs_WT.txt','condition_KR274_vs_WT.txt','condition_KR365_vs_WT.txt','condition_KR692_vs_WT.txt','condition_KR76_vs_WT.txt')

glist_up <- list()
glist_dw <- list()
for (i in flist){
  df<-read.table(paste('/home/wus/2023_3_19-009-cxl_RNA/7-25/',i,sep=''))
  genes_up<-df %>% dplyr::filter(pvalue < 0.05,log2FoldChange > 1)%>%arrange(desc(log2FoldChange))%>%dplyr::select(Row.names)
  genes_up.list<-as.list(genes_up)
  glist_up <-append(glist_up,genes_up.list)
  print(glist_up)


  genes_dw<-df %>% dplyr::filter(pvalue < 0.05,log2FoldChange < -1)%>%arrange(desc(log2FoldChange))%>%dplyr::select(Row.names)
  genes_dw.list<-as.list(genes_dw)
  glist_dw <-append(glist_dw,genes_dw.list)
  print(glist_dw)
}


#统计所有基因出现的次数
freq_up=as.data.frame(table(unlist(glist_up)))
​
#应用RRA算法，对基因进行整合排序
ag=aggregateRanks(glist_up)
#添加基因出现的次数
ag$freq_up=freq_up[match(ag$Name,freq_up$Var1),2]
ag

write.table(ag,'/home/wus/2023_3_19-009-cxl_RNA/11-3_GWCNA-RRA/RRA_6KR_1WT_up.txt')


#统计所有基因出现的次数
freq_dw=as.data.frame(table(unlist(glist_dw)))
​
#应用RRA算法，对基因进行整合排序
ag=aggregateRanks(glist_dw)
#添加基因出现的次数
ag$freq_dw=freq_up[match(ag$Name,freq_dw$Var1),2]
ag

write.table(ag,'/home/wus/2023_3_19-009-cxl_RNA/11-3_GWCNA-RRA/RRA_6KR_1WT_dw.txt')