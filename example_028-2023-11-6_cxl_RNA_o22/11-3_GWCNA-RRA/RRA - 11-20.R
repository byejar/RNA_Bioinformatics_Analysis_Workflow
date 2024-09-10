# conda activate DESeq2
library(RobustRankAggreg)
library(dplyr)



setwd('/scratch/wus/028-2023-11-6_cxl_RNA_o22/8_DEseq2/')

flist<-c('condition_E64_vs_WT.txt','condition_E32_vs_WT.txt','condition_E258_vs_WT.txt','condition_E232_vs_WT.txt')

glist_up <- list()
glist_dw <- list()
for (i in flist){
  df<-read.table(paste('/scratch/wus/028-2023-11-6_cxl_RNA_o22/8_DEseq2/',i,sep=''))
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
#应用RRA算法，对基因进行整合排序
ag=aggregateRanks(glist_up)
#添加基因出现的次数
ag$freq_up=freq_up[match(ag$Name,freq_up$Var1),2]
ag

write.table(ag,'RRA_Emut_1WT_up.txt',sep='\t')


#统计所有基因出现的次数
freq_dw=as.data.frame(table(unlist(glist_dw)))
#应用RRA算法，对基因进行整合排序
ag=aggregateRanks(glist_dw)
#添加基因出现的次数
ag$freq_dw=freq_up[match(ag$Name,freq_dw$Var1),2]
ag

write.table(ag,'RRA_Emut_1WT_dw.txt',sep='\t')


###########
flist<-c('condition_A65_vs_WT.txt','condition_A299_vs_WT.txt','condition_A120_vs_WT.txt')

glist_up <- list()
glist_dw <- list()
for (i in flist){
  df<-read.table(paste('/scratch/wus/028-2023-11-6_cxl_RNA_o22/8_DEseq2/',i,sep=''))
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
#应用RRA算法，对基因进行整合排序
ag=aggregateRanks(glist_up)
#添加基因出现的次数
ag$freq_up=freq_up[match(ag$Name,freq_up$Var1),2]
ag

write.table(ag,'RRA_Amut_1WT_up.txt',sep='\t')


#统计所有基因出现的次数
freq_dw=as.data.frame(table(unlist(glist_dw)))
#应用RRA算法，对基因进行整合排序
ag=aggregateRanks(glist_dw)
#添加基因出现的次数
ag$freq_dw=freq_up[match(ag$Name,freq_dw$Var1),2]
ag

write.table(ag,'RRA_Amut_1WT_dw.txt',sep='\t')