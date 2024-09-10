#conda activate DESeq2
library(DESeq2)
library(dplyr)


setwd('/home/wus/023_2023-9-7_cxl_RNA/023-clx-X101SC23062450-Z01-J004/cell_line_analysis/8_DEseq2/')
df<-read.table('../7_featureCounts/h2v2',header=1,row.names=1)
df2<-df

m_base <- df2 %>% select(order(colnames(df2)))
m_base[1:4,]
#重命名列
# colnames(m_base)<-c('N_C-1','N_P-1','N_C-2','N_P-2','N_C-3','N_P-3','N_C-4','N_P-4','S1_C-1','S1_P-1','S1_C-2','S1_P-2','S1_C-3','S1_P-3','S1_C-4','S1_P-4','S2_C-1','S2_P-1','S2_C-2','S2_P-2','S2_C-3','S2_P-3','S2_C-4','S2_P-4')
# m_base <- m_base %>% select(order(colnames(m_base)))
# m_base[1:4,]

condition <- factor(c(rep("306",6), rep("EV",4), rep("WT",4)))
coldata <- data.frame(row.names = colnames(m_base), condition)
#coldata$condition <- relevel(coldata$condition)
coldata

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

keep <- rowSums(counts(dds) >= 10) >= 5  #过滤低表达基因，至少有2个样品都满足10个以上的reads数

# # 预处理，过滤低丰度的数据
# countData <- count[apply(count, 1, sum) > 0 , ]

dds <- dds[keep, ]

# dds <- DESeq(dds, betaPrior = FALSE)
# resultsNames(dds)


library('factoextra')
vsd <- vst(dds,blind = TRUE)    #参数blind=TRUE是为了不把样本分组信息考虑在内——即以无偏的方式进行转换
sampleDists <- dist(t(assay(vsd)))
res1 <- hcut(sampleDists, k = 6, stand = FALSE,hc_method ="average" )
fviz_dend(res1,
          rect_fill = T,
          # 字体大小
          cex = 1,
          # 字体颜色
          color_labels_by_k=T,
          # 平行放置
          horiz=T)

rld <- vst(dds, blind=FALSE)     #vst()标准化。
plotPCA(rld, intgroup="condition",ntop=500)
# p<-plotPCA(rld, intgroup="condition",ntop=500)
# ggsave('/home/wus/019_2023-07-20_zc_BMDM/202307BMDM-si/cell_line/10_analysis/pca-8-10.png',p,dpi=600)


# p+theme_bw()+
# theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#   axis.ticks.y=element_line(color="black",size=0.8,lineend = 2),
#   #axis.line.y=element_line(linetype=1,color="black",size=0.4),
#   axis.ticks.x=element_line(color="black",size=0.8,lineend = 2),
#   panel.border = element_rect(color = "black", size = 0.6, fill = NA))


# pca_data<-plotPCA(rld, intgroup="condition",ntop=500,returnData=T)

# pca<-ggplot(pca_data,aes(x=PC1,y=PC2,color=condition))+ geom_point(size=2)+theme_bw()+
# theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#   axis.ticks.y=element_line(color="black",size=0.8,lineend = 2),
#   #axis.line.y=element_line(linetype=1,color="black",size=0.4),
#   axis.ticks.x=element_line(color="black",size=0.8,lineend = 2),
#   panel.border = element_rect(color = "black", size = 0.6, fill = NA))

# ggsave('/home/wus/019_2023-07-20_zc_BMDM/202307BMDM-si/cell_line/10_analysis/pca.png',pca,dpi=600)


########################################################################
coldata <- data.frame(row.names = colnames(m_base), condition)
coldata$condition <- relevel(coldata$condition, ref = "EV")
coldata

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

keep <- rowSums(counts(dds) >= 10) >= 5  #过滤低表达基因，至少有2个样品都满足10个以上的reads数

# # 预处理，过滤低丰度的数据
# countData <- count[apply(count, 1, sum) > 0 , ]

dds <- dds[keep, ]

dds <- DESeq(dds, betaPrior = FALSE)
resultsNames(dds)

length(resultsNames(dds))

for (i in resultsNames(dds)){
     print(i)
     res <- results(dds, name=i)
     summary(res)          #看一下结果的概要信息，p值默认小于0.1。

     table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。
     res <- res[order(res$padj),]  #按照padj 进行升序排列
     resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)#保存时保存的是count标准化结果
     #write.table(resdata,file = paste(i,'.txt',sep=''),sep='\t') 
}


res <- results(dds, name='condition_WT_vs_EV')
     summary(res,alpha=0.7)

#######################################################################
coldata <- data.frame(row.names = colnames(m_base), condition)
coldata$condition <- relevel(coldata$condition, ref = "WT")
coldata

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

keep <- rowSums(counts(dds) >= 10) >= 5  #过滤低表达基因，至少有2个样品都满足10个以上的reads数

# # 预处理，过滤低丰度的数据
# countData <- count[apply(count, 1, sum) > 0 , ]

dds <- dds[keep, ]

dds <- DESeq(dds, betaPrior = FALSE)
resultsNames(dds)

length(resultsNames(dds))

for (i in resultsNames(dds)){
     print(i)
     res <- results(dds, name=i)
     summary(res)          #看一下结果的概要信息，p值默认小于0.1。

     table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。
     res <- res[order(res$padj),]  #按照padj 进行升序排列
     resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)#保存时保存的是count标准化结果
     write.table(resdata,file = paste(i,'.txt',sep=''),sep='\t') 
}


### 改变输出的顺序
# results(dds, contrast=c("condition","KR293", "parental"))
# results(dds, contrast=c("condition", "parental","KR293"))
### 不方便，全改成处理组 vs con ？
### 所以这里的名字都应该改成KR vs con 




# coldata <- data.frame(row.names = colnames(m_base), condition)
# coldata$condition <- relevel(coldata$condition, ref = "N-C")
# coldata

# dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

# keep <- rowSums(counts(dds) >= 10) >= 5  #过滤低表达基因，至少有2个样品都满足10个以上的reads数

# # # 预处理，过滤低丰度的数据
# # countData <- count[apply(count, 1, sum) > 0 , ]

# dds <- dds[keep, ]

# dds <- DESeq(dds, betaPrior = FALSE)
# resultsNames(dds)

# length(resultsNames(dds))

# for (i in resultsNames(dds)){
#      print(i)
#      res <- results(dds, name=i)
#      summary(res)          #看一下结果的概要信息，p值默认小于0.1。

#      table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。
#      res <- res[order(res$padj),]  #按照padj 进行升序排列
#      resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)#保存时保存的是count标准化结果
#      write.table(resdata,file = paste(i,'.txt',sep=''),sep='\t') 
# }

# #######################################################################
# coldata <- data.frame(row.names = colnames(m_base), condition)
# coldata$condition <- relevel(coldata$condition, ref = "S1-C")
# coldata

# dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

# keep <- rowSums(counts(dds) >= 10) >= 5  #过滤低表达基因，至少有2个样品都满足10个以上的reads数

# # # 预处理，过滤低丰度的数据
# # countData <- count[apply(count, 1, sum) > 0 , ]

# dds <- dds[keep, ]

# dds <- DESeq(dds, betaPrior = FALSE)
# resultsNames(dds)

# length(resultsNames(dds))

# for (i in resultsNames(dds)){
#      print(i)
#      res <- results(dds, name=i)
#      summary(res)          #看一下结果的概要信息，p值默认小于0.1。

#      table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。
#      res <- res[order(res$padj),]  #按照padj 进行升序排列
#      resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)#保存时保存的是count标准化结果
#      write.table(resdata,file = paste(i,'.txt',sep=''),sep='\t') 
# }

# #######################################################################
# coldata <- data.frame(row.names = colnames(m_base), condition)
# coldata$condition <- relevel(coldata$condition, ref = "S2-C")
# coldata

# dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

# keep <- rowSums(counts(dds) >= 10) >= 5  #过滤低表达基因，至少有2个样品都满足10个以上的reads数

# # # 预处理，过滤低丰度的数据
# # countData <- count[apply(count, 1, sum) > 0 , ]

# dds <- dds[keep, ]

# dds <- DESeq(dds, betaPrior = FALSE)
# resultsNames(dds)

# length(resultsNames(dds))

# for (i in resultsNames(dds)){
#      print(i)
#      res <- results(dds, name=i)
#      summary(res)          #看一下结果的概要信息，p值默认小于0.1。

#      table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。
#      res <- res[order(res$padj),]  #按照padj 进行升序排列
#      resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)#保存时保存的是count标准化结果
#      write.table(resdata,file = paste(i,'.txt',sep=''),sep='\t') 
# }



# n=0
# for (i in resultsNames(dds)[c(2,3)]){
#      df<-read.table(paste(i,'.txt',sep=''),header=1,row.names=1)
#      df[1:4,1:4]
#      a<-filter(df, (log2FoldChange > 1 | log2FoldChange < -1) & (pvalue<0.05))[,c("Row.names","log2FoldChange","pvalue")]
#      colnames(a)<-c("Row.names",paste(i,'log2FC',sep='_'),paste(i,'pvalue',sep='_'))     
#           if(n==0){m=a}
#           else{m<-merge(m,a,by='Row.names')}
#      n=n+1

# }

# write.table(m,file = 'MC38_tumor_parental_VS_KR.txt',sep='\t',row.names=F) 


# n=0
# for (i in resultsNames(dds)[c(2,3)]){
#      df<-read.table(paste(i,'.txt',sep=''),header=1,row.names=1)
#      df[1:4,1:4]
#      a<-filter(df, (log2FoldChange > -1 | log2FoldChange < 1) | (pvalue > 0.05))[,c("Row.names","log2FoldChange","pvalue")]
#      colnames(a)<-c("Row.names",paste(i,'log2FC',sep='_'),paste(i,'pvalue',sep='_'))     
#           if(n==0){m=a}
#           else{m<-merge(m,a,by='Row.names')}
#      n=n+1

# }

# write.table(m,file = 'MC38_tumor_parental_VS_KR_non-change.txt',sep='\t',row.names=F) 



# ########################
# df<-read.table('../7_featureCounts/m2',header=1,row.names=1)
# df2<-df

# m_base <- df2 %>% select(order(colnames(df2)))
# m_base[1:4,]

# condition <- factor(c(rep("KR293",4), rep("KR83",4), rep("parental",4), rep("WT",4)))
# coldata <- data.frame(row.names = colnames(m_base), condition)
# coldata$condition <- relevel(coldata$condition, ref = "WT")
# coldata

# dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

# keep <- rowSums(counts(dds) >= 10) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数
# dds <- dds[keep, ]


# dds <- DESeq(dds, betaPrior = FALSE)
# resultsNames(dds)

# length(resultsNames(dds))

# for (i in resultsNames(dds)){
#      print(i)
#      res <- results(dds, name=i)
#      summary(res)          #看一下结果的概要信息，p值默认小于0.1。

#      table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。
#      res <- res[order(res$padj),]  #按照padj 进行升序排列
#      resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)#保存时保存的是count标准化结果
#      write.table(resdata,file = paste(i,'.txt',sep='')) 
# }

# n=0
# for (i in resultsNames(dds)[c(2,3)]){
#      df<-read.table(paste(i,'.txt',sep=''),header=1,row.names=1)
#      df[1:4,1:4]
#      a<-filter(df, (log2FoldChange > 1 | log2FoldChange < -1) & (pvalue<0.05))[,c("Row.names","log2FoldChange","pvalue")]
#      colnames(a)<-c("Row.names",paste(i,'log2FC',sep='_'),paste(i,'pvalue',sep='_'))     
#           if(n==0){m=a}
#           else{m<-merge(m,a,by='Row.names')}
#      n=n+1

# }


# write.table(m,file = 'MC38_tumor_WT_VS_KR.txt',sep='\t',row.names=F) 


# n=0
# for (i in resultsNames(dds)[c(2,3)]){
#      df<-read.table(paste(i,'.txt',sep=''),header=1,row.names=1)
#      df[1:4,1:4]
#      a<-filter(df, (log2FoldChange > -1 | log2FoldChange < 1) | (pvalue > 0.05))[,c("Row.names","log2FoldChange","pvalue")]
#      colnames(a)<-c("Row.names",paste(i,'log2FC',sep='_'),paste(i,'pvalue',sep='_'))     
#           if(n==0){m=a}
#           else{m<-merge(m,a,by='Row.names')}
#      n=n+1

# }

# write.table(m,file = 'MC38_tumor_WT_VS_KR_non-change.txt',sep='\t',row.names=F) 


# ############# 合并 wt parental
# ## 取变化方向相同的
# WT_VS_KR<-read.table('MC38_tumor_WT_VS_KR.txt',header=1,row.names=1)
# head(WT_VS_KR)
# WT_VS_KR<-WT_VS_KR[,!grepl('pvalue',colnames(WT_VS_KR))]
# WT_VS_KR$sum <- rowSums(sign(WT_VS_KR))
# WT_VS_KR<-WT_VS_KR[WT_VS_KR$sum != 0,]
# WT_VS_KR<-WT_VS_KR[,1:2]
# write.table(WT_VS_KR,file = 'MC38_tumor_WT_VS_KR-consistent.txt',sep='\t',row.names=T) 

# #### 
# parental_VS_KR<-read.table('MC38_tumor_parental_VS_KR.txt',header=1,row.names=1)
# head(parental_VS_KR)
# parental_VS_KR<-parental_VS_KR[,!grepl('pvalue',colnames(parental_VS_KR))]
# parental_VS_KR$sum <- rowSums(sign(parental_VS_KR))
# parental_VS_KR<-parental_VS_KR[parental_VS_KR$sum != 0,]
# parental_VS_KR<-parental_VS_KR[,1:2]
# head(parental_VS_KR)
# write.table(parental_VS_KR,file = 'MC38_tumor_parental_VS_KR-consistent.txt',sep='\t',row.names=T) 

# ####
# m<-merge(WT_VS_KR,parental_VS_KR, by='row.names',all=TRUE)
# write.table(m,file = 'MC38_tumor_WT&parental_VS_KR-consistent.txt',sep='\t',row.names=F) 

# ##### 添加注释
# bio<-read.table('../9_merge_file/src/Homo_sapiens.GRCh38.84_gene_name-gene_biotype.txt',header=1)
# df2<-merge(bio,m,by.y='Row.names',by.x='gene_name')
# write.table(df2,'MC38_tumor_WT&parental_VS_KR-consistent-biotype.txt',sep='\t')



setwd('/home/wus/022-2023-9-4-gmk-RNA/20230703-mTOR_K2066_mice_Thr_6h/cell_line_analysis/8_DEseq2/')
df<-read.table('../7_featureCounts/9-7/m2',header=1,row.names=1)
df2<-df
df2[1:4,]
m_base <- df2 %>% select(order(colnames(df2)))
m_base[1:4,]

#重命名列
colnames(m_base)<-c('Con-WT-1','Thr-WT-2','Thr-WT-3','Thr-WT-4','Thr-KR-1','Thr-KR-2','Thr-KR-3','Thr-KR-4','Con-WT-2','Con-WT-3','Con-WT-4','Con-KR-1','Con-KR-2','Con-KR-3','Con-KR-4','Thr-WT-1')
m_base <- m_base %>% select(order(colnames(m_base)))
m_base[1:4,]

condition <- factor(c(rep("Con-KR",4), rep("Con-WT",4), rep("Thr-KR",4), rep("Thr-WT",4)))
coldata <- data.frame(row.names = colnames(m_base), condition)
#coldata$condition <- relevel(coldata$condition)
coldata

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

keep <- rowSums(counts(dds) >= 10) >= 5  #过滤低表达基因，至少有2个样品都满足10个以上的reads数

# # 预处理，过滤低丰度的数据
# countData <- count[apply(count, 1, sum) > 0 , ]

dds <- dds[keep, ]

# dds <- DESeq(dds, betaPrior = FALSE)
# resultsNames(dds)


library('factoextra')
vsd <- vst(dds,blind = TRUE)    #参数blind=TRUE是为了不把样本分组信息考虑在内——即以无偏的方式进行转换
sampleDists <- dist(t(assay(vsd)))
res1 <- hcut(sampleDists, k = 4, stand = FALSE,hc_method ="average" )
fviz_dend(res1,
          rect_fill = T,
          # 字体大小
          cex = 1,
          # 字体颜色
          color_labels_by_k=T,
          # 平行放置
          horiz=T)

rld <- vst(dds, blind=FALSE)     #vst()标准化。
plotPCA(rld, intgroup="condition",ntop=500)


########################################################################
coldata <- data.frame(row.names = colnames(m_base), condition)
coldata$condition <- relevel(coldata$condition, ref = "Con.KR")
coldata

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

keep <- rowSums(counts(dds) >= 10) >= 5  #过滤低表达基因，至少有2个样品都满足10个以上的reads数

# # 预处理，过滤低丰度的数据
# countData <- count[apply(count, 1, sum) > 0 , ]

dds <- dds[keep, ]

dds <- DESeq(dds, betaPrior = FALSE)
resultsNames(dds)

length(resultsNames(dds))

for (i in resultsNames(dds)){
     print(i)
     res <- results(dds, name=i)
     summary(res)          #看一下结果的概要信息，p值默认小于0.1。

     table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。
     res <- res[order(res$padj),]  #按照padj 进行升序排列
     resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)#保存时保存的是count标准化结果
     write.table(resdata,file = paste(i,'.txt',sep=''),sep='\t') 
}

#######################################################################
# coldata <- data.frame(row.names = colnames(m_base), condition)
# coldata$condition <- relevel(coldata$condition, ref = "Con.WT")
# coldata

# dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

# keep <- rowSums(counts(dds) >= 10) >= 5  #过滤低表达基因，至少有2个样品都满足10个以上的reads数

# # # 预处理，过滤低丰度的数据
# # countData <- count[apply(count, 1, sum) > 0 , ]

# dds <- dds[keep, ]

# dds <- DESeq(dds, betaPrior = FALSE)
# resultsNames(dds)

# length(resultsNames(dds))

# for (i in resultsNames(dds)){
#      print(i)
#      res <- results(dds, name=i)
#      summary(res)          #看一下结果的概要信息，p值默认小于0.1。

#      table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。
#      res <- res[order(res$padj),]  #按照padj 进行升序排列
#      resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)#保存时保存的是count标准化结果
#      write.table(resdata,file = paste(i,'.txt',sep=''),sep='\t') 
# }

#######################################################################
coldata <- data.frame(row.names = colnames(m_base), condition)
coldata$condition <- relevel(coldata$condition, ref = "Thr.KR")
coldata

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

keep <- rowSums(counts(dds) >= 10) >= 5  #过滤低表达基因，至少有2个样品都满足10个以上的reads数

# # 预处理，过滤低丰度的数据
# countData <- count[apply(count, 1, sum) > 0 , ]

dds <- dds[keep, ]

dds <- DESeq(dds, betaPrior = FALSE)
resultsNames(dds)

length(resultsNames(dds))

for (i in resultsNames(dds)){
     print(i)
     res <- results(dds, name=i)
     summary(res)          #看一下结果的概要信息，p值默认小于0.1。

     table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。
     res <- res[order(res$padj),]  #按照padj 进行升序排列
     resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)#保存时保存的是count标准化结果
     write.table(resdata,file = paste(i,'.txt',sep=''),sep='\t') 
}

