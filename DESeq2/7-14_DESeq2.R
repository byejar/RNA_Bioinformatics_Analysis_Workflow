############################# 
library(DESeq2) 
library(dplyr)

#/home/wus/017-2023-6-28_cxl/dish/Projrct/pipline_2023-3-21_RNA-analysis

setwd('/home/wus/017-2023-6-28_cxl/dish/Projrct/pipline_2023-3-21_RNA-analysis/')
df<-read.table('./6_featureCounts/h2',header=1,row.names=1)


# #df2<-df[,!grepl('WT_450|KR_26_3',colnames(df))]
# df2<-df
# m_base <- df2 %>% select(order(colnames(df2)))
# m_base[1:4,]

# condition <- factor(c(rep("KR234",3), rep("KR26",3), rep("KR274",3), rep("KR365",3), rep("KR692",3), rep("KR76",3), rep("WT",3), rep("WT450",3), rep("WT460",3))) 

# #condition <- factor(c(rep("KR",23),  rep("WT",4)))   
# coldata <- data.frame(row.names = colnames(m_base), condition)
# coldata 

# dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)


# keep <- rowSums(counts(dds) >= 10) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数  
# dds <- dds[keep, ] 

# ####跑一下层次聚类和pca
# library('factoextra')          
# vsd <- vst(dds,blind = TRUE)    #参数blind=TRUE是为了不把样本分组信息考虑在内——即以无偏的方式进行转换
# sampleDists <- dist(t(assay(vsd))) 
# res1 <- hcut(sampleDists, k = 8, stand = FALSE,hc_method ="average" ) 
# fviz_dend(res1,
#           rect_fill = T,
#           # 字体大小
#           cex = 1,
#           # 字体颜色
#           color_labels_by_k=T,
#           # 平行放置
#           horiz=T)

# rld <- vst(dds, blind=FALSE)     #vst()标准化。
# plotPCA(rld, intgroup="condition",ntop=500)     

 

# ############################# 删除WT450
# df<-read.table('../6_featureCounts/h2',header=1,row.names=1)

# df2<-df[,!grepl('h_WT450',colnames(df))]
# m_base <- df2 %>% select(order(colnames(df2)))
# m_base[1:4,]

# #condition <- factor(c(rep("KR234",3), rep("KR26",3), rep("KR274",3), rep("KR365",3), rep("KR692",3), rep("KR76",3), rep("WT",3), rep("WT450",3), rep("WT460",3))) 

# condition <- factor(c(rep("KR",18),  rep("WT&WT460",6)))   
# coldata <- data.frame(row.names = colnames(m_base), condition)
# coldata 

# dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)


# keep <- rowSums(counts(dds) >= 10) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数  
# dds <- dds[keep, ] 

# ####跑一下层次聚类和pca
# library('factoextra')          
# vsd <- vst(dds,blind = TRUE)    #参数blind=TRUE是为了不把样本分组信息考虑在内——即以无偏的方式进行转换
# sampleDists <- dist(t(assay(vsd))) 
# res1 <- hcut(sampleDists, k = 8, stand = FALSE,hc_method ="average" ) 
# fviz_dend(res1,
#           rect_fill = T,
#           # 字体大小
#           cex = 1,
#           # 字体颜色
#           color_labels_by_k=T,
#           # 平行放置
#           horiz=T)

# rld <- vst(dds, blind=FALSE)     #vst()标准化。
# plotPCA(rld, intgroup="condition",ntop=500)     



# dds1 <- DESeq(dds)    # 做DESeq分析，输入的只能是原始矩阵，虽然有多种标准化方法，但是DESeq2在做差异表达分析的时候会用它自己的方法进行标准化。估计方差稳定变换只是为了做聚类或PCA
# resultsNames(dds1)    # 查看结果的名称。
# dds1$condition        #默认后者的处理组比前面的对照组。
# res <- results(dds1)  # 必要步骤！！！
# summary(res)          #看一下结果的概要信息，p值默认小于0.1。

# table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。
# res <- res[order(res$padj),]  #按照padj 进行升序排列
# resdata <- merge(as.data.frame(res), as.data.frame(counts(dds1, normalized=TRUE)),by="row.names",sort=FALSE)#保存时保存的是count标准化结果
# write.csv(resdata,file = "h_WT&WT460_VS_KR.csv") 
# #############################  删除wt wt460组

# df2<-df[,!grepl('h_P|WT460',colnames(df))]

# m_base <- df2 %>% select(order(colnames(df2)))
# m_base[1:4,]


# condition <- factor(c(rep("KR",18),  rep("WT450",3)))   
# coldata <- data.frame(row.names = colnames(m_base), condition)
# coldata 

# dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)


# keep <- rowSums(counts(dds) >= 10) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数  
# dds <- dds[keep, ] 

# ####跑一下层次聚类和pca
# library('factoextra')          
# vsd <- vst(dds,blind = TRUE)    #参数blind=TRUE是为了不把样本分组信息考虑在内——即以无偏的方式进行转换
# sampleDists <- dist(t(assay(vsd))) 
# res1 <- hcut(sampleDists, k = 8, stand = FALSE,hc_method ="average" ) 
# fviz_dend(res1,
#           rect_fill = T,
#           # 字体大小
#           cex = 1,
#           # 字体颜色
#           color_labels_by_k=T,
#           # 平行放置
#           horiz=T)

# rld <- vst(dds, blind=FALSE)     #vst()标准化。
# plotPCA(rld, intgroup="condition",ntop=500)     



# dds1 <- DESeq(dds)    # 做DESeq分析，输入的只能是原始矩阵，虽然有多种标准化方法，但是DESeq2在做差异表达分析的时候会用它自己的方法进行标准化。估计方差稳定变换只是为了做聚类或PCA
# resultsNames(dds1)    # 查看结果的名称。
# dds1$condition        #默认后者的处理组比前面的对照组。
# res <- results(dds1)  # 必要步骤！！！
# summary(res)          #看一下结果的概要信息，p值默认小于0.1。

# table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。
# res <- res[order(res$padj),]  #按照padj 进行升序排列
# resdata <- merge(as.data.frame(res), as.data.frame(counts(dds1, normalized=TRUE)),by="row.names",sort=FALSE)#保存时保存的是count标准化结果
# write.csv(resdata,file = "h_WT_450_VS_KR.csv")  

# #############################  删除WT460 WT450 组

# df2<-df[,!grepl('WT450|WT460',colnames(df))]

# m_base <- df2 %>% select(order(colnames(df2)))
# m_base[1:4,]

# #condition <- factor(c(rep("KR234",4), rep("KR26",3), rep("KR274",4), rep("KR365",4), rep("KR692",4), rep("KR76",4), rep("WT",4), rep("WT_450",4))) 

# condition <- factor(c(rep("KR",18),  rep("WT",3)))   
# coldata <- data.frame(row.names = colnames(m_base), condition)
# coldata 

# dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)


# keep <- rowSums(counts(dds) >= 10) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数  
# dds <- dds[keep, ] 

# ####跑一下层次聚类和pca
# library('factoextra')          
# vsd <- vst(dds,blind = TRUE)    #参数blind=TRUE是为了不把样本分组信息考虑在内——即以无偏的方式进行转换
# sampleDists <- dist(t(assay(vsd))) 
# res1 <- hcut(sampleDists, k = 8, stand = FALSE,hc_method ="average" ) 
# fviz_dend(res1,
#           rect_fill = T,
#           # 字体大小
#           cex = 1,
#           # 字体颜色
#           color_labels_by_k=T,
#           # 平行放置
#           horiz=T)

# rld <- vst(dds, blind=FALSE)     #vst()标准化。
# plotPCA(rld, intgroup="condition",ntop=500)     



# dds1 <- DESeq(dds)    # 做DESeq分析，输入的只能是原始矩阵，虽然有多种标准化方法，但是DESeq2在做差异表达分析的时候会用它自己的方法进行标准化。估计方差稳定变换只是为了做聚类或PCA
# resultsNames(dds1)    # 查看结果的名称。
# dds1$condition        #默认后者的处理组比前面的对照组。
# res <- results(dds1)  # 必要步骤！！！
# summary(res)          #看一下结果的概要信息，p值默认小于0.1。

# table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。
# res <- res[order(res$padj),]  #按照padj 进行升序排列
# resdata <- merge(as.data.frame(res), as.data.frame(counts(dds1, normalized=TRUE)),by="row.names",sort=FALSE)#保存时保存的是count标准化结果
# write.csv(resdata,file = "h_WT_VS_KR.csv")  


# #############################  删除wt wt450组

# df2<-df[,!grepl('h_P|WT450',colnames(df))]

# m_base <- df2 %>% select(order(colnames(df2)))
# m_base[1:4,]


# condition <- factor(c(rep("KR",18),  rep("WT460",3)))   
# coldata <- data.frame(row.names = colnames(m_base), condition)
# coldata 

# dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)


# keep <- rowSums(counts(dds) >= 10) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数  
# dds <- dds[keep, ] 

# ####跑一下层次聚类和pca
# library('factoextra')          
# vsd <- vst(dds,blind = TRUE)    #参数blind=TRUE是为了不把样本分组信息考虑在内——即以无偏的方式进行转换
# sampleDists <- dist(t(assay(vsd))) 
# res1 <- hcut(sampleDists, k = 8, stand = FALSE,hc_method ="average" ) 
# fviz_dend(res1,
#           rect_fill = T,
#           # 字体大小
#           cex = 1,
#           # 字体颜色
#           color_labels_by_k=T,
#           # 平行放置
#           horiz=T)

# rld <- vst(dds, blind=FALSE)     #vst()标准化。
# plotPCA(rld, intgroup="condition",ntop=500)     



# dds1 <- DESeq(dds)    # 做DESeq分析，输入的只能是原始矩阵，虽然有多种标准化方法，但是DESeq2在做差异表达分析的时候会用它自己的方法进行标准化。估计方差稳定变换只是为了做聚类或PCA
# resultsNames(dds1)    # 查看结果的名称。
# dds1$condition        #默认后者的处理组比前面的对照组。
# res <- results(dds1)  # 必要步骤！！！
# summary(res)          #看一下结果的概要信息，p值默认小于0.1。

# table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。
# res <- res[order(res$padj),]  #按照padj 进行升序排列
# resdata <- merge(as.data.frame(res), as.data.frame(counts(dds1, normalized=TRUE)),by="row.names",sort=FALSE)#保存时保存的是count标准化结果
# write.csv(resdata,file = "h_WT_460_VS_KR.csv")  


#################### 两两比较
df2<-df

m_base <- df2 %>% select(order(colnames(df2)))
m_base[1:4,]

condition <- factor(c(rep("KR234",3), rep("KR26",3), rep("KR274",3), rep("KR365",3), rep("KR692",3), rep("KR76",3), rep("WT",3), rep("WT450",3), rep("WT460",3))) 
coldata <- data.frame(row.names = colnames(m_base), condition)
coldata$condition <- relevel(coldata$condition, ref = "WT")
coldata 

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)


keep <- rowSums(counts(dds) >= 10) >= 10  #过滤低表达基因，至少有n个样品都满足10个以上的reads数  
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
     write.table(resdata,file = paste(i,'.txt',sep='')) 
}

n=0
for (i in resultsNames(dds)[c(2,3,4,5,6,7)]){
     df<-read.table(paste(i,'.txt',sep=''),header=1,row.names=1)
     df[1:4,1:4]
     a<-filter(df, (log2FoldChange > 1 | log2FoldChange < -1) & (pvalue<0.05))[,c("Row.names","log2FoldChange","pvalue")]
     colnames(a)<-c("Row.names",paste(i,'log2FC',sep='_'),paste(i,'pvalue',sep='_'))     
          if(n==0){m=a}
          else{m<-merge(m,a,by='Row.names')}
     n=n+1

}

write.table(m,file = 'dish_WT_VS_KR.txt',sep='\t',row.names=F) 


n=0
for (i in resultsNames(dds)[c(2,3,4,5,6,7)]){
     df<-read.table(paste(i,'.txt',sep=''),header=1,row.names=1)
     df[1:4,1:4]
     a<-filter(df, (log2FoldChange > -1 | log2FoldChange < 1) | (pvalue > 0.05))[,c("Row.names","log2FoldChange","pvalue")]
     colnames(a)<-c("Row.names",paste(i,'log2FC',sep='_'),paste(i,'pvalue',sep='_'))     
          if(n==0){m=a}
          else{m<-merge(m,a,by='Row.names')}
     n=n+1

}

write.table(m,file = 'dish_WT_VS_KR_non-change.txt',sep='\t',row.names=F) 

#################### 7-14 两两


df2<-df

m_base <- df2 %>% select(order(colnames(df2)))
m_base[1:4,]

condition <- factor(c(rep("KR234",3), rep("KR26",3), rep("KR274",3), rep("KR365",3), rep("KR692",3), rep("KR76",3), rep("WT",3), rep("WT450",3), rep("WT460",3))) 
coldata <- data.frame(row.names = colnames(m_base), condition)
coldata$condition <- relevel(coldata$condition, ref = "WT460")
coldata 

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)


keep <- rowSums(counts(dds) >= 10) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数  
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
     write.table(resdata,file = paste(i,'.txt',sep='')) 
}

n=0
for (i in resultsNames(dds)[c(2,3,4,5,6,7)]){
     df<-read.table(paste(i,'.txt',sep=''),header=1,row.names=1)
     df[1:4,1:4]
     a<-filter(df, (log2FoldChange > 1 | log2FoldChange < -1) & (pvalue<0.05))[,c("Row.names","log2FoldChange","pvalue")]
     colnames(a)<-c("Row.names",paste(i,'log2FC',sep='_'),paste(i,'pvalue',sep='_'))     
          if(n==0){m=a}
          else{m<-merge(m,a,by='Row.names')}
     n=n+1
}
write.table(m,file = 'dish_WT460_VS_KR.txt',sep='\t',row.names=F) 


n=0
for (i in resultsNames(dds)[c(2,3,4,5,6,7)]){
     df<-read.table(paste(i,'.txt',sep=''),header=1,row.names=1)
     df[1:4,1:4]
     a<-filter(df, (log2FoldChange > -1 | log2FoldChange < 1) | (pvalue > 0.05))[,c("Row.names","log2FoldChange","pvalue")]
     colnames(a)<-c("Row.names",paste(i,'log2FC',sep='_'),paste(i,'pvalue',sep='_'))     
          if(n==0){m=a}
          else{m<-merge(m,a,by='Row.names')}
     n=n+1
}
write.table(m,file = 'dish_WT460_VS_KR_non-change.txt',sep='\t',row.names=F) 




####
df2<-df

m_base <- df2 %>% select(order(colnames(df2)))
m_base[1:4,]

condition <- factor(c(rep("KR234",3), rep("KR26",3), rep("KR274",3), rep("KR365",3), rep("KR692",3), rep("KR76",3), rep("WT",3), rep("WT450",3), rep("WT460",3))) 
coldata <- data.frame(row.names = colnames(m_base), condition)
coldata$condition <- relevel(coldata$condition, ref = "WT450")
coldata 

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)


keep <- rowSums(counts(dds) >= 10) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数  
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
     write.table(resdata,file = paste(i,'.txt',sep='')) 
}


# "Row.names","log2FoldChange","pvalue"

# df<-read.csv('condition_KR76_vs_WT460.txt',header=1,row.names=1)
# a1<-filter(df, (log2FoldChange > 1 | log2FoldChange < -1) & (pvalue<0.05))[,c("Row.names","log2FoldChange","pvalue")]


# df<-read.csv('condition_KR692_vs_WT460.txt',header=1,row.names=1)
# a2<-filter(df, (log2FoldChange > 1 | log2FoldChange < -1) & (pvalue<0.05))[,c("Row.names","log2FoldChange","pvalue")]
# m1<-merge(a1,a2,by='Row.names')

# df<-read.csv('condition_KR365_vs_WT460.txt',header=1,row.names=1)
# a3<-filter(df, (log2FoldChange > 1 | log2FoldChange < -1) & (pvalue<0.05))[,c("Row.names","log2FoldChange","pvalue")]
# m2<-merge(m1,a3,by='Row.names')

# df<-read.csv('condition_KR274_vs_WT460.txt',header=1,row.names=1)
# a4<-filter(df, (log2FoldChange > 1 | log2FoldChange < -1) & (pvalue<0.05))[,c("Row.names","log2FoldChange","pvalue")]
# m3<-merge(m2,a4,by='Row.names')

# df<-read.csv('condition_KR26_vs_WT460.txt',header=1,row.names=1)
# a5<-filter(df, (log2FoldChange > 1 | log2FoldChange < -1) & (pvalue<0.05))[,c("Row.names","log2FoldChange","pvalue")]
# m4<-merge(m3,a5,by='Row.names')

# df<-read.csv('condition_KR234_vs_WT460.txt',header=1,row.names=1)
# a6<-filter(df, (log2FoldChange > 1 | log2FoldChange < -1) & (pvalue<0.05))[,c("Row.names","log2FoldChange","pvalue")]
# m5<-merge(m4,a6,by='Row.names')

