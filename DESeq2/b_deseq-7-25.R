
# #envirment:DESeq2
library(DESeq2)
library(dplyr)


# #重排列
# #dataframe %>% select(order(colnames(dataframe)))

# ############ merge

# WT460<-read.table('/home/wus/2023_3_19-009-cxl_RNA/HCT116_WT460-20230522/CleanData/pipline_2023-3-21_RNA-analysis/pipline_2023-3-21_RNA-analysis/6_featureCounts/h.genename.matrix',header=1)
# others<-read.table('/home/wus/2023_3_19-009-cxl_RNA/TPL202303714+TPL202303545/pipline_2023-3-23_RNA-analysis/6_featureCounts/h.genename.matrix',header=1)
# merge_h=merge(WT460,others)

# write.csv(merge_h,'/home/wus/2023_3_19-009-cxl_RNA/merge_file-6-11/merge_h.csv',row.names = FALSE)


# rm(list=ls())
# library(DESeq2) 
# library(dplyr)

# df<-read.csv('/home/wus/2023_3_19-009-cxl_RNA/merge_file-6-11/3.csv',header=1,row.names=1)
# #重排列
# m_base <- df %>% select(order(colnames(df)))

# m_base[1:4,]

# condition <- factor(c(rep("KR234",4), rep("KR26",4), rep("KR274",4), rep("KR365",4), rep("KR692",4), rep("KR76",4), rep("WT",4), rep("WT450",4),rep("WT460",4))) 

   
# coldata <- data.frame(row.names = colnames(m_base), condition)
# coldata 

# dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)


# keep <- rowSums(counts(dds) >= 10) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数  
# dds <- dds[keep, ] 

# ####跑一下层次聚类和pca
# library('factoextra')          
# vsd <- vst(dds,blind = TRUE)    #参数blind=TRUE是为了不把样本分组信息考虑在内——即以无偏的方式进行转换
# sampleDists <- dist(t(assay(vsd))) 
# res1 <- hcut(sampleDists, k = 9, stand = FALSE,hc_method ="average" ) 
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


# #################### wt460 vs KR(del KR26)

# library(DESeq2) 
# library(dplyr)

# df<-read.csv('/home/wus/2023_3_19-009-cxl_RNA/merge_file-6-11/3.csv',header=1,row.names=1)

# df2<-df[,!grepl('WT_450|WT_[1-4]_|KR_26_3',colnames(df))]

# m_base <- df2 %>% select(order(colnames(df2)))
# m_base[1:4,]

# #condition <- factor(c(rep("KR234",4), rep("KR274",4), rep("KR365",4), rep("KR692",4), rep("KR76",4), rep("WT460",4))) 

# condition <- factor(c(rep("KR",23),  rep("WT460",4)))   
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
# write.csv(resdata,file = "/home/wus/2023_3_19-009-cxl_RNA/merge_file-6-11/WT460_WT/h_WT460_VS_KR(delete_KR26_3).csv")  
# write.table(resdata,file = "/home/wus/2023_3_19-009-cxl_RNA/merge_file-6-11/WT460_WT/h_WT460_VS_KR(delete_KR26_3).txt",sep='\t')

# ########################  

# # df<-read.csv('/home/wus/2023_3_19-009-cxl_RNA/merge_file-6-11/2.csv',header=1,row.names=1)

# # WT460<-df[,grepl('WT_460',colnames(df))]

# # WT<-df[,grepl('WT_[1-4]_',colnames(df))]

# # One sample t test : 单一样本t检验
# # Comparison of an observed mean with a
# # a theoretical mean
# #t.test(x, mu=0) 
# # Paired t test ：配对样本t检验
# #t.test(x, y, paired=TRUE)
# # Independent t test：独立样本t检验
# # Comparison of the means of two independent samples (x & y)
# #t.test(x, y)
# # Paired t test
# #t.test(x, y, paired=TRUE)


# #t.test(WT, WT460)

# ################################### wt wt460
# df<-read.csv('/home/wus/2023_3_19-009-cxl_RNA/merge_file-6-11/3.csv',header=1,row.names=1)

# df2<-df[,!grepl('WT_450|KR',colnames(df))]

# m_base <- df2 %>% select(order(colnames(df2)))
# m_base[1:4,]

# #condition <- factor(c(rep("KR234",4), rep("KR274",4), rep("KR365",4), rep("KR692",4), rep("KR76",4), rep("WT460",4))) 

# condition <- factor(c(rep("WT",4),  rep("WT460",4)))   
# coldata <- data.frame(row.names = colnames(m_base), condition)
# coldata 

# dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

# ############## 过滤方法
# keep <- rowSums(counts(dds) >= 10) >= 3  #过滤低表达基因，至少有3个样品都满足10个以上的reads数  
# dds <- dds[keep, ] 


# ###

# dds1 <- DESeq(dds)    # 做DESeq分析，输入的只能是原始矩阵，虽然有多种标准化方法，但是DESeq2在做差异表达分析的时候会用它自己的方法进行标准化。估计方差稳定变换只是为了做聚类或PCA
# resultsNames(dds1)    # 查看结果的名称。
# dds1$condition        #默认后者的处理组比前面的对照组。
# res <- results(dds1)  # 必要步骤！！！
# summary(res)          #看一下结果的概要信息，p值默认小于0.1。

# table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。
# res <- res[order(res$padj),]  #按照padj 进行升序排列
# resdata <- merge(as.data.frame(res), as.data.frame(counts(dds1, normalized=TRUE)),by="row.names",sort=FALSE)#保存时保存的是count标准化结果
# write.csv(resdata,file = "/home/wus/2023_3_19-009-cxl_RNA/merge_file-6-11/WT460_WT/h_WT460_VS_WT.csv")  
# write.table(resdata,file = "/home/wus/2023_3_19-009-cxl_RNA/merge_file-6-11/WT460_WT/h_WT460_VS_WT.txt",sep='\t') 


# #################### wt460-and-WT vs KR(del KR26)

# library(DESeq2) 
# library(dplyr)

# df<-read.csv('/home/wus/2023_3_19-009-cxl_RNA/merge_file-6-11/3.csv',header=1,row.names=1)

# df2<-df[,!grepl('WT_450|KR_26_3',colnames(df))]

# m_base <- df2 %>% select(order(colnames(df2)))
# m_base[1:4,]

# #condition <- factor(c(rep("KR234",4), rep("KR274",4), rep("KR365",4), rep("KR692",4), rep("KR76",4), rep("WT460",4))) 

# condition <- factor(c(rep("KR",23),  rep("WT",8)))   
# coldata <- data.frame(row.names = colnames(m_base), condition)
# coldata 

# dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)


# keep <- rowSums(counts(dds) >= 10) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数  
# dds <- dds[keep, ] 

# #####################

# dds1 <- DESeq(dds)    # 做DESeq分析，输入的只能是原始矩阵，虽然有多种标准化方法，但是DESeq2在做差异表达分析的时候会用它自己的方法进行标准化。估计方差稳定变换只是为了做聚类或PCA
# resultsNames(dds1)    # 查看结果的名称。
# dds1$condition        #默认后者的处理组比前面的对照组。
# res <- results(dds1)  # 必要步骤！！！
# summary(res)          #看一下结果的概要信息，p值默认小于0.1。

# table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。
# res <- res[order(res$padj),]  #按照padj 进行升序排列
# resdata <- merge(as.data.frame(res), as.data.frame(counts(dds1, normalized=TRUE)),by="row.names",sort=FALSE)#保存时保存的是count标准化结果
# write.csv(resdata,file = "/home/wus/2023_3_19-009-cxl_RNA/merge_file-6-11/WT460_WT/h_WT460-and-WT_VS_KR(delete_KR26_3).csv")  
# write.table(resdata,file = "/home/wus/2023_3_19-009-cxl_RNA/merge_file-6-11/WT460_WT/h_WT460-and-WT_VS_KR(delete_KR26_3).txt",sep='\t')  



# ########################## 两两
# library(DESeq2) 
# library(dplyr)

# setwd('/home/wus/2023_3_19-009-cxl_RNA/7-14/')

# df<-read.csv('/home/wus/2023_3_19-009-cxl_RNA/merge_file-6-11/3.csv',header=1,row.names=1)

# df2<-df[,!grepl('KR_26_3',colnames(df))]

# m_base <- df2 %>% select(order(colnames(df2)))
# m_base[1:4,]

# condition <- factor(c(rep("KR234",4), rep("KR26",3), rep("KR274",4), rep("KR365",4), rep("KR692",4), rep("KR76",4), rep("WT",4), rep("WT450",4), rep("WT460",4))) 
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
# for (i in resultsNames(dds)[c(2,3,4,5,6,7)]){
#      df<-read.table(paste(i,'.txt',sep=''),header=1,row.names=1)
#      df[1:4,1:4]
#      a<-filter(df, (log2FoldChange > 1 | log2FoldChange < -1) & (pvalue<0.05))[,c("Row.names","log2FoldChange","pvalue")]
#      colnames(a)<-c("Row.names",paste(i,'log2FC',sep='_'),paste(i,'pvalue',sep='_'))     
#           if(n==0){m=a}
#           else{m<-merge(m,a,by='Row.names')}
#      n=n+1

# }

# write.table(m,file = 'tumour_WT_VS_KR.txt',sep='\t',row.names=F) 

# #####
# setwd('/home/wus/2023_3_19-009-cxl_RNA/7-14/')

# df<-read.csv('/home/wus/2023_3_19-009-cxl_RNA/merge_file-6-11/3.csv',header=1,row.names=1)

# df2<-df[,!grepl('KR_26_3',colnames(df))]

# m_base <- df2 %>% select(order(colnames(df2)))
# m_base[1:4,]

# condition <- factor(c(rep("KR234",4), rep("KR26",3), rep("KR274",4), rep("KR365",4), rep("KR692",4), rep("KR76",4), rep("WT",4), rep("WT450",4), rep("WT460",4))) 
# coldata <- data.frame(row.names = colnames(m_base), condition)
# coldata$condition <- relevel(coldata$condition, ref = "WT460")
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
# for (i in resultsNames(dds)[c(2,3,4,5,6,7)]){
#      df<-read.table(paste(i,'.txt',sep=''),header=1,row.names=1)
#      df[1:4,1:4]
#      a<-filter(df, (log2FoldChange > 1 | log2FoldChange < -1) & (pvalue<0.05))[,c("Row.names","log2FoldChange","pvalue")]
#      colnames(a)<-c("Row.names",paste(i,'log2FC',sep='_'),paste(i,'pvalue',sep='_'))     
#           if(n==0){m=a}
#           else{m<-merge(m,a,by='Row.names')}
#      n=n+1

# }

# write.table(m,file = 'tumour_WT460_VS_KR.txt',sep='\t',row.names=F) 





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



########################
setwd('/home/wus/017-2023-6-28_cxl/dish/Projrct/pipline_2023-3-21_RNA-analysis/7_DESeq2/7-25/')

df<-read.table('../../6_featureCounts/h2',header=1,row.names=1)

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
for (i in resultsNames(dds)[c(2:7)]){
     df<-read.table(paste(i,'.txt',sep=''),header=1,row.names=1)
     df[1:4,1:4]
     print(i)
     a<-filter(df, (log2FoldChange > 1 | log2FoldChange < -1) & (pvalue<0.05))[,c("Row.names","log2FoldChange","pvalue")]
     colnames(a)<-c("Row.names",paste(i,'log2FC',sep='_'),paste(i,'pvalue',sep='_'))     
          if(n==0){m=a
                    print(n)}
          else{m<-merge(m,a,by='Row.names')
               print(n)}
     n=n+1

}


write.table(m,file = 'HCT116_cellline_KR_VS_WT460.txt',sep='\t',row.names=F) 


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

############
setwd('/home/wus/017-2023-6-28_cxl/dish/Projrct/pipline_2023-3-21_RNA-analysis/7_DESeq2/7-25/')
df<-read.table('../../6_featureCounts/h2',header=1,row.names=1)

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
for (i in resultsNames(dds)[c(2:7)]){
     df<-read.table(paste(i,'.txt',sep=''),header=1,row.names=1)
     df[1:4,1:4]
     print(i)
     a<-filter(df, (log2FoldChange > 1 | log2FoldChange < -1) & (pvalue<0.05))[,c("Row.names","log2FoldChange","pvalue")]
     colnames(a)<-c("Row.names",paste(i,'log2FC',sep='_'),paste(i,'pvalue',sep='_'))     
          if(n==0){m=a
                    print(n)}
          else{m<-merge(m,a,by='Row.names')
               print(n)}
     n=n+1

}


write.table(m,file = 'HCT116_cellline_KR_VS_WT.txt',sep='\t',row.names=F) 






############# 合并 wt parental
## 取变化方向相同的
WT_VS_KR<-read.table('HCT116_cellline_KR_VS_WT460.txt',header=1,row.names=1)
head(WT_VS_KR)
WT_VS_KR<-WT_VS_KR[,!grepl('pvalue',colnames(WT_VS_KR))]
WT_VS_KR$sum <- rowSums(sign(WT_VS_KR))
WT_VS_KR<-WT_VS_KR[WT_VS_KR$sum != 0,]
WT_VS_KR <- select(WT_VS_KR,-sum)
write.table(WT_VS_KR,file = 'HCT116_cellline_KR_VS_WT460-consistent.txt',sep='\t',row.names=T) 

#### 
parental_VS_KR<-read.table('HCT116_cellline_KR_VS_WT.txt',header=1,row.names=1)
head(parental_VS_KR)
parental_VS_KR<-parental_VS_KR[,!grepl('pvalue',colnames(parental_VS_KR))]
parental_VS_KR$sum <- rowSums(sign(parental_VS_KR))
parental_VS_KR<-parental_VS_KR[parental_VS_KR$sum != 0,]
parental_VS_KR <- select(parental_VS_KR,-sum)
head(parental_VS_KR)
write.table(parental_VS_KR,file = 'HCT116_cellline_KR_VS_WT-consistent.txt',sep='\t',row.names=T) 

####
m<-merge(WT_VS_KR,parental_VS_KR, by='row.names',all=TRUE)
write.table(m,file = 'HCT116_cellline_KR_VS_WT-or-WT460-consistent.txt',sep='\t',row.names=F) 