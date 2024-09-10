
#envirment:DESeq2



#重排列
#dataframe %>% select(order(colnames(dataframe)))

############ merge

# WT460<-read.table('/home/wus/2023_3_19-009-cxl_RNA/HCT116_WT460-20230522/CleanData/pipline_2023-3-21_RNA-analysis/pipline_2023-3-21_RNA-analysis/6_featureCounts/h.genename.matrix',header=1)
# others<-read.table('/home/wus/2023_3_19-009-cxl_RNA/TPL202303714+TPL202303545/pipline_2023-3-23_RNA-analysis/6_featureCounts/h.genename.matrix',header=1)
# merge_h=merge(WT460,others)

# write.csv(merge_h,'/home/wus/2023_3_19-009-cxl_RNA/merge_file-6-11/merge_h.csv',row.names = FALSE)




rm(list=ls())
library(DESeq2) 
library(dplyr)

df<-read.table('/scratch/wus/030-2023-12-14-cxl-mus-liver/TPL2023113386/TPL2023113386/analysis/7_featureCounts/fc/m2',header=1,row.names=1)
#重排列
# m_base <- df %>% select(order(colnames(df)))

# m_base[1:4,]
colnames(df)
m_base <- df
# condition <- factor(c('Thr_m_HO16','Complete_m_HO28','Thr_m_HO29','Complete_m_HO30','Complete_m_HO32','Thr_m_HO48',
#      'Complete_m_HO50','Complete_m_HO85','Thr_m_HO86','Complete_m_HO87','Complete_m_HO88','Thr_m_HO9',
#      'Thr_m_WT24','Complete_m_WT33','Thr_m_WT46','Thr_m_WT49','Complete_m_WT51','Complete_m_WT54',
#      'Thr_m_WT77','Complete_m_WT79','Thr_m_WT80','Thr_m_WT81'))
colnames(m_base) <-c('Thr_m_HO16','Complete_m_HO28','Thr_m_HO29','Complete_m_HO30','Complete_m_HO32','Thr_m_HO48',
     'Complete_m_HO50','Complete_m_HO85','Thr_m_HO86','Thr_m_HO87','Complete_m_HO88','Thr_m_HO9',
     'Thr_m_WT24','Complete_m_WT33','Thr_m_WT46','Thr_m_WT49','Complete_m_WT51','Complete_m_WT54',
     'Thr_m_WT77','Complete_m_WT79','Thr_m_WT80','Thr_m_WT81')
m_base <- m_base %>% select(order(colnames(m_base)))
condition <- factor(c(rep("Complete_m_HO",6), rep("Complete_m_WT",4), rep("Thr_m_HO",6), rep("Thr_m_WT",6)) )
#condition <- factor(c(rep("A120",4), rep("A299",4), rep("A65",4), rep("E232",4), rep("E258",4), rep("E32",4),rep("E64",4), rep("WT",4))) 
# condition <- factor(colnames(df))
# condition <- factor(c(rep("HO",12), rep("WT",10)) )
# condition <- factor(c(rep("Complete_WT",4), rep("A299",4), rep("A65",4), rep("E232",4), rep("E258",4), rep("E32",4),rep("E64",4), rep("WT",4))) 
   
coldata <- data.frame(row.names = colnames(m_base), condition)
coldata 

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)


# keep <- rowSums(counts(dds) >= 4) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数  
# dds <- dds[keep, ] 

####跑一下层次聚类和pca
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

# rld <- rlog(dds, blind=FALSE)     #vst()标准化。
# plotPCA(rld, intgroup="condition",ntop=500)     
plotPCA(vsd, intgroup="condition",ntop=500)  


pcaData <- plotPCA(rld, intgroup=c("condition"),returnData = T)
#这里按照condition排序了，原因见下。
pcaData <- pcaData[order(pcaData$condition,decreasing=F),]
#知道每一个组有多少样本
table(pcaData$condition)
# II III  IV 
# 20  11  10 
#根据上面的结果，设置每一个组的数量，方便加颜色；这一步是把每一个样本根据分组情况画出来，效果见图1
# plot(pcaData[,1:2],pch = 19,col= c(rep("red",20),rep("green",11),rep("blue",10)),
#      cex=1)
plot(pcaData[,1:2])+
#加上样本名字，效果见图2
text(pcaData[,1],pcaData[,2],row.names(pcaData),cex=1, font = 1)
#加上图例，效果见图3
# legend(-70,43,inset = .02,pt.cex= 1.5,title = "Grade",legend = c("II", "III","IV"), 
#        col = c( "red","green","blue"),pch = 19, cex=0.75,bty="n")


#######################

coldata <- data.frame(row.names = colnames(m_base), condition)
coldata$condition <- relevel(coldata$condition, ref = "Complete_m_HO")
coldata

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

# keep <- rowSums(counts(dds) >= 4) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数

# # # 预处理，过滤低丰度的数据
# # countData <- count[apply(count, 1, sum) > 0 , ]

# dds <- dds[keep, ]

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
     write.table(resdata,file = paste('/scratch/wus/030-2023-12-14-cxl-mus-liver/TPL2023113386/TPL2023113386/analysis/8_DEseq2/',i,'.txt',sep=''),sep='\t') 
}

#######################

coldata <- data.frame(row.names = colnames(m_base), condition)
coldata$condition <- relevel(coldata$condition, ref = "Complete_m_WT")
coldata

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

# keep <- rowSums(counts(dds) >= 4) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数

# # # 预处理，过滤低丰度的数据
# # countData <- count[apply(count, 1, sum) > 0 , ]

# dds <- dds[keep, ]

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
     write.table(resdata,file = paste('/scratch/wus/030-2023-12-14-cxl-mus-liver/TPL2023113386/TPL2023113386/analysis/8_DEseq2/',i,'.txt',sep=''),sep='\t') 
}

####################### Thr_m_HO
coldata <- data.frame(row.names = colnames(m_base), condition)
coldata$condition <- relevel(coldata$condition, ref = "Thr_m_HO")
coldata

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

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
     write.table(resdata,file = paste('/scratch/wus/030-2023-12-14-cxl-mus-liver/TPL2023113386/TPL2023113386/analysis/8_DEseq2/',i,'.txt',sep=''),sep='\t') 
}
####################### Thr_m_WT
coldata <- data.frame(row.names = colnames(m_base), condition)
coldata$condition <- relevel(coldata$condition, ref = "Thr_m_WT")
coldata

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

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
     write.table(resdata,file = paste('/scratch/wus/030-2023-12-14-cxl-mus-liver/TPL2023113386/TPL2023113386/analysis/8_DEseq2/',i,'.txt',sep=''),sep='\t') 
}




rm(list=ls())
library(DESeq2) 
library(dplyr)

df<-read.table('/scratch/wus/028-2023-11-6_cxl_RNA_o22/7_featureCounts/m2',header=1,row.names=1)
#重排列
# m_base <- df %>% select(order(colnames(df)))

# m_base[1:4,]

m_base <- df
condition <- factor(c(rep("A_mut",12), rep("E_mut",16), rep("WT",4))) 

   
coldata <- data.frame(row.names = colnames(m_base), condition)
coldata 

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)


keep <- rowSums(counts(dds) >= 4) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数  
dds <- dds[keep, ] 




##############

coldata <- data.frame(row.names = colnames(m_base), condition)
coldata$condition <- relevel(coldata$condition, ref = "E_mut")
coldata

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

keep <- rowSums(counts(dds) >= 4) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数

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
     write.table(resdata,file = paste('/scratch/wus/028-2023-11-6_cxl_RNA_o22/8_DEseq2/',i,'.txt',sep=''),sep='\t') 
}

