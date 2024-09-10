
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

df<-read.table('/scratch/wus/028-2023-11-6_cxl_RNA_o22/7_featureCounts/m2',header=1,row.names=1)
#重排列
# m_base <- df %>% select(order(colnames(df)))

# m_base[1:4,]

m_base <- df
condition <- factor(c(rep("A120",4), rep("A299",4), rep("A65",4), rep("E232",4), rep("E258",4), rep("E32",4),rep("E64",4), rep("WT",4))) 

   
coldata <- data.frame(row.names = colnames(m_base), condition)
coldata 

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)


keep <- rowSums(counts(dds) >= 4) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数  
dds <- dds[keep, ] 

####跑一下层次聚类和pca
library('factoextra')          
vsd <- vst(dds,blind = TRUE)    #参数blind=TRUE是为了不把样本分组信息考虑在内——即以无偏的方式进行转换
sampleDists <- dist(t(assay(vsd))) 
res1 <- hcut(sampleDists, k = 8, stand = FALSE,hc_method ="average" ) 
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


#######################

coldata <- data.frame(row.names = colnames(m_base), condition)
coldata$condition <- relevel(coldata$condition, ref = "WT")
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
