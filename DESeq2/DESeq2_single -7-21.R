library(DESeq2)
library(dplyr)


setwd('/home/wus/018-2023-7-18_cxl_MC38/7-18_RNA-pipline/tumor/')
df<-read.table('./7_featureCounts/m2',header=1,row.names=1)
df2<-df

m_base <- df2 %>% select(order(colnames(df2)))
m_base[1:4,]

condition <- factor(c(rep("KR293",4), rep("KR83",4), rep("parental",4), rep("WT",4)))
coldata <- data.frame(row.names = colnames(m_base), condition)
coldata$condition <- relevel(coldata$condition, ref = "parental")
coldata

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

keep <- rowSums(counts(dds) >= 10) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数
dds <- dds[keep, ]

dds <- DESeq(dds, betaPrior = FALSE)
resultsNames(dds)

### 改变输出的顺序
results(dds, contrast=c("condition","KR293", "parental"))
results(dds, contrast=c("condition", "parental","KR293"))
### 不方便，全改成处理组 vs con ？
### 所以这里的名字都应该改成KR vs con 

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
for (i in resultsNames(dds)[c(2,3)]){
     df<-read.table(paste(i,'.txt',sep=''),header=1,row.names=1)
     df[1:4,1:4]
     a<-filter(df, (log2FoldChange > 1 | log2FoldChange < -1) & (pvalue<0.05))[,c("Row.names","log2FoldChange","pvalue")]
     colnames(a)<-c("Row.names",paste(i,'log2FC',sep='_'),paste(i,'pvalue',sep='_'))     
          if(n==0){m=a}
          else{m<-merge(m,a,by='Row.names')}
     n=n+1

}

write.table(m,file = 'MC38_tumor_parental_VS_KR.txt',sep='\t',row.names=F) 


n=0
for (i in resultsNames(dds)[c(2,3)]){
     df<-read.table(paste(i,'.txt',sep=''),header=1,row.names=1)
     df[1:4,1:4]
     a<-filter(df, (log2FoldChange > -1 | log2FoldChange < 1) | (pvalue > 0.05))[,c("Row.names","log2FoldChange","pvalue")]
     colnames(a)<-c("Row.names",paste(i,'log2FC',sep='_'),paste(i,'pvalue',sep='_'))     
          if(n==0){m=a}
          else{m<-merge(m,a,by='Row.names')}
     n=n+1

}

write.table(m,file = 'MC38_tumor_parental_VS_KR_non-change.txt',sep='\t',row.names=F) 



########################
df<-read.table('../7_featureCounts/m2',header=1,row.names=1)
df2<-df

m_base <- df2 %>% select(order(colnames(df2)))
m_base[1:4,]

condition <- factor(c(rep("KR293",4), rep("KR83",4), rep("parental",4), rep("WT",4)))
coldata <- data.frame(row.names = colnames(m_base), condition)
coldata$condition <- relevel(coldata$condition, ref = "WT")
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
for (i in resultsNames(dds)[c(2,3)]){
     df<-read.table(paste(i,'.txt',sep=''),header=1,row.names=1)
     df[1:4,1:4]
     a<-filter(df, (log2FoldChange > 1 | log2FoldChange < -1) & (pvalue<0.05))[,c("Row.names","log2FoldChange","pvalue")]
     colnames(a)<-c("Row.names",paste(i,'log2FC',sep='_'),paste(i,'pvalue',sep='_'))     
          if(n==0){m=a}
          else{m<-merge(m,a,by='Row.names')}
     n=n+1

}


write.table(m,file = 'MC38_tumor_WT_VS_KR.txt',sep='\t',row.names=F) 


n=0
for (i in resultsNames(dds)[c(2,3)]){
     df<-read.table(paste(i,'.txt',sep=''),header=1,row.names=1)
     df[1:4,1:4]
     a<-filter(df, (log2FoldChange > -1 | log2FoldChange < 1) | (pvalue > 0.05))[,c("Row.names","log2FoldChange","pvalue")]
     colnames(a)<-c("Row.names",paste(i,'log2FC',sep='_'),paste(i,'pvalue',sep='_'))     
          if(n==0){m=a}
          else{m<-merge(m,a,by='Row.names')}
     n=n+1

}

write.table(m,file = 'MC38_tumor_WT_VS_KR_non-change.txt',sep='\t',row.names=F) 


############# 合并 wt parental
## 取变化方向相同的
WT_VS_KR<-read.table('MC38_tumor_WT_VS_KR.txt',header=1,row.names=1)
head(WT_VS_KR)
WT_VS_KR<-WT_VS_KR[,!grepl('pvalue',colnames(WT_VS_KR))]
WT_VS_KR$sum <- rowSums(sign(WT_VS_KR))
WT_VS_KR<-WT_VS_KR[WT_VS_KR$sum != 0,]
WT_VS_KR<-WT_VS_KR[,1:2]
write.table(WT_VS_KR,file = 'MC38_tumor_WT_VS_KR-consistent.txt',sep='\t',row.names=T) 

#### 
parental_VS_KR<-read.table('MC38_tumor_parental_VS_KR.txt',header=1,row.names=1)
head(parental_VS_KR)
parental_VS_KR<-parental_VS_KR[,!grepl('pvalue',colnames(parental_VS_KR))]
parental_VS_KR$sum <- rowSums(sign(parental_VS_KR))
parental_VS_KR<-parental_VS_KR[parental_VS_KR$sum != 0,]
parental_VS_KR<-parental_VS_KR[,1:2]
head(parental_VS_KR)
write.table(parental_VS_KR,file = 'MC38_tumor_parental_VS_KR-consistent.txt',sep='\t',row.names=T) 

####
m<-merge(WT_VS_KR,parental_VS_KR, by='row.names',all=TRUE)
write.table(m,file = 'MC38_tumor_WT&parental_VS_KR-consistent.txt',sep='\t',row.names=F) 

##### 添加注释
bio<-read.table('../9_merge_file/src/Homo_sapiens.GRCh38.84_gene_name-gene_biotype.txt',header=1)
df2<-merge(bio,m,by.y='Row.names',by.x='gene_name')
write.table(df2,'MC38_tumor_WT&parental_VS_KR-consistent-biotype.txt',sep='\t')



