### 合并看
df<-read.table('../7_featureCounts/m2',header=1,row.names=1)
df2<-df
m_base <- df2 %>% select(order(colnames(df2)))
m_base[1:4,]

condition <- factor(c(rep("KR293",4), rep("KR83",4), rep("parental",4), rep("WT",4)))
condition <- factor(c(rep("KR",8),  rep("WT",8)))   
coldata <- data.frame(row.names = colnames(m_base), condition)

coldata

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)

keep <- rowSums(counts(dds) >= 10) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数
dds <- dds[keep, ]

dds1 <- DESeq(dds)    # 做DESeq分析，输入的只能是原始矩阵，虽然有多种标准化方法，但是DESeq2在做差异表达分析的时候会用它自己的方法进行标准化。估计方差稳定变换只是为了做聚类或PCA
resultsNames(dds1)    # 查看结果的名称。
dds1$condition        #默认后者的处理组比前面的对照组。
res <- results(dds1)  # 必要步骤！！！
summary(res)          #看一下结果的概要信息，p值默认小于0.1。

table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。
res <- res[order(res$padj),]  #按照padj 进行升序排列
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds1, normalized=TRUE)),by="row.names",sort=FALSE)#保存时保存的是count标准化结果
write.table(resdata,file = "m_WT&parental_VS_KR.txt",sep='\t') 
