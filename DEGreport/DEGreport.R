library(DEGreport)
library(DESeq2) 
library(dplyr)

df<-read.csv('/home/wus/2023_3_19-009-cxl_RNA/merge_file-6-11/3.csv',header=1,row.names=1)
#重排列


df2<-df[,!grepl('WT_450|WT_460|KR_26_3',colnames(df))]

m_base <- df2 %>% select(order(colnames(df2)))
m_base[1:4,]

condition <- factor(c(rep("KR234",4), rep("KR26",3), rep("KR274",4), rep("KR365",4), rep("KR692",4), rep("KR76",4), rep("WT",4))) 

   
coldata <- data.frame(row.names = colnames(m_base), condition)
coldata 

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)


keep <- rowSums(counts(dds) >= 10) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数  
dds <- dds[keep, ] 



# dds <- DESeqDataSetFromTximport(m_base, colData=coldata, design=~condition)

# dds <- DESeq(dds)
dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)
# ddsTC <- DESeqDataSet(dds, ~ condition  )

# #作为比较的模型，以 1 作为对照
# ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ 'WT')
# resTC <- results(ddsTC)
# resTC$symbol <- mcols(ddsTC)$symbol
# head(resTC[order(resTC$padj),], 4)

rld <- rlog(dds, blind=TRUE)
rld_mat <- assay(rld)

degPatterns(rld_mat , metadata = coldata, time ="condition", col=NULL)
# rld <- vst(dds, blind=FALSE)     #vst()标准化。
# dds_rlg<-rlog(dds)
# ma = assay(rlog(dds)) 

# clusters <- degPatterns(ma, metadata = meta, time ="condition", col=NULL)

 
# res <- degPatterns(dds_rlg, condition, time = "condition")

# ddsTC <- DESeqDataSet(fission, ~ minute )

# #作为比较的模型，以 1 作为对照
# ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ 1)
# resTC <- results(ddsTC)
# resTC$symbol <- mcols(ddsTC)$symbol
# head(resTC[order(resTC$padj),], 4)


# Extract results
res_LRT <- results(dds_lrt)



# Subset the LRT results to return genes with padj < 0.05
padj.cutoff=0.00001
sig_res_LRT <- res_LRT %>%
               data.frame() %>%
               tibble::rownames_to_column(var="gene") %>% 
               as_tibble() %>% 
               filter(padj < padj.cutoff)
 
# Get sig gene lists
sigLRT_genes <- sig_res_LRT %>% 
                pull(gene)
                
length(sigLRT_genes)

# Compare to numbers we had from Wald test
# nrow(sigOE)
# nrow(sigKD)

rld <- rlog(dds, blind=TRUE)
rld_mat <- assay(rld)


clustering_sig_genes <- sig_res_LRT %>%
                  arrange(padj) %>%
                  head(n=1000)


# Obtain rlog values for those significant genes
cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]

# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
degPatterns(rld_mat , metadata = coldata, time ="condition", col=NULL)
clusters <- degPatterns(cluster_rlog, metadata = coldata, time ="condition", col=NULL)


################# A WT 
library(DEGreport)
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


dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)


rld <- rlog(dds, blind=TRUE)
rld_mat <- assay(rld)

res_LRT <- results(dds_lrt)



# Subset the LRT results to return genes with padj < 0.05
padj.cutoff=0.00001
sig_res_LRT <- res_LRT %>%
               data.frame() %>%
               tibble::rownames_to_column(var="gene") %>% 
               as_tibble() %>% 
               filter(padj < padj.cutoff)
 
# Get sig gene lists
sigLRT_genes <- sig_res_LRT %>% 
                pull(gene)
                
length(sigLRT_genes)

# Compare to numbers we had from Wald test
# nrow(sigOE)
# nrow(sigKD)

rld <- rlog(dds, blind=TRUE)
rld_mat <- assay(rld)


clustering_sig_genes <- sig_res_LRT %>%
                  arrange(padj) %>%
                  head(n=1000)


# Obtain rlog values for those significant genes
cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]

# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups

clusters <- degPatterns(cluster_rlog, metadata = coldata, time ="condition", col=NULL)



################### 12-4
df<-read.table('/scratch/wus/028-2023-11-6_cxl_RNA_o22/7_featureCounts/m2',header=1,row.names=1)
#重排列
# m_base <- df %>% select(order(colnames(df)))

# m_base[1:4,]

df2<-df[,!grepl('E',colnames(df))]

m_base <- df2 %>% select(order(colnames(df2)))


#condition <- factor(c(rep("A120",4), rep("A299",4), rep("A65",4), rep("E232",4), rep("E258",4), rep("E32",4),rep("E64",4), rep("WT",4))) 
condition <- factor(c(rep("A",12),  rep("WT",4))) 

   
coldata <- data.frame(row.names = colnames(m_base), condition)
coldata 

dds <- DESeqDataSetFromMatrix(countData=m_base, colData=coldata, design=~condition)


keep <- rowSums(counts(dds) >= 4) >= 10  #过滤低表达基因，至少有2个样品都满足10个以上的reads数  
dds <- dds[keep, ] 
