##############
#1.2	
# 58831	27

#cbind

# Geneid	Description
#############
setwd('/home/wus/023_2023-9-7_cxl_RNA/023-clx-X101SC23062450-Z01-J004/cell_line_analysis/8_DEseq2/')
df<-read.table('../7_featureCounts/h2v2',header=1)
head(df)

nrow(df)
ncol(df)
df$Description<-'NNN'
test2 = df[,c(1,ncol(df),3:ncol(df)-1)]
head(test2)
tail(test2)
#### 
setwd('/home/wus/023_2023-9-7_cxl_RNA/023-clx-X101SC23062450-Z01-J004/cell_line_analysis/11_GSEA/')

V1 <- c('#1.2', nrow(df)) 
V2 <- c('', ncol(df)-2) 
EmpData <- data.frame(V1=V1, V2=V2)
print(EmpData)
write.table(EmpData,'express.gct',row.names=F,sep='\t',col.names =F,quote=FALSE)


write.table(test2,'express.gct',row.names=F,append=T,sep='\t',quote=FALSE)

#######   306_EV
setwd('/home/wus/023_2023-9-7_cxl_RNA/023-clx-X101SC23062450-Z01-J004/cell_line_analysis/8_DEseq2/')
df<-read.table('../7_featureCounts/h2v2',header=1)
head(df)
df2<-df[,!grepl('WT',colnames(df))]
df2$Description<-'NNN'
test2 = df2[,c(1,ncol(df2),3:ncol(df2)-1)]
head(test2)

setwd('/home/wus/023_2023-9-7_cxl_RNA/023-clx-X101SC23062450-Z01-J004/cell_line_analysis/11_GSEA/')
V1 <- c('#1.2', nrow(df2)) 
V2 <- c('', ncol(df2)-2) 
EmpData <- data.frame(V1=V1, V2=V2)
print(EmpData)
write.table(EmpData,'306_EV_express.gct',row.names=F,sep='\t',col.names =F,quote=FALSE)

write.table(test2,'306_EV_express.gct',row.names=F,append=T,sep='\t',quote=FALSE)

### 306_WT
setwd('/home/wus/023_2023-9-7_cxl_RNA/023-clx-X101SC23062450-Z01-J004/cell_line_analysis/8_DEseq2/')
df<-read.table('../7_featureCounts/h2v2',header=1)
head(df)
df2<-df[,!grepl('EV',colnames(df))]
df2$Description<-'NNN'
test2 = df2[,c(1,ncol(df2),3:ncol(df2)-1)]
head(test2)

setwd('/home/wus/023_2023-9-7_cxl_RNA/023-clx-X101SC23062450-Z01-J004/cell_line_analysis/11_GSEA/')
V1 <- c('#1.2', nrow(df2)) 
V2 <- c('', ncol(df2)-2) 
EmpData <- data.frame(V1=V1, V2=V2)
print(EmpData)
write.table(EmpData,'306_WT_express.gct',row.names=F,sep='\t',col.names =F,quote=FALSE)

write.table(test2,'306_WT_express.gct',row.names=F,append=T,sep='\t',quote=FALSE)

### EV_WT
setwd('/home/wus/023_2023-9-7_cxl_RNA/023-clx-X101SC23062450-Z01-J004/cell_line_analysis/8_DEseq2/')
df<-read.table('../7_featureCounts/h2v2',header=1)
head(df)
df2<-df[,!grepl('306',colnames(df))]
df2$Description<-'NNN'
test2 = df2[,c(1,ncol(df2),3:ncol(df2)-1)]
head(test2)

setwd('/home/wus/023_2023-9-7_cxl_RNA/023-clx-X101SC23062450-Z01-J004/cell_line_analysis/11_GSEA/')
V1 <- c('#1.2', nrow(df2)) 
V2 <- c('', ncol(df2)-2) 
EmpData <- data.frame(V1=V1, V2=V2)
print(EmpData)
write.table(EmpData,'EV_WT_express.gct',row.names=F,sep='\t',col.names =F,quote=FALSE)

write.table(test2,'EV_WT_express.gct',row.names=F,append=T,sep='\t',quote=FALSE)