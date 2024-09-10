
bio<-read.table('../src/Homo_sapiens.GRCh38.84_gene_name-gene_biotype.txt',header=1)

df<-read.csv('../../7_DESeq2/h_WT_460_VS_KR.csv')
df2<-merge(bio,df,by.y='Row.names',by.x='gene_name')
write.table(df2,'h_WT_460_VS_KR-biotype.txt',sep='\t')



df<-read.csv('../../7_DESeq2/h_WT_VS_KR.csv')
df2<-merge(bio,df,by.y='Row.names',by.x='gene_name')
write.table(df2,'h_WT_VS_KR-biotype.txt',sep='\t')






df<-read.csv('../../7_DESeq2/h_WT_450_VS_KR.csv')
df2<-merge(bio,df,by.y='Row.names',by.x='gene_name')
write.table(df2,'h_WT_450_VS_KR-biotype.txt',sep='\t')




df<-read.csv('../../7_DESeq2/h_WT&WT460_VS_KR.csv')
df2<-merge(bio,df,by.y='Row.names',by.x='gene_name')
write.table(df2,'h_WT&WT460_VS_KR.csv.txt',sep='\t')

