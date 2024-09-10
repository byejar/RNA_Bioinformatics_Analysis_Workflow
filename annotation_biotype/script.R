
bio<-read.table('../src/Homo_sapiens.GRCh38.84_gene_name-gene_biotype.txt',header=1)
df<-read.table('../h_WT460_VS_WT.txt')

df2<-merge(df,bio,by.x='Row.names',by.y='gene_name')

write.table(df2,'h_WT460_VS_WT-biotype.txt',sep='\t')



df<-read.table('../h_WT460_VS_KR(delete_KR26_3).txt')

df2<-merge(df,bio,by.x='Row.names',by.y='gene_name')

write.table(df2,'h_WT460_VS_KR(delete_KR26_3)-biotype.txt',sep='\t')






df<-read.table('../h_WT460-and-WT_VS_KR(delete_KR26_3).txt')

df2<-merge(df,bio,by.x='Row.names',by.y='gene_name')

write.table(df2,'h_WT460-and-WT_VS_KR(delete_KR26_3)-biotype.txt',sep='\t')
