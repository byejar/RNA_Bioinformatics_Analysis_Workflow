
bio<-read.table('/home/wus/2023_3_19-009-cxl_RNA/merge_file-6-11/WT460_WT/src/Homo_sapiens.GRCh38.84_gene_name-gene_biotype.txt',header=1)
#df<-read.table('/scratch/wus/028-2023-11-6_cxl_RNA_o22/8_DEseq2/condition_A_mut_vs_WT.txt')

#df2<-merge(df,bio,by.x='Row.names',by.y='gene_name')

#write.table(df2,'A_mut_vs_WT-biotype.txt',sep='\t')



df<-read.table('/scratch/wus/028-2023-11-6_cxl_RNA_o22/8_DEseq2/RRA_Amut_1WT_dw.txt')

df2<-merge(df,bio,by.x='Name',by.y='gene_name',all.y=F)

write.table(df2,'RRA_Amut_1WT_dw-biotype.txt',sep='\t')


df<-read.table('/scratch/wus/028-2023-11-6_cxl_RNA_o22/8_DEseq2/RRA_Amut_1WT_up.txt')

df2<-merge(df,bio,by.x='Name',by.y='gene_name',all.y=F)

write.table(df2,'RRA_Amut_1WT_up-biotype.txt',sep='\t')




#df<-read.table('../h_WT460-and-WT_VS_KR(delete_KR26_3).txt')

#df2<-merge(df,bio,by.x='Row.names',by.y='gene_name')

#write.table(df2,'h_WT460-and-WT_VS_KR(delete_KR26_3)-biotype.txt',sep='\t')
