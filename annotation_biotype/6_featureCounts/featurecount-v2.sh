#R3.7
#path
#H_path='/home/wus/2023_3_19-009-cxl_RNA/TPL202303714+TPL202303545/pipline_2023-3-23_RNA-analysis/5_XenofilteR/human/'
M_path='/home/wus/2023_3_19-009-cxl_RNA/TPL202303714+TPL202303545/pipline_2023-3-23_RNA-analysis/5_XenofilteR/mouse/Filtered_bams/'

#gff
#H_gff='/home/wus/database/Homo_sapiens.GRCh38.84.gtf'
M_gff='/home/wus/cellranger_database/refdata-gex-mm10-2020-A/genes/genes.gtf'



#featureCounts -T 16 -f -t gene -g gene_id -a ${H_gff} -o h.count ${H_path}h_C1.sorted_Filtered.bam ${H_path}h_C2.sorted_Filtered.bam  ${H_path}h_C3.sorted_Filtered.bam ${H_path}h_T1.sorted_Filtered.bam ${H_path}h_T2.sorted_Filtered.bam ${H_path}h_T3.sorted_Filtered.bam


#featureCounts -T 16 -f -t gene -g gene_id -a ${M_gff} -o m.count ${M_path}234_1_clean.sorted.bam  ${M_path}234_2_clean.sorted.bam ${M_path}234_3_clean.sorted.bam ${M_path}234_4_clean.sorted.bam  ${M_path}WT_1_clean.sorted.bam ${M_path}WT_2_clean.sorted.bam ${M_path}WT_3_clean.sorted.bam ${M_path}WT_4_clean.sorted.bam



#featureCounts -T 32 -p -t exon -g gene_id -a ${M_gff} -o m.geneid ${M_path}*.sorted.bam

featureCounts -T 32 -p -t exon -g gene_name -a ${M_gff} -o m.genename ${M_path}*.sorted_Filtered.bam
#featureCounts -T 16 -p -t exon -g gene_name -a ${H_gff} -o h.genename ${H_path}*.sorted_Filtered.bam
