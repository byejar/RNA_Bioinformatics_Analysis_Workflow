#R3.7
#conda activate R3.7
#path
H_path='../4-1_samtools/'
#M_path='/home/wus/12-15_004_zc_TPL2022112807_RNA/TPL2022112807/analysis/4-1_samtools/'

#gff
H_gff='/home/wus/database/Homo_sapiens.GRCh38.84.gtf'
#M_gff='/home/wus/cellranger_database/refdata-gex-mm10-2020-A/genes/genes.gtf'





#featureCounts -T 32 -p -t exon -g gene_name -a ${M_gff} -o m.genename ${M_path}*.sorted.bam
featureCounts -T 32 -p -t exon -g gene_name -a ${H_gff} -o h.genename ${H_path}*.sorted.bam
