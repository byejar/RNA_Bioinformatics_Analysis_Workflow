#base
#path 
Path='/home/wus/017-2023-6-28_cxl/dish/Data/X101SC23062450-Z01-J001/00.CleanData'
Path2='/home/wus/017-2023-6-28_cxl/dish/Data/X101SC23062450-Z01-J002/00.CleanData'
#database
#m_Data='/home/wus/scRNA_program_blackmouseRNA/Nude/hisat2_database/genome'
#h_Data='/home/wus/grch38/genome'

echo '>>>>> fastqc'

for i in `cat file_name.txt `
do

	nohup fastqc -t 4  $Path/$i/$i\_*.fq.gz -o ./fastQC &
        #nohup hisat2 -p 10  -x $h_Data  -1 $Path2/$i/$i\_1.clean.fq.gz  -2 $Path2/$i/$i\_1.clean.fq.gz -S h_$i\.sam &

done

echo '>>>>> wait'

wait

echo '>>>>> multiqc'
cd ./fastQC
multiqc .
