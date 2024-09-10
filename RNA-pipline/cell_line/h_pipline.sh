###
# time: 2023-7-20
# wus
###

#### path (change)
wd_path=`pwd`
################ change1 input filefold pathway
fq_file_pwd='/home/wus/018-2023-7-18_cxl_MC38/X101SC23062450-Z01-J004/00.CleanData'
#####


#### database
#m_Data='/home/wus/scRNA_program_blackmouseRNA/Nude/hisat2_database/genome'
h_Data='/home/wus/grch38/genome'
#gff
H_gff='/home/wus/database/Homo_sapiens.GRCh38.84.gtf'
#M_gff='/home/wus/cellranger_database/refdata-gex-mm10-2020-A/genes/genes.gtf'


#### file
#           $fq_file_pwd/$i/$i\_1.fq.gz
#           $fq_file_pwd/$i/$i\_2.fq.gz
####


#### write file
ls $fq_file_pwd  > ./2_FastQC_MultiQC/file_name.txt 
echo 'write file_name.txt'
####

source activate DESeq2


echo '>>>>>>>>>>>>>>>>>>>>> start '
### 
cd ./2_FastQC_MultiQC/
#./fastQC.sh

echo '>>>>>>>>> 2.fastqc'

for i in `cat file_name.txt `
do

	nohup fastqc -t 4  $fq_file_pwd/$i/*.fq.gz -o ./fastQC &

done

echo '>>>>> wait'

wait

echo '>>>>> multiqc'
cd ./fastQC
multiqc .

cd ../..
##

##### 
#conda activate base
echo '>>>>>>>>> 4.hisat2'

cd ./4_hisat2/

#################### change2 fastq file subffix

for i in `cat ../2_FastQC_MultiQC/file_name.txt `

do

	#nohup hisat2 -p 10  -x $m_Data  -1 $fq_file_pwd/$i/$i\_1.clean.fq.gz  -2 $fq_file_pwd/$i/$i\_2.clean.fq.gz -S m_$i\.sam &
        nohup hisat2 -p 10  -x $h_Data  -1 $fq_file_pwd/$i/$i\_1.clean.fq.gz  -2 $fq_file_pwd/$i/$i\_2.clean.fq.gz -S h_$i\.sam &

done

cd ..

echo '>>>>> wait'
wait

####

echo '>>>>>>>>> 5.samtools'

cd ./5_samtools/

ls ../4_hisat2/ | egrep *.sam  > sam.list

for i in `cat sam.list`
do
        j=${i%.*}
        nohup samtools view -@ 5  -bS ../4_hisat2/$i > $j\.bam &
done

echo '>>>>> wait'
wait


for i in `cat sam.list`
do
        j=${i%.*}
        nohup samtools sort   $j\.bam -@ 5 -m 5G -o $j\.sorted.bam &
done

cd ..

echo '>>>>> wait'
wait



####

source activate R3.7
#echo '>>>>>>>>> 6_XenofilteR'
#cd ./6_XenofilteR

## file
#ls $wd_path/5_samtools/* | grep h_ | grep sorted.bam$ > h_sorted.bam.txt
#ls $wd_path/5_samtools/* | grep m_ | grep sorted.bam$ > m_sorted.bam.txt

#paste h_sorted.bam.txt m_sorted.bam.txt -d, >  h_sample.list
#paste m_sorted.bam.txt h_sorted.bam.txt -d, >  m_sample.list


#mkdir mouse
#mkdir human

#nohup Rscript h_XenofilteR.R &
#nohup Rscript m_XenofilteR.R &

#echo '>>>>> wait'
#wait

#cd ..


###
echo '>>>>>>>>> 7_featureCounts'
cd ./7_featureCounts

#H_path='../6_XenofilteR/human/Filtered_bams/'
#M_path='../6_XenofilteR/mouse/Filtered_bams/'
H_path='../5_samtools/'
#M_path='../5_samtools/'

#featureCounts -T 32 -p -t exon -g gene_name -a ${M_gff} -o m.genename ${M_path}*.sorted.bam
featureCounts -T 32 -p -t exon -g gene_name -a ${H_gff} -o h.genename ${H_path}*.sorted.bam


H_path_trans=${H_path//\//\\\/}
#M_path_trans=${M_path//\//\\\/}



awk ' {$2=$3=$4=$5=$6 = ""; print $0}' h.genename  > h.genename.matrix
sed -i '1d' h.genename.matrix
sed "1s/${H_path_trans}//g" h.genename.matrix  > h1
sed '1s/.sorted_Filtered.bam//g' h1 > h2


#awk ' {$2=$3=$4=$5=$6 = ""; print $0}' m.genename  > m.genename.matrix
#sed -i '1d' m.genename.matrix
#sed "1s/${M_path_trans}//g" m.genename.matrix  > m1
#sed '1s/.sorted_Filtered.bam//g' m1 > m2


cd ..


#### DESeq2

