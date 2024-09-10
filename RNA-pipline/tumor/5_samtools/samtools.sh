ls ../4_hisat2/ | grep P  > sam.list

for i in `cat sam.list`
do
        j=${i%.*}
        #nohup samtools view -@ 5  -bS ../4_hisat2/$i > $j\.bam &
        nohup samtools sort   $j\.bam -@ 5 -m 5G -o $j\.sorted.bam &
done

#for i in `cat sam.list`
#do
#        j=${i%.*}
        #nohup samtools view -@ 5  -bS ../4_hisat2/$i > $j\.bam &
#        nohup samtools sort   $j\.bam -@ 5 -m 5G -o $j\.sorted.bam &
#done










#samtools view -@ 10  -bS ../4_hisat2/h_C1.sam > h_C1.bam
#samtools sort   h_C1.bam -@ 10 -m 20G -o h_C1.sorted.bam

