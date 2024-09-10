#awk -OSF '\t' '{ $2=$3=$4=$5=$6=null;print }' m.genename  > m.genename.matrix
awk ' {$2=$3=$4=$5=$6 = ""; print $0}' h.genename  > h.genename.matrix
sed -i '1d' h.genename.matrix

sed '1s/..\/4-1_samtools\///g' h.genename.matrix  > h1
sed '1s/.sorted.bam//g' h1 > h2
