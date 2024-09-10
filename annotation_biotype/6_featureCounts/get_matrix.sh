#awk ' {$2=$3=$4=$5=$6 = ""; print $0}' h.genename  > h.genename.matrix
#sed -i '1d' h.genename.matrix
#sed '1s/\/home\/wus\/2023_3_19-009-cxl_RNA\/TPL202303714+TPL202303545\/pipline_2023-3-23_RNA-analysis\/5_XenofilteR\/human\///g' h.genename.matrix  > 1
#sed '1s/\.sorted_Filtered\.bam//g' 1 > 2


awk ' {$2=$3=$4=$5=$6 = ""; print $0}' m.genename  > m.genename.matrix
sed -i '1d' m.genename.matrix
sed '1s/\/home\/wus\/2023_3_19-009-cxl_RNA\/TPL202303714+TPL202303545\/pipline_2023-3-23_RNA-analysis\/5_XenofilteR\/mouse\/Filtered_bams\///g' m.genename.matrix  > m1
sed '1s/\.sorted_Filtered\.bam//g' m1 > m2

