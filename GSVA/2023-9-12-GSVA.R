


gmt<-gmtPathways("/home/wus/019_2023-07-20_zc_BMDM/202307BMDM-si/cell_line/10_analysis/HALLMARK_MTORC1_SIGNALING.v2023.1.Hs.gmt")

######## GO KEGG
conda activate irGSEA
# library(Seurat)
# library(msigdbr）
library(GSVA)
# library(tidyverse)
# library(clusterProfiler)
library(patchwork)
library(limma)
library(ggplot2)
library(dplyr)
# rm(list=1s())

# scRNA_aneuploid<-subset(scRNA_gT,subset=idents=='aneuploid')

# DefaultAssay(scRNA_aneuploid)<-"RNA"
# scRNA_aneuploid <-NormalizeData(scRNA_aneuploid)

# save(scRNA_aneuploid,file='/home/wus/021-wus-USP39-mTOR/scRNA_aneuploid.RData')

# load('/home/wus/021-wus-USP39-mTOR/scRNA_aneuploid.RData')

setwd('/home/wus/023_2023-9-7_cxl_RNA/023-clx-X101SC23062450-Z01-J004/cell_line_analysis/8_DEseq2/')
df<-read.table('../7_featureCounts/h2v2',header=1,row.names=1)


df2<-df[,!grepl('306',colnames(df))]
m_base <- df2 %>% select(order(colnames(df2)))
m_base[1:4,]

# head(scRNA_aneuploid)

# scRNA_aneuploid$USP39_ident <-paste(scRNA_aneuploid$USP39_express, scRNA_aneuploid$orig.ident, sep="_") 
# Idents(scRNA_aneuploid)<-"USP39_ident"
# expr <-AverageExpression(scRNA_aneuploid,assays= "RNA",slot= "data")[[1]]
# expr<-expr[rowSums(expr)>0,]#选取非零基因
# expr <-as.matrix(expr)
# head(expr)
expr <-as.matrix(m_base)
head(expr)
expr<-expr[rowSums(expr)>0,]#选取非零基因
head(expr)

gmt <- GSEABase::getGmt("/home/wus/MSigDB/msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/h.all.v2023.1.Hs.symbols.gmt")

gsva.res <-gsva(expr,gmt,method="gsva",kcdf = 'Poisson')
# saveRDS(gsva.res,"gsva.res.rds")
gsva.df <-data.frame(Genesets=rownames(gsva.res),gsva.res,check.names = F)
# write.csv(gsva.df,"gsva_res.csv",row.names F)
pheatmap::pheatmap(gsva.res,show_colnames =T,
scale ="row",angle_col ="45",cluster_cols=F,
color= colorRampPalette(c("navy","white","firebrick3"))(50))



sort(colnames(gsva.df)[-1])

group_list <-data.frame(sample =sort(colnames(gsva.df)[-1]),group =c(rep("EV",4),rep("WT",4)))
group_list

design <-model.matrix(~0 +factor(group_list$group))
colnames(design)<-levels(factor(group_list$group))
rownames(design)<-colnames(gsva.res)
# rownames(design)<-c('aneuploid','aneuploid')
design


#构建差异比较矩阵
contrast.matrix <-makeContrasts(EV-WT,levels = design)
#差异分析，case vs.con
fit <-lmFit(gsva.res,design)
fit2 <-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
x <-topTable(fit2,coef =1,n =Inf,adjust.method ="BH",sort.by ="P")
head(x)
#

# 把通路的1imma分析结果保存到文件
# write.csv(x,"gsva_limma.csv",quote F)
# 输出t值，用做绘图的输入数据
# pathway <-str_replace(row.names(x),"HALLMARK_","")
pathway <-row.names(x)
df <-data.frame(ID =pathway,score =x$t)
head(df)

# write.csv(df,"enrich_bar.csv",quote F,row.names F)

# 按照scoref的值分组
cutoff <-0
df$group <-cut(df$score,breaks =c(-Inf,cutoff,Inf),labels =c(1,2))
head(df)

# 按照score排序
sortdf <-df[order(df$score),]
sortdf$ID <-factor(sortdf$ID,levels =sortdf$ID)
head(sortdf)
#

ggplot(sortdf,aes(ID,score,fill =group))+geom_bar(stat ='identity')+
coord_flip()+
scale_fill_manual(values =c('palegreen3','dodgerblue4'),guide= 'none')+
#画2条虚线
geom_hline(yintercept= c(-1,1),
color="white",
linetype=2,#画虚线
linewidth=0.3)+#线的粗细
#写1abel
geom_text(data= subset(df,score > 0),
aes(x=ID,y=-0.01,label=ID,color =group),
size=3,#字的大小
hjust="outward")+#字的对齐方式
geom_text(data =subset(df,score <0),
aes(x=ID,y=0.01,label=paste0("",ID),color=group),#bar跟坐标轴间留出间隙
size =3,hjust ="inward")+
scale_colour_manual(values =c("black","black"),guide= FALSE)+
xlab("")+ylab("t value of GSVA score")+
theme_bw()+#去除背景色
theme(panel.grid=element_blank())+#去除网格线
theme(panel.border=element_rect(size=0.6))+#边框粗细
theme(axis.line.y =element_blank(),axis.ticks.y =element_blank(),axis.text.y =element_blank())
# ggsave("gsva.pdf",width 6,height 8)



########################## kegg
gmt <- GSEABase::getGmt("/home/wus/MSigDB/msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
gmt <- GSEABase::getGmt("/home/wus/MSigDB/msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/c5.go.bp.v2023.1.Hs.symbols.gmt")

gsva.res <-gsva(expr,gmt,method="gsva",kcdf = 'Poisson')
# saveRDS(gsva.res,"gsva.res.rds")
gsva.df <-data.frame(Genesets=rownames(gsva.res),gsva.res,check.names = F)
# write.csv(gsva.df,"gsva_res.csv",row.names F)
# pheatmap::pheatmap(gsva.res,show_colnames =T,
# scale ="row",angle_col ="45",cluster_cols=F,
# color= colorRampPalette(c("navy","white","firebrick3"))(50))

pheatmap::pheatmap(gsva.res,show_colnames =T,show_rownames =F,
scale ="row",angle_col ="45",cluster_cols=F,
color= colorRampPalette(c("navy","white","firebrick3"))(50))


sort(colnames(gsva.df)[-1])

group_list <-data.frame(sample =sort(colnames(gsva.df)[-1]),group =c(rep("EV",4),rep("WT",4)))
group_list

design <-model.matrix(~0 +factor(group_list$group))
colnames(design)<-levels(factor(group_list$group))
rownames(design)<-colnames(gsva.res)
# rownames(design)<-c('aneuploid','aneuploid')
design


#构建差异比较矩阵
contrast.matrix <-makeContrasts(EV-WT,levels = design)
#差异分析，case vs.con
fit <-lmFit(gsva.res,design)
fit2 <-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
x <-topTable(fit2,coef =1,n =Inf,adjust.method ="BH",sort.by ="P")
head(x)

#
# 把通路的1imma分析结果保存到文件
# write.csv(x,"gsva_limma.csv",quote F)
# 输出t值，用做绘图的输入数据
# pathway <-str_replace(row.names(x),"HALLMARK_","")
pathway <-row.names(x)
df <-data.frame(ID =pathway,score =x$t)
head(df)

# write.csv(df,"enrich_bar.csv",quote F,row.names F)

# 按照scoref的值分组
cutoff <-0
df$group <-cut(df$score,breaks =c(-Inf,cutoff,Inf),labels =c(1,2))
head(df)

# 按照score排序
sortdf <-df[order(df$score),]
sortdf$ID <-factor(sortdf$ID,levels =sortdf$ID)
head(sortdf)
#

# ggplot(sortdf,aes(ID,score,fill =group))+geom_bar(stat ='identity')+
# coord_flip()+
# scale_fill_manual(values =c('palegreen3','dodgerblue4'),guide= FALSE)+
# #画2条虚线
# geom_hline(yintercept= c(-1,1),
# color="white",
# linetype=2,#画虚线
# size=0.3)+#线的粗细
# #写1abel
# geom_text(data= subset(df,score < 0),
# aes(x=ID,y=0.1,label=ID,color =group),
# size=3,#字的大小
# hjust="outward")+#字的对齐方式
# geom_text(data =subset(df,score >0),
# aes(x=ID,y=-0.1,label=paste0("",ID),color=group),#bar跟坐标轴间留出间隙
# size =3,hjust ="inward")+
# scale_colour_manual(values =c("black","black"),guide= FALSE)+
# xlab("")+ylab("t value of GSVA score")+
# theme_bw()+#去除背景色
# theme(panel.grid=element_blank())+#去除网格线
# theme(panel.border=element_rect(size=0.6))+#边框粗细
# theme(axis.line.y =element_blank(),axis.ticks.y =element_blank(),axis.text.y =element_blank())
# ggsave("gsva.pdf",width 6,height 8)

########################## 前后50
h<-head(sortdf,25)
t<-tail(sortdf,25)
ht<-rbind(h,t)
ggplot(ht,aes(ID,score,fill =group))+geom_bar(stat ='identity')+
coord_flip()+
scale_fill_manual(values =c('palegreen3','dodgerblue4'),guide= FALSE)+
#画2条虚线
geom_hline(yintercept= c(-1,1),
color="white",
linetype=2,#画虚线
size=0.3)+#线的粗细
#写1abel
geom_text(data= subset(ht,score > 0),
aes(x=ID,y=-0.01,label=ID,color =group),
size=3,#字的大小
hjust="outward")+#字的对齐方式
geom_text(data =subset(ht,score <0),
aes(x=ID,y=0.01,label=paste0("",ID),color=group),#bar跟坐标轴间留出间隙
size =3,hjust ="inward")+
scale_colour_manual(values =c("black","black"),guide= FALSE)+
xlab("")+ylab("t value of GSVA score")+
theme_bw()+#去除背景色
theme(panel.grid=element_blank())+#去除网格线
theme(panel.border=element_rect(size=0.6))+#边框粗细
theme(axis.line.y =element_blank(),axis.ticks.y =element_blank(),axis.text.y =element_blank())