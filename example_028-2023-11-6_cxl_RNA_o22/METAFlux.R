library(METAFlux)
data("human_gem")
data("human_blood")# for human derived samples
bulk_test_example<-read.table('/scratch/wus/028-2023-11-6_cxl_RNA_o22/8_DEseq2/condition_A_mut_vs_E_mut.txt',row.names='Row.names')

head(bulk_test_example[,-c(1:7)])

bulk_test_example<-bulk_test_example[,-c(1:7)]

scores<-calculate_reaction_score(bulk_test_example)
flux<-compute_flux(mras=scores,medium = human_blood)#if data are human derived samples
cbrt <- function(x) {
    sign(x) * abs(x)^(1/3)
}

flux=cbrt(flux)
library(ggplot2)
#compute pathway level activity for all samples
pathway<-unique(unlist(human_gem$SUBSYSTEM))
pathway_score<-list()
for (i in pathway){
path=i
activity_score<-c()
for (d in 1:ncol(flux)){
activity_score[d]<-mean(abs(flux[which(unlist(human_gem$SUBSYSTEM)==i),d]))
} 
pathway_score[[i]]<-activity_score
}

all_pathway_score<-as.data.frame(do.call(rbind,pathway_score))

colnames(all_pathway_score)<-colnames(flux)

#heatmap 
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
pheatmap::pheatmap(all_pathway_score,cluster_cols = F,color = rev(mapal),scale = "row")


####
# all_pathway_score[,!grepl('E',colnames(all_pathway_score))]
all_pathway_score[,!grepl('E',colnames(all_pathway_score))]
pheatmap::pheatmap(all_pathway_score[,!grepl('E',colnames(all_pathway_score))],cluster_cols = F,color = rev(mapal),scale = "row")


# wilcoxon_test
# i=1
# wilcox.test(as.vector(unlist(all_pathway_score[i,1:12])), as.vector(unlist(all_pathway_score[i,13:16])), paired = FALSE)

all_pathway_score<-all_pathway_score[,!grepl('E',colnames(all_pathway_score))]

p_value=c()
for (i in 1:nrow(all_pathway_score)){
	#a<-wilcox.test(as.vector(unlist(all_pathway_score[i,1:12])), as.vector(unlist(all_pathway_score[i,13:16])), paired = FALSE)
	a<-wilcox.test(as.vector(unlist(all_pathway_score[i,1:12])), as.vector(unlist(all_pathway_score[i,13:32])), paired = FALSE)

p_value[i]<-a$p.value
}
all_pathway_score$p<-p_value
all_pathway_score_p<-all_pathway_score[all_pathway_score$p<0.05,]

all_pathway_score_p[,1:32]
pheatmap::pheatmap(all_pathway_score_p[,1:32],cluster_cols = F,color = rev(mapal),scale = "row")
