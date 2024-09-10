# conda activate DESeq2
library(WGCNA)
library(dplyr)
library(stringr )
df<-read.table('/home/wus/2023_3_19-009-cxl_RNA/7-25/condition_KR234_vs_WT.txt',sep='')


## 按索引删除
df <- select(df,c(-2,-3,-4,-5,-6,-7))

femData<-df
dim(femData)
names(femData)
#行为基因，列为不同样本的基因表达量或其他信息
#提取出表达量的数据 ，删去不需要的数据重新生成矩阵

datExpr0 = as.data.frame(t(femData[, -c(1:1)]))  #提取加转置
names(datExpr0) = femData$Row.names #基因名字
# rownames(datExpr0) = names(femData)[-c(2:36)]  #样品名字
datExpr0 #就是一个以每行为样本，一列为一个基因的数据框


sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(12,9) #视图
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


# 删除离群样本
abline(h = 4e+05, col = "red") #划定需要剪切的枝长
clust = cutreeStatic(sampleTree, cutHeight = 4e+05, minSize = 10)
# 这时候会从高度为15这里横切，把离群样本分开
table(clust)   
keepSamples = (clust==1)  #保留非离群(clust==1)的样本
datExpr = datExpr0[keepSamples, ]  #去除离群值后的数据
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()  #开启多线程
allowWGCNAThreads()
ALLOW_WGCNA_THREADS=24
memory.limit(size = 200000)


dataExpr<-datExpr
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, networkType="unsigned")

# *- 出图 -*
# 设置绘制界面分为2个子图，1行2列布局
par(mfrow = c(1,2))
cex1 = 0.9

# -左侧子图：基因连通性和其数量的线性相关系数图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))

# 对应位置填充β取值
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")

# 添加筛选基线，R-square=0.85
abline(h=0.85,col="red")

# - 右侧子图：平均连通性图
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))

# 对应位置填充β取值
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,cex=cex1, col="red")



net = blockwiseModules(datExpr, power = 14,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)


sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

################# 表型
traitData = read.csv("/home/wus/2023_3_19-009-cxl_RNA/11-3_GWCNA-RRA/datTraits.csv",row.names=1);
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE) #用颜色代表关联度
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
datTraits=traitData


nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# 重新计算带有颜色标签的模块
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# 通过相关值对每个关联进行颜色编码
sizeGrWindow(10,6)
# 展示模块与表型数据的相关系数和 P值
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# 用热图的形式展示相关系数
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#colors = greenWhiteRed(50)不适用于红绿色盲患者，建议用 blueWhiteRed代替.
#该分析确定了几个重要的模块-特征关联。我们将体重作为感兴趣的特征来研究。

names(datExpr)[moduleColors=="pink"]
write.csv(names(datExpr)[moduleColors=="pink"],'/home/wus/2023_3_19-009-cxl_RNA/11-3_GWCNA-RRA/pink.csv')

names(datExpr)[moduleColors=="yellow"]
write.csv(names(datExpr)[moduleColors=="yellow"],'/home/wus/2023_3_19-009-cxl_RNA/11-3_GWCNA-RRA/yellow.csv')

######### 表型 可以是肿瘤大小
# datTraits<-datExpr
# datTraits$sample=rownames(datExpr)
# datTraits<-datTraits[,'sample',drop=FALSE]
# write.csv(datTraits,'/home/wus/2023_3_19-009-cxl_RNA/11-3_GWCNA-RRA/datTraits.csv')

# # sampleTree2 = hclust(dist(datExpr), method = "average")
# # traitColors = numbers2colors(datTraits, signed = FALSE) #用颜色代表关联度
# # plotDendroAndColors(sampleTree2, traitColors,
# #                     groupLabels = names(datTraits),
# #                     main = "Sample dendrogram and trait heatmap")


# nGenes = ncol(datExpr);
# nSamples = nrow(datExpr);
# # 重新计算带有颜色标签的模块
# MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
# MEs = orderMEs(MEs0)
# moduleTraitCor = cor(MEs, datTraits, use = "p");
# moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# # 通过相关值对每个关联进行颜色编码
# sizeGrWindow(10,6)
# # 展示模块与表型数据的相关系数和 P值
# textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
#                    signif(moduleTraitPvalue, 1), ")", sep = "");
# dim(textMatrix) = dim(moduleTraitCor)
# par(mar = c(6, 8.5, 3, 3));
# # 用热图的形式展示相关系数
# labeledHeatmap(Matrix = moduleTraitCor,
#                xLabels = names(datTraits),
#                yLabels = names(MEs),
#                ySymbols = names(MEs),
#                colorLabels = FALSE,
#                colors = greenWhiteRed(50),
#                textMatrix = textMatrix,
#                setStdMargins = FALSE,
#                cex.text = 0.5,
#                zlim = c(-1,1),
#                main = paste("Module-trait relationships"))
# #colors = greenWhiteRed(50)不适用于红绿色盲患者，建议用 blueWhiteRed代替.
# #该分析确定了几个重要的模块-特征关联。我们将体重作为感兴趣的特征来研究。


# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
# MEs = net$MEs

# ### 不需要重新计算，改下列名字就好
# ### 官方教程是重新计算的，起始可以不用这么麻烦
# MEs_col = MEs
# colnames(MEs_col) = paste0("ME", labels2colors(
#   as.numeric(str_replace_all(colnames(MEs),"ME",""))))
# MEs_col = orderMEs(MEs_col)

# # 根据基因间表达量进行聚类所得到的各模块间的相关性图
# # marDendro/marHeatmap 设置下、左、上、右的边距
# plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
#                       marDendro = c(3,3,2,4),
#                       marHeatmap = c(3,4,2,2), plotDendrograms = T, 
#                       xLabelsAngle = 90)


# # 如果采用分步计算，或设置的blocksize>=总基因数，直接load计算好的TOM结果
# # 否则需要再计算一遍，比较耗费时间
# type = "unsigned"

# # 相关性计算
# # 官方推荐 biweight mid-correlation & bicor
# # corType: pearson or bicor
# # 为与原文档一致，故未修改
# corType = "pearson"

# # TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
# load(net$TOMFiles[1], verbose=T)

# ## Loading objects:
# ##   TOM

# TOM <- as.matrix(TOM)

# dissTOM = 1-TOM
# # Transform dissTOM with a power to make moderately strong 
# # connections more visible in the heatmap
# plotTOM = dissTOM^7
# # Set diagonal to NA for a nicer plot
# diag(plotTOM) = NA
# # Call the plot function

# # 这一部分特别耗时，行列同时做层级聚类
# TOMplot(plotTOM, net$dendrograms, moduleColors, 
#         main = "Network heatmap plot, all genes")



# probes = colnames(datExpr)
# dimnames(TOM) <- list(probes, probes)

# # Export the network into edge and node list files Cytoscape can read
# # threshold 默认为0.5, 可以根据自己的需要调整，也可以都导出后在
# # cytoscape中再调整
# cyt = exportNetworkToCytoscape(TOM,
#              edgeFile = paste(exprMat, ".edges.txt", sep=""),
#              nodeFile = paste(exprMat, ".nodes.txt", sep=""),
#              weighted = TRUE, threshold = 0,
#              nodeNames = probes, nodeAttr = moduleColors)




# weight = as.data.frame(datTraits$KR);
# names(weight) = "KR";
# modNames = substring(names(MEs), 3)
# geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
# MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
# names(geneModuleMembership) = paste("MM", modNames, sep="");
# names(MMPvalue) = paste("p.MM", modNames, sep="");
# geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));#和体重性状的关联
# GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
# names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
# names(GSPvalue) = paste("p.GS.", names(weight), sep="");


# 运行以下代码可视化GS和MM
# module = "pink"
# column = match(module, modNames);
# moduleGenes = moduleColors==module;
# sizeGrWindow(7, 7);
# par(mfrow = c(1,1));
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                    abs(geneTraitSignificance[moduleGenes, 1]),
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = "Gene significance for body weight",
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



# #########################
