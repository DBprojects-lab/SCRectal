setwd("/Projects/deng/Rectal/SCENIC")
library(SCENIC) 
packageVersion("SCENIC") #‘1.2.4’

hg38_dbs <- list('500bp'= 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather', 
                 '10kb' = 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather')
scenicOptions <- initializeScenic(org='hgnc',
                                datasetTitle='Rectal', 
                                dbDir="databases",
                                dbs=hg38_dbs,
                                nCores=20) 
EPC=readRDS("/Projects/deng/Rectal/EPC.rds")
table(EPC$seurat_clusters)
EPC=subset(EPC,idents=c(0:17))
exprMat  <-  as.matrix(EPC@assays$RNA@counts)
dim(exprMat)
exprMat[1:4,1:4] 
cellInfo <-  EPC@meta.data
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,minCountsPerGene=3*.01*ncol(exprMat),minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered, scenicOptions)

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 

scenicOptions <- initializeScenic(org="hgnc",dbDir="databases" , dbs=hg38_dbs, datasetTitle='Rectal', nCores=1)

scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat) #(Non-detected genes are shuffled at the end of the ranking. Keep it in mind when choosing the threshold for calculating the AUC)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 


hg38_dbs <- list('500bp'= 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather', 
                 '10kb' = 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather')
scenicOptions <- initializeScenic(org='hgnc',
                                datasetTitle='Rectal', 
                                dbDir="databases",
                                dbs=hg38_dbs,
                                nCores=1) 
EPC=readRDS("/Projects/deng/Rectal/EPC.rds")
EPC=subset(EPC,idents=c(0:17))
cellInfo <-  EPC@meta.data
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
table(EPC$seurat_clusters)

cellInfo <- data.frame(seuratCluster=EPC$seurat_clusters)
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster), function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType=regulonActivity_byCellType[,c(1:16)]
write.table(regulonActivity_byCellType,file="regulonActivity_byCellType.txt",sep="\t",quote=F)


library(viridis)
TFScore=read.table("D:/Application/Rectal/SCENIC/regulonActivity_byCellType.txt",header=T,row.names=1,sep="\t",check.names=F)
TFScore=TFScore[rowMax(TFScore)>0.01,]
pheatmap(TFScore,scale="row",clustering_method="ward.D2",color = viridis(10),show_rownames = F)
pdf("D:/Application/Rectal/SCENIC/heatmapDetail.pdf",height=15)
pheatmap(TFScore,scale="row",clustering_method="ward.D2",color = colorRampPalette(c( "white","white","white","pink","red"))(50),show_rownames = T)
dev.off()

##Identify the targets fo Myc
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF %in% c("MYC") & highConfAnnot==TRUE,]
tarGenes=tableSubset$gene
exprMat  <-  AverageExpression(EPC)
exprMatTmp  <-  exprMat$RNA[tarGenes,]
write.table(exprMatTmp,file="MYCtarget.txt",sep="\t",quote=F)
exprMatTmp=read.table("D:/Application/Rectal/SCENIC/MYCtarget.txt",header=T,row.names=1,check.names=F)
pdf("MYCTargetExpr.pdf")
pheatmap(exprMatTmp ,scale="row",clustering_method="ward.D2",color = colorRampPalette(c( "navy","white","red"))(50),show_rownames = F)
dev.off()

#check the expression of Myc
pdf("MYCExprAll.pdf",width=5)
VlnPlot(EPC,features=c("MYC","ELF1","ELF3","HNF4A","PPARG"),pt.size=0,ncol=1)&
theme(plot.title=element_text(size=10),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y=element_text(size=8))
dev.off()


