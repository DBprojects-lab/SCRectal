library(Seurat)
library(ggplot2)
setwd("/Projects/deng/Rectal")
N.data <- Read10X(data.dir = "/Projects/deng/Rectal/data/N/filtered_feature_bc_matrix")
N <- CreateSeuratObject(counts = N.data, project = "SCN1")
N$group="Ctrl"
N1.data <- Read10X(data.dir = "/Projects/deng/Rectal/data/N1/filtered_feature_bc_matrix")
N1 <- CreateSeuratObject(counts = N1.data, project = "SCN2")
N1$group="Ctrl"
T.data <- Read10X(data.dir = "/Projects/deng/Rectal/data/T/filtered_feature_bc_matrix")
T <- CreateSeuratObject(counts = T.data, project = "SC1")
T$group="SC"
T1.data <- Read10X(data.dir = "/Projects/deng/Rectal/data/T1/filtered_feature_bc_matrix")
T1 <- CreateSeuratObject(counts = T1.data, project = "SC2")
T1$group="SC"
R18072055.data<- Read10X(data.dir = "/Projects/deng/Rectal/data/R18072055")
R18072055 <- CreateSeuratObject(counts = R18072055.data, project = "PC1")
R18072055$group="PC"
R18072056.data<- Read10X(data.dir = "/Projects/deng/Rectal/data/R18072056/")
R18072056 <- CreateSeuratObject(counts = R18072056.data, project = "PN1")
R18072056$group="Ctrl"
R18072059.data<- Read10X(data.dir = "/Projects/deng/Rectal/data/R18072059")
R18072059 <- CreateSeuratObject(counts = R18072059.data, project = "PC2")
R18072059$group="PC"
R18072060.data<- Read10X(data.dir = "/Projects/deng/Rectal/data/R18072060")
R18072060 <- CreateSeuratObject(counts = R18072060.data, project = "PN2")
R18072060$group="Ctrl"

RectalB1 <- merge(N, y = c(N1,T,T1), add.cell.ids = c("Normal", "Normal1","Tumor", "Tumor1"), project = "FirstBatch")
RectalB2 <- merge(R18072055, y = c(R18072056,R18072059,R18072060), add.cell.ids = c("R55", "R56","R59", "R60"), project = "SecondBatch")
RectalAll <- merge(RectalB1,RectalB2)
RectalAll[["percent.mt"]] <- PercentageFeatureSet(RectalAll, pattern = "^MT-")
RectalAll[["percent.ribo"]] <- PercentageFeatureSet(RectalAll, pattern = "^RP[SL]")

#check the QC
VlnPlot(RectalAll, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4, pt.size=0)

dim(RectalAll) #36601 32156
#filtered by gene epxression number
#only consider cells with at least 200 detected genes and genes need to be expressed in at least 3 cells. the cells expressed more than 7500 genes were also removed, indicate doublets. 
selected_c <- WhichCells(RectalAll, expression = nFeature_RNA > 200 & nFeature_RNA < 7500)
selected_f <- rownames(RectalAll)[Matrix::rowSums(RectalAll) > 3]
data.filt <- subset(RectalAll, features = selected_f, cells = selected_c)
dim(data.filt) #24985 28986

#filtered by mitrochondril percentage and ribosome percentage
selected_mito <- WhichCells(data.filt, expression = percent.mt < 25)
selected_ribo <- WhichCells(data.filt, expression = percent.ribo > 3)
# and subset the object to only keep those cells
data.filt <- subset(data.filt, cells = selected_mito)
data.filt <- subset(data.filt, cells = selected_ribo)
dim(data.filt) #24985 16334

# Filter Mitocondrial
data.filt <- data.filt[!grepl("^MT-", rownames(data.filt)), ]
dim(data.filt)  #24972 16334
data.filt <- data.filt[!grepl("^RP[SL]", rownames(data.filt)), ]
dim(data.filt)  #24872 16334


Rectal <- NormalizeData(data.filt, normalization.method = "LogNormalize", scale.factor = 10000)
Rectal <- FindVariableFeatures(Rectal, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Rectal)
Rectal <- ScaleData(Rectal, features = all.genes)
Rectal <- RunPCA(Rectal, features = VariableFeatures(object = Rectal),npcs = 100)
#DimHeatmap(Rectal , dims = 20:45, cells = 500, balanced = TRUE)
pdf("RectalElbowPlot.pdf")
ElbowPlot(Rectal,ndims = 100)
dev.off()
# t-SNE and Clustering
Rectal<- RunTSNE(Rectal, reduction = "pca", dims = 1:50)
Rectal<- FindNeighbors(Rectal, reduction = "pca", dims = 1:50)
Rectal<- FindClusters(Rectal, resolution = 0.5)
# Visualization
pdf("Visualization/clusterLabel.pdf",height=4,width=5)
TSNEPlot(Rectal, label=TRUE,label.size = 6)&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("Visualization/sampleLabel.pdf",height=4,width=6)
TSNEPlot(Rectal, group.by="orig.ident")&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("Visualization/GroupLabel.pdf",height=4,width=6)
TSNEPlot(Rectal, group.by="group",cols = c("LightSkyBlue","DarkViolet","Magenta"))&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

####check the expression of makers identified by previous study#################
FeaturePlot(Rectal,features=c("VIL1","KRT20", "CLDN7", "CDH1"))&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank()) #Epithelial cells
FeaturePlot(Rectal,features=c("CD38", "MZB1", "DERL3","MS4A1","CD79A","CD79B"))&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank()) #B cells
FeaturePlot(Rectal,features=c("ENG")) #endothelial cells
FeaturePlot(Rectal,features=c("COL3A1", "DCN","TSB2","KIT"))&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank()) #Epithelial cells 
FeaturePlot(Rectal,features=c("TRBC2", "CD3D", "CD3E", "CD3G"))&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank()) #T cells
FeaturePlot(Rectal,features=c("KIT", "TPSB2","ENG"))&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank())  #mast cells & endothelial cells
FeaturePlot(Rectal,features=c("ITGAX", "CD68", "CD14", "CCL3"))&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank()) #myeloid cells
FeaturePlot(Rectal,c("HLA-DQB1","HLA-DQA1","MS4A1","HLA-DRA"))&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank()) #Dendritic cell

DCMarker=c("CD8A","ITGAE","CD4","ITGAM","BANK1");
pdf("Visualization/DCMarker.pdf",height=8,width=10)
FeaturePlot(Rectal,DCMarker,pt.size=0.3)&NoLegend()&
theme(plot.title=element_text(size=10),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


####check the expression of makers specific for each cluster#################
marker=c("CD3D","CLDN7","SPINK1","CD79A","HLA-DRA","BANK1","COL3A1","ENG","G0S2","CD14","TYMS","KIT");
pdf("Visualization/Marker.pdf",height=8,width=10)
FeaturePlot(Rectal,marker,pt.size=0.3)&NoLegend()&
theme(plot.title=element_text(size=10),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


####quantify the relatinship between clusters to help to identify the cell types#################
RectalAveExpr <- AverageExpression(Rectal)[["RNA"]]
pdf("Visualization/RectalCellTypeRelation.pdf",height=5)
plot(hclust(as.dist(1-cor(RectalAveExpr)),method="ward.D2"),hang=-1)
dev.off()
Idents(Rectal)=Rectal$seurat_clusters
RectalAveExpr <- AverageExpression(Rectal)[["RNA"]]
pdf("Visualization/RectalClusterRelation.pdf",height=5)
plot(hclust(as.dist(1-cor(RectalAveExpr))),hang=-1)
dev.off()

####Identification the noval markers and validation the previous markers for each cluster#################
for (i in 0:20){
DEG <- FindMarkers(Rectal, ident.1 = i, verbose = FALSE,min.pct = 0.25)
DEG=DEG[DEG$p_val<0.01,]
write.table(DEG,file=paste("/Projects/deng/Rectal/Markers/Cluster",i,"Deg.txt",sep=""),sep="\t",quote=F)
}

#Remove the clusters with extremly lower cells
Rectal=subset(Rectal,idents=c(0:20))

new.cluster.ids <- c("T cells","T cells","T cells","Epithelial cells","Plasma B cells","Epithelial cells","Epithelial cells","Macrophages","T cells","Epithelial cells","Epithelial cells","Epithelial cells","Naive B cells","Fibroblasts","Plasma B cells","Endothelial cells","Neutrophils","Progenitor cell","Mast cells","Epithelial cells","Fibroblasts")
names(new.cluster.ids) <- levels(Rectal)
Rectal <- RenameIdents(Rectal, new.cluster.ids)
Rectal$cellType=Idents(Rectal)

pdf("Visualization/CellTypeName.pdf",height=4,width=7)
TSNEPlot(Rectal)&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

marker=c("CD3D","CLDN7","SPINK1","MZB1","CD68","BANK1","COL3A1","ENG","G0S2","CLSPN","MKI67","KIT");
pdf("Visualization/MarkerVlnPlot.pdf",height=10,width=7)
VlnPlot(Rectal,features=marker,pt.size=0,ncol=1)&
theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=8))
dev.off()

pdf("Visualization/MarkerFeaturePlot.pdf",height=8,width=10)
FeaturePlot(Rectal,marker,pt.size=0.3)&NoLegend()&
theme(plot.title=element_text(size=10),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

Rectal$cellType=Idents(Rectal)
t=as.matrix(table(Rectal$cellType))
tmp=paste0(Rectal$orig.ident,"_",Rectal$cellType,sep="")
tmp=as.matrix(table(tmp))
TmpInfo=as.matrix(do.call(rbind, strsplit(as.character(rownames(tmp)),'_')))
result=data.frame("Sample"=TmpInfo[,1],"CellType"=TmpInfo[,2],"Percent"=tmp[,1])


Sample_order=factor(result$Sample,levels=c("PN1","PN2","SCN1","SCN2","PC1","PC2","SC1","SC2"))
CellType_order=factor(result$CellType,levels=c("Neutrophils","T cells","Epithelial cells","Plasma B cells","Macrophages","B cells","Premature B cell","Fibroblasts","Endothelial cells","Progenitor cell","Mast cells"))
g=ggplot(result, aes(Sample_order, Percent, fill=CellType_order)) +
  geom_bar(stat="identity",position="fill") +
  scale_y_continuous(expand=c(0,0))+
  theme(axis.text.x = element_text(angle = 0))+
  theme(axis.title.x = element_blank())
pdf("Visualization/cellTypeDistribution.pdf",height=6,width=10)
print(g)
dev.off()

saveRDS(Rectal,"/Projects/deng/Rectal/Rectal.rds")


