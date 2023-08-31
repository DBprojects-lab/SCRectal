library(Seurat)
library(clustree)
library(ggplot2)
setwd("/Projects/deng/Rectal")
Idents(Rectal)=Rectal$cellType
Tcell=subset(Rectal,idents="T cells")
Tcell <- RunPCA(Tcell, features = VariableFeatures(object = Tcell),npcs = 100)
pdf("Visualization/TcellElbowPlot.pdf")
ElbowPlot(Tcell,ndims = 50)
dev.off()
Tcell<- RunTSNE(Tcell, reduction = "pca", dims = 1:20)
Tcell<- FindNeighbors(Tcell, reduction = "pca", dims = 1:20)

#check the resolutions
TcellTmp <- FindClusters(Tcell, resolution = c(seq(0.1,1.5,0.1)))
clus.tree.out <- clustree(TcellTmp)
pdf(file = "Visualization/TcellTree.pdf", width = 15, height = 20)
print(clus.tree.out)
dev.off()
pdf(file = "Visualization/TcellClusterLabel1.0.pdf",height=4,width=5)
TSNEPlot(TcellTmp, group="RNA_snn_res.1",label=TRUE)&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

Tcell<- FindClusters(Tcell, resolution = 1)

pdf("Visualization/TcellClusterLabel.pdf",height=4,width=5)
TSNEPlot(Tcell, label=TRUE)&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("Visualization/TcellSampleLabel.pdf",height=4,width=6)
TSNEPlot(Tcell, group.by="orig.ident")&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("Visualization/TcellGroupLabel.pdf",height=4,width=6)
TSNEPlot(Tcell, group.by="group",cols = c("LightSkyBlue","MediumPurple","Magenta"))&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

###Check the expression of the previous validated markers
pdf("Visualization/TcellCD4-CD8.pdf",height=4,width=9)
FeaturePlot(Tcell,c("CD4","CD8A"))&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),,panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid")) #Dendritic cell
dev.off()
pdf("TcellCD4-CD8VlnPlot.pdf",height=4,width=8)
VlnPlot(Tcell,features=c("CD4","CD8A"),pt.size=0.1,ncol=1)&NoLegend()&
theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("TcellTregVlnPlot.pdf",height=4,width=8)
VlnPlot(Tcell,features=c("HAVCR2","CD274","TIGIT"),pt.size=0.1,ncol=1)&NoLegend()&
theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("TcellTregGene.pdf")
FeaturePlot(Tcell,c("HAVCR2","CD274","TIGIT","CTLA4"))&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),,panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid")) #Dendritic cell
dev.off()
pdf("TcellGZMK.pdf")
FeaturePlot(Tcell,c("GZMK","KLRG1"))&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),,panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid")) #Dendritic cell
dev.off()
pdf("Visualization/TcellActivation.pdf")
FeaturePlot(Rectal,c("CD28","CD80","CD86","CD40","CD40LG","ICOS","ICOSLG"))&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),,panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid")) #Dendritic cell
dev.off()
pdf("Visualization/TcellRelatedMarker.pdf")
FeaturePlot(Tcell,c("LAG3","PDCD1","TIGIT","IFNG","ZFP36","SG15","IFIT1","OAS1","MKI67","AREG"))&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),,panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid")) #Dendritic cell
dev.off()


###obtian the cell cluster specific markers
Tcell.markers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Tcell.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) -> top3
library(RColorBrewer)
display.brewer.all()
mapal <- colorRampPalette(RColorBrewer::brewer.pal(9,"PuRd"))(256)
DoHeatmap(Tcell, features = top3$gene) + NoLegend()+  scale_fill_gradientn(colours = mapal)
pdf("Visualization/TcellMarker.pdf",height=7,width=9)
DoHeatmap(Tcell, features = top3$gene) + NoLegend()+ scale_fill_viridis(option="heat")
dev.off()
write.table(Tcell.markers,file="TcellMarker.txt",sep="\t",quote=F)


###obtian the cell ratio for each cluster
tmp=paste0(Tcell$orig.ident,"_",Tcell$seurat_clusters,sep="")
tmp=as.matrix(table(tmp))
TmpInfo=as.matrix(do.call(rbind, strsplit(as.character(rownames(tmp)),'_')))
result=data.frame("Sample"=TmpInfo[,1],"CellClusters"=TmpInfo[,2],"Percent"=tmp[,1])
Cluster_order=factor(result$CellClusters,levels=c(0:15))
Sample_order=factor(result$Sample,levels=c("PN1","PN2","SCN1","SCN2","PC1","PC2","SC1","SC2"))
cbbPalette <- c("MediumTurquoise", " MediumAquamarine", "Turquoise", "LightSeaGreen", "MediumPurple", "BlueViolet", "Violet", "Magenta")
g=ggplot(result, aes(Cluster_order, Percent, fill=Sample_order)) +
  geom_bar(stat="identity",position="fill") +
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=cbbPalette)+
  theme(axis.text.x = element_text(angle = 0))+
  theme(axis.title.x = element_blank())
pdf("Visualization/TcellSampleDis.pdf",width=14,height=5)
print(g)
dev.off()



####Perfrom the DEG analysis between SCs/PCs and control
Tcell=readRDS("/Projects/deng/Rectal/Tcell.rds")
Tcell$Category=paste0(Tcell$seurat_clusters,"_",Tcell$group,sep="")
Idents(Tcell)=Tcell$seurat_clusters
Tcell=subset(Tcell,orig.ident %in% c("SC1","SCN1"))
table(Tcell$orig.ident)

setwd("/Projects/deng/Rectal/PairDEG/Tcell")
  for(i in unique(Tcell$seurat_clusters)){
   CaseGroup=paste(i,"_SC",sep="");
   CtlGroup=paste(i,"_Ctrl",sep="");
   if(is.na(table(Tcell$Category)[CaseGroup])) next;
   if(is.na(table(Tcell$Category)[CtlGroup])) next;
   if(table(Tcell$Category)[CaseGroup]<10) next;
   if(table(Tcell$Category)[CtlGroup]<10) next;
   DEG <- FindMarkers(Tcell, ident.1 = CaseGroup, ident.2 = CtlGroup, verbose = FALSE,min.pct = 0.25)
   DEG=DEG[DEG$ p_val_adj<0.01,]
   write.table(DEG,file=paste("SC1/",i,"Deg.txt",sep=""),sep="\t",quote=F)
}

Tcell=readRDS("/Projects/deng/Rectal/Tcell.rds")
Tcell$Category=paste0(Tcell$seurat_clusters,"_",Tcell$group,sep="")
Idents(Tcell)=Tcell$Category
table(Tcell$orig.ident)
setwd("/Projects/deng/Rectal/PairDEG/Tcell")

for(i in unique(Tcell$seurat_clusters)){
   CaseGroup=paste(i,"_SC",sep="");
   CtlGroup=paste(i,"_Ctrl",sep="");
   if(is.na(table(Tcell$Category)[CaseGroup])) next;
   if(is.na(table(Tcell$Category)[CtlGroup])) next;
   if(table(Tcell$Category)[CaseGroup]<10) next;
   if(table(Tcell$Category)[CtlGroup]<10) next;
   DEG <- FindMarkers(Tcell, ident.1 = CaseGroup, ident.2 = CtlGroup, verbose = FALSE,min.pct = 0.25)
   DEG=DEG[DEG$ p_val_adj<0.01,]
   write.table(DEG,file=paste("SCGroup/",i,"Deg.txt",sep=""),sep="\t",quote=F)
}

saveRDS(Tcell,"/Projects/deng/Rectal/Tcell.rds")
