library(Seurat)
library(ggplot2)
setwd("/Projects/deng/Rectal")
Rectal=readRDS("/Projects/deng/Rectal/Rectal.rds")
Myeloid=subset(Rectal,cellType %in% c("Macrophages","Neutrophils","Mast cells"))
Myeloid <- RunPCA(Myeloid, features = VariableFeatures(object = Myeloid),npcs = 100)
pdf("Visualization/MyeloidElbowPlot.pdf")
ElbowPlot(Myeloid,ndims = 50)
dev.off()
Myeloid<- RunTSNE(Myeloid, reduction = "pca", dims = 1:20)
Myeloid<- FindNeighbors(Myeloid, reduction = "pca", dims = 1:20)
Myeloid<- FindClusters(Myeloid, resolution = 0.5)


Myeloid.markers <- FindAllMarkers(Myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Myeloid.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top3
library(RColorBrewer)
display.brewer.all()
mapal <- colorRampPalette(RColorBrewer::brewer.pal(9,"PuRd"))(256)
DoHeatmap(Myeloid, features = top3$gene) + NoLegend()+  scale_fill_gradientn(colours = mapal)
pdf("Myeloid/MyeloidMarker.pdf",height=7,width=9)
DoHeatmap(Myeloid, features = top3$gene) + NoLegend()+ scale_fill_viridis(option="heat")
dev.off()
write.table(Myeloid.markers,file="MyeloidMarker.txt",sep="\t",quote=F)

#Subytpe marker visualization
pdf("MyeloidSubTypeMarker.pdf",width=9)
FeaturePlot(Myeloid, c("GPNMB","LGMN","C1QA","C1QB","C1QC","VCAN","FCN1","EREG","APOC1","HLA-DPB1","CD74","S100A9","S100A8"))&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


#DEG analysis
setwd("/Projects/deng/Rectal/Myeloid/DEG")
Myeloid$Category=paste0(Myeloid$seurat_clusters,"_",Myeloid$group,sep="")
Idents(Myeloid)=Myeloid$Category
MyeloidScDEG <- FindMarkers(Myeloid, ident.1 = "0_SC", ident.2 = "0_Ctrl", verbose = FALSE,min.pct = 0.25)
MyeloidScDEG=MyeloidScDEG[MyeloidScDEG$ p_val_adj<0.01,]
write.table(MyeloidScDEG,file="0_Myeloids_SC.txt",quote=F,sep="\t")
MyeloidPcDEG <- FindMarkers(Myeloid, ident.1 = "0_PC", ident.2 = "0_Ctrl", verbose = FALSE,min.pct = 0.25)
MyeloidPcDEG=MyeloidPcDEG[MyeloidPcDEG$ p_val_adj<0.01,]
write.table(MyeloidPcDEG,file="0_Myeloids_PC.txt",quote=F,sep="\t")


saveRDS(Myeloid,file="/Projects/deng/Rectal/Myeloid/Myeloid.rds")



MyeloidScDEG <- FindMarkers(Myeloid, ident.1 = "0_SC", ident.2 = "0_Ctrl", verbose = FALSE,min.pct = 0.25)
MyeloidPcDEG <- FindMarkers(Myeloid, ident.1 = "0_PC", ident.2 = "0_Ctrl", verbose = FALSE,min.pct = 0.25)
library(ggrepel)
hs_data=MyeloidScDEG
hs_data=hs_data[!rownames(hs_data)%in%c("XIST"),]
hs_data$threshold = as.factor(ifelse(hs_data$p_val_adj < 0.01 & abs(hs_data$avg_log2FC)>0, ifelse(hs_data$avg_log2FC > 0,"Up",'Down'),'Others'))
table(hs_data$threshold)
hs_data$ID=rownames(hs_data)
t=ggplot(data = hs_data, aes(x = avg_log2FC, y = -log10(p_val_adj), colour=threshold, label =ID )) +
  geom_point(alpha=0.4, size=3) +
  theme_bw() + 
  ylim(0,20)+xlim(-2,2)+
  scale_color_manual(values=c("blue","grey","red")) +
  geom_vline(xintercept=c(-log2(1.5),log2(1.5)),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2FC",y="-log10 (p-value)",title="Differential Expressed Genes") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
    geom_text_repel(
    data = subset(hs_data, abs(hs_data$avg_log2FC)>0.58 & hs_data$p_val_adj<0.01),
    aes(label = ID),
    size = 3,
    max.overlaps=10,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
pdf("Myeloid/MyeloidScDEG.pdf",width=5,height=4)
print(t)
dev.off()

geneCombine=read.table("D:/Application/Rectal/Myeloid/geneCombine.txt",sep="\t",header=T)
DEGNumber=data.frame(table(geneCombine$subType,geneCombine$Group))
colnames(DEGNumber)=c("SubType","Group","Number")
g=ggplot(DEGNumber, aes(SubType, Number, fill=Group)) +
  geom_bar(stat="identity",position=position_dodge()) +
  scale_y_continuous(expand=c(0,0))+theme_bw()+
  scale_fill_manual(values=c("MediumPurple","Magenta"))+
  theme(axis.text.x = element_text(angle = 0))+
  theme(axis.title.x = element_blank())
pdf("D:/Application/Rectal/Myeloid/geneCombineDEGDis.pdf",width=4,height=4)
print(g)
dev.off()


Myeloid=subset(Myeloid,idents=c(0:4))
tmp=paste0(Myeloid$orig.ident,"_",Myeloid$seurat_clusters,sep="")
tmp=as.matrix(table(tmp))
TmpInfo=as.matrix(do.call(rbind, strsplit(as.character(rownames(tmp)),'_')))
result=data.frame("Sample"=TmpInfo[,1],"CellClusters"=TmpInfo[,2],"Percent"=tmp[,1])
Cluster_order=factor(result$CellClusters,levels=c(0:4))
Sample_order=factor(result$Sample,levels=c("PN1","PN2","SCN1","SCN2","PC1","PC2","SC1","SC2"))
cbbPalette <- c("MediumTurquoise", " MediumAquamarine", "Turquoise", "LightSeaGreen", "MediumPurple", "BlueViolet", "Violet", "Magenta")
g=ggplot(result, aes(Cluster_order, Percent, fill=Sample_order)) +
  geom_bar(stat="identity",position="fill") +
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=cbbPalette)+
  theme(axis.text.x = element_text(angle = 0))+
  theme(axis.title.x = element_blank())
pdf("Visualization/MyeloidSampleDis.pdf",width=5,height=5)
print(g)
dev.off()



