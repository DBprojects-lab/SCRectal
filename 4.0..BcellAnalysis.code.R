library(Seurat)
library(ggplot2)
setwd("/Projects/deng/Rectal/Bcell")
Bcell=subset(Rectal,cellType %in% c("Plasma B cells","Naive B cells"))
pdf("BcellElbowPlot.pdf")
ElbowPlot(Bcell,ndims = 50)
dev.off()
Bcell<- RunTSNE(Bcell, reduction = "pca", dims = 1:50)
Bcell<- FindNeighbors(Bcell, reduction = "pca", dims = 1:50)
Bcell<- FindClusters(Bcell, resolution = 0.1)
pdf("BcellClusterLabel.pdf",height=4,width=5)
TSNEPlot(Bcell, label=TRUE)&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("BcellsMarker.pdf",width=9)
FeaturePlot(Bcell, c("MS4A1","MZB1","IGHG1","IGHA1"))&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
https://zh.wikipedia.org/wiki/%E6%8A%97%E4%BD%93
https://www.cell.com/trends/immunology/fulltext/S1471-4906(22)00003-5
Four B cell clusters identified:
Naïve B cells: MS4A1+IGHG1–;
Plasma cells: MZB1+IGHG1+;
Cycling plasma cells: MZB1+IGHG1+MKI67+;
Memory B cells: MS4A1+IGHG1+

pdf("BcellsGroup.pdf",height=4,width=6)
TSNEPlot(Bcell, group.by="group",cols = c("LightSkyBlue","MediumPurple","Magenta"))&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
new.cluster.ids <- c("IgA B cells","Naive B cells","IgG B cells")
names(new.cluster.ids) <- levels(Bcell)
Bcell <- RenameIdents(Bcell, new.cluster.ids)
Bcell$cellType=Idents(Bcell)
pdf("BcellType.pdf",height=4,width=5)
TSNEPlot(Bcell, group.by="cellType",label=TRUE)&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()




Bcell$Category=paste0(Bcell$cellType,"_",Bcell$group,sep="")
Idents(Bcell)=Bcell$Category
BcellScDEG <- FindMarkers(Bcell, ident.1 = "Naive B cells_SC", ident.2 = "Naive B cells_Ctrl", verbose = FALSE,min.pct = 0.25)
BcellScDEG=BcellScDEG[BcellScDEG$ p_val_adj<0.01,]
write.table(BcellScDEG,file="DEG/Naive_Bcells_SC.txt",quote=F,sep="\t")
BcellPcDEG <- FindMarkers(Bcell, ident.1 = "Naive B cells_PC", ident.2 = "Naive B cells_Ctrl", verbose = FALSE,min.pct = 0.25)
BcellPcDEG=BcellPcDEG[BcellPcDEG$ p_val_adj<0.01,]
write.table(BcellPcDEG,file="DEG/Naive_Bcells_PC.txt",quote=F,sep="\t")

geneCombine=read.table("D:/Application/Rectal/Bcell/geneCombine.txt",sep="\t",header=T)
DEGNumber=data.frame(table(geneCombine$subType,geneCombine$Group))
colnames(DEGNumber)=c("SubType","Group","Number")
g=ggplot(DEGNumber, aes(SubType, Number, fill=Group)) +
  geom_bar(stat="identity",position=position_dodge()) +
  scale_y_continuous(expand=c(0,0))+theme_bw()+
  scale_fill_manual(values=c("MediumPurple","Magenta"))+
  theme(axis.text.x = element_text(angle = 0))+
  theme(axis.title.x = element_blank())
pdf("D:/Application/Rectal/Bcell/geneCombineDEGDis.pdf",width=4,height=4)
print(g)
dev.off()


saveRDS(Bcell,file="/Projects/deng/Rectal/Bcell/BCell.rds")