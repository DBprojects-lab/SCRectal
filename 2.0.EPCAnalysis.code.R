library(Seurat)
library(clustree)
library(ggplot2)
setwd("/Projects/deng/Rectal")

#==========================EPC analysis==========================================
Rectal=readRDS("/Projects/deng/Rectal/Rectal.rds")

setwd("/Projects/deng/Rectal")
EPC=subset(Rectal,idents="Epithelial cells")
EPC <- RunPCA(EPC, features = VariableFeatures(object = EPC),npcs = 100)
pdf("Visualization/EPCElbowPlot.pdf")
ElbowPlot(EPC,ndims = 100)
dev.off()
EPC<- RunTSNE(EPC, reduction = "pca", dims = 1:50)
EPC<- FindNeighbors(EPC, reduction = "pca", dims = 1:50)
EPC<- FindClusters(EPC, resolution = 0.7)

#check the resolution###############
EPCTmp <- FindClusters(EPC, resolution = c(seq(0.1,1.5,0.1)))
clus.tree.out <- clustree(EPCTmp)
pdf(file = "Visualization/EPCTree.pdf", width = 15, height = 20)
print(clus.tree.out)
dev.off()
pdf(file = "Visualization/EPCClusterLabel0.7.pdf",height=4,width=5)
TSNEPlot(EPCTmp, group="RNA_snn_res.0.7",label=TRUE)&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


pdf("Visualization/EPCClusterLabel.pdf",height=4,width=5)
TSNEPlot(EPC, label=TRUE)&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("Visualization/EPCSampleLabel.pdf",height=4,width=6)
TSNEPlot(EPC, group.by="orig.ident")&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("Visualization/EPCGroupLabel.pdf",height=4,width=6)
TSNEPlot(EPC, group.by="group",cols = c("LightSkyBlue","MediumPurple","Magenta"))&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("GroupLabel.tiff",height=250,width=260)
TSNEPlot(EPC, group.by="group",cols = c("LightSkyBlue","MediumPurple","Magenta"))&NoLegend()&theme(axis.text.x = element_blank(),axis.text.y = element_blank(),plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

saveRDS(EPC,"/Projects/deng/Rectal/EPC.rds")

EPC=subset(EPC,idents=c(0:17))
EPCAveExpr <- AverageExpression(EPC)[["RNA"]]
pdf("Visualization/EPCClusterRelation.pdf",height=5)
plot(hclust(as.dist(1-cor(EPCAveExpr)),method="ward.D2"),hang=-1)
dev.off()

##detect the cell ratio from normal sample
Ctrl=subset(EPC,group=="Ctrl")
t=as.matrix(table(Ctrl$seurat_clusters))/as.matrix(table(EPC$seurat_clusters))
tmp=data.frame(cbind(paste0("C",c(0:19)),"Number"=t))
colnames(tmp)=c("Cluster","CtrlRatio")
tmp=tmp[c(1:18),]
ClusterOrder=factor(tmp$Cluster,levels=paste0("C",c(0:17)))
g=ggplot(tmp, aes(ClusterOrder, as.numeric(CtrlRatio), fill=Cluster)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 0))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
pdf("Visualization/EPCSampleCtrlDis.pdf",width=14,height=5)
print(g)
dev.off()

##detect the cell ratio
t=as.matrix(table(EPC$seurat_clusters))
tmp=paste0(EPC$orig.ident,"_",EPC$seurat_clusters,sep="")
tmp=as.matrix(table(tmp))
TmpInfo=as.matrix(do.call(rbind, strsplit(as.character(rownames(tmp)),'_')))
result=data.frame("Sample"=TmpInfo[,1],"CellClusters"=TmpInfo[,2],"Percent"=tmp[,1])
Cluster_order=factor(result$CellClusters,levels=c(17,2,8,6,11,16,15,1,7,12,3,14,9,13,4,5,0,10))
Sample_order=factor(result$Sample,levels=c("PN1","PN2","SCN1","SCN2","PC1","PC2","SC1","SC2"))
cbbPalette <- c("MediumTurquoise", " MediumAquamarine", "Turquoise", "LightSeaGreen", "MediumPurple", "BlueViolet", "Violet", "Magenta")
g=ggplot(result, aes(Cluster_order, Percent, fill=Sample_order)) +
  geom_bar(stat="identity",position="fill") +
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=cbbPalette)+
  theme(axis.text.x = element_text(angle = 0))+
  theme(axis.title.x = element_blank())
pdf("/Projects/deng/Rectal/Visualization/EPCSampleDis.pdf",width=14,height=5)
print(g)
dev.off()


setwd("/Projects/deng/Rectal/EPC")
EPC=readRDS("/Projects/deng/Rectal/EPC.rds")

####DEG analysis for the tumor cells
C0DEG=FindMarkers(EPC,ident.1=c(0),ident.2=c(17,2,8,6,11,16,15),pct.min=0.25)
C0DEG=C0DEG[C0DEG$p_val_adj<0.01,]
C0High=C0DEG[C0DEG$avg_log2FC>1,]
C0Low=C0DEG[C0DEG$avg_log2FC<-1,]
C0HighId=bitr(rownames(C0High),"SYMBOL","ENTREZID",OrgDb="org.Hs.eg.db")[,2]
C0LowId=bitr(rownames(C0Low),"SYMBOL","ENTREZID",OrgDb="org.Hs.eg.db")[,2]

C10DEG=FindMarkers(EPC,ident.1=c(10),ident.2=c(17,2,8,6,11,16,15),pct.min=0.25)
C10DEG=C10DEG[C10DEG$p_val_adj<0.01,]
C10High=C10DEG[C10DEG$avg_log2FC>1,]
C10Low=C10DEG[C10DEG$avg_log2FC<-1,]
C10HighId=bitr(rownames(C10High),"SYMBOL","ENTREZID",OrgDb="org.Hs.eg.db")[,2]
C10LowId=bitr(rownames(C10Low),"SYMBOL","ENTREZID",OrgDb="org.Hs.eg.db")[,2]

C13DEG=FindMarkers(EPC,ident.1=13,ident.2=c(17,2,8,6,11,16,15),pct.min=0.25)
C13DEG=C13DEG[C13DEG$p_val_adj<0.01,]
C13High=C13DEG[C13DEG$avg_log2FC>1,]
C13Low=C13DEG[C13DEG$avg_log2FC<-1,]
C13HighId=bitr(rownames(C13High),"SYMBOL","ENTREZID",OrgDb="org.Hs.eg.db")[,2]
C13LowId=bitr(rownames(C13Low),"SYMBOL","ENTREZID",OrgDb="org.Hs.eg.db")[,2]

C3DEG=FindMarkers(EPC,ident.1=3,ident.2=c(17,2,8,6,11,16,15),pct.min=0.25)
C3DEG=C3DEG[C3DEG$p_val_adj<0.01,]
C3DEG=C3DEG[abs(C3DEG$avg_log2FC)>1,]
C3High=C3DEG[C3DEG$avg_log2FC>0,]
C3Low=C3DEG[C3DEG$avg_log2FC<0,]
C3HighId=bitr(rownames(C3High),"SYMBOL","ENTREZID",OrgDb="org.Hs.eg.db")[,2]
C3LowId=bitr(rownames(C3Low),"SYMBOL","ENTREZID",OrgDb="org.Hs.eg.db")[,2]

#functional enrichment for the DEG lists
upDEGId=list("C0High"=C0HighId,"C10High"=C10HighId,"C13High"=C13HighId,"C3High"=C3HighId)
library(UpSetR)
pdf("upDEGOverlap.pdf")
upset(fromList(upDEGId))
dev.off()
DegKEGG <- compareCluster(geneCluster = upDEGId, fun = "enrichKEGG",minGSSize = 5)

pdf("CancerCellKEGG.pdf",width=10,height=10)
dotplot(DegKEGG,showCategory=20,label_format = 50)
dev.off()
y <- setReadable(DegKEGG, 'org.Hs.eg.db', 'ENTREZID')
write.table(data.frame(y),file="CancerCellKEGG.txt",sep="\t",quote=F)

upDEGIdGO <- compareCluster(geneCluster = upDEGId, fun = "enrichGO",OrgDb="org.Hs.eg.db",minGSSize = 5,ont="BP")
bp2 <- simplify(upDEGIdGO, cutoff=0.7, by="p.adjust", select_fun=min)
pdf("CancerCellUpGOBPbp2.pdf",width=10,height=5)
dotplot(bp2,showCategory=5,label_format=100)
dev.off()
y <- setReadable(bp2 , 'org.Hs.eg.db', 'ENTREZID')
write.table(data.frame(y),file="CancerCellUpGOBPbp2.txt",sep="\t",quote=F)

downDEGId=list("C0Low"=C0LowId,"C10Low"=C10LowId,"C13Low"=C13LowId,"C3Low"=C3LowId)
downDEGIdGO <- compareCluster(geneCluster = downDEGId, fun = "enrichGO",OrgDb="org.Hs.eg.db",minGSSize = 5,ont="BP")
bp2 <- simplify(downDEGIdGO)
pdf("CancerCellDownGOBPbp2.pdf",width=10,height=5)
dotplot(downDEGIdGO,showCategory=5,label_format=100)
dev.off()


C0DEG=FindMarkers(EPC,ident.1=c(0),ident.2=c(17,2,8,6,11,16,15),pct.min=0.25)

#Quantify the affect of TFs on DEG list
#TF analysis were performed using SCENIC: TranscriptionFactorBySCENIC.code.R
ELF1=read.table("/Projects/deng/Rectal/SCENIC/TargetGene/ELF1.txt")
ELF3=read.table("/Projects/deng/Rectal/SCENIC/TargetGene/ELF3.txt")
MYC=read.table("/Projects/deng/Rectal/SCENIC/TargetGene/MYC.txt")
PPARG=read.table("/Projects/deng/Rectal/SCENIC/TargetGene/PPARG.txt")
HNF4A=read.table("/Projects/deng/Rectal/SCENIC/TargetGene/HNF4A.txt")
TargetGene=unique(c(ELF1[,1],ELF3[,1],MYC[,1],PPARG[,1],HNF4A[,1]))
library(ggrepel)
hs_data=C0DEG
hs_data$ID=rownames(hs_data)
hs_data$Target = as.factor(ifelse(hs_data$ID %in% TargetGene, "Target",'NoTarget'))
table(hs_data$Target)
hs_data$threshold = as.factor(ifelse(hs_data$p_val_adj < 0.01 & abs(hs_data$avg_log2FC)>1, ifelse(hs_data$avg_log2FC > 1,"Up",'Down'),'Others'))
table(hs_data$threshold)
t=ggplot(data = hs_data, aes(x = avg_log2FC, y = -log10(p_val_adj), colour=Target, label =ID )) +
  geom_point(alpha=0.4) +
  theme_bw() + 
  scale_color_manual(values=c("grey","red")) +
  #geom_vline(xintercept=c(-log2(1.5),log2(1.5)),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2FC",y="-log10 (p-value)",title="Differential Expressed Genes") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
    geom_text_repel(
    data = subset(hs_data, abs(hs_data$avg_log2FC)>1 & hs_data$p_val_adj<0.01),
    aes(label = ID),
    size = 3,
    max.overlaps=10,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
pdf("Cluster0DEGTarget.pdf",width=5,height=4)
print(t)
dev.off()


TFs=c("ELF3","CTNNB1","PPARG","ESRRA","SREBF1","HNF4A","ELF1","MYC","POLE4")
for(tf in TFs){
tableSubset <- regulonTargetsInfo[TF %in% tf & highConfAnnot==TRUE,]
tarGenes=tableSubset$gene
#overlap=intersect(tarGenes,immune[,1])
overlap=intersect(tarGenes,rownames(SCPCDEG))
print(paste0(tf,"_",length(overlap),sep=""))
}
[1] "ELF3_8"
[1] "CTNNB1_0"
[1] "PPARG_0"
[1] "ESRRA_0"
[1] "SREBF1_0"
[1] "HNF4A_1"
[1] "ELF1_18"
[1] "MYC_2"
[1] "POLE4_0"
for(tf in TFs){
tableSubset <- regulonTargetsInfo[TF %in% tf & highConfAnnot==TRUE,]
tarGenes=tableSubset$gene
write.table(tarGenes,file=paste0("TargetGene/",tf,".txt"),sep="\t",quote=F,row.names=F,col.names=F)
}


##visualization the specific genes between clusters
library(viridis)
NeutrophilsGene=c("FTL","CTSH","PKM","CPNE1","CEACAM6","PPIA","HSP90AB1","PA2G4","MGST1","NME2","S100A11","PSMB7","MIF","SERPINB1","EEF2","S100A9","LCN2","GOLGA7","LAMP2","PDXK","CSTB","CCT2","PSMA2","NIT2","RAB7A","VAT1","PSMD11","PSMD12","XRCC5","ATP6V1D","TMBIM1","PSMC2","ATP6AP2","CCT8","CNN2","KPNB1","CD47","RAB14","ARPC5","PSMD2","RAB5C","FUCA2","ILF2","GNS","LGALS3","CYB5R3","PRDX4","RHOA","CREG1","PSMD14","KCMF1","ERP44","CD55","SDCBP","PSMB1","RAC1","HLA-B","TIMP2","ANXA2","ANXA3","CXCL1")
t=DotPlot(EPC,features=NeutrophilsGene)
data=t$data
clusterOrder=factor(data$id,levels=c(0,10,13,3,17,2,8,6,11,16,15,1,7,12,14,9,4,5))
g=ggplot(data, aes(clusterOrder,factor(features.plot,levels=rev(unique(features.plot))), size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_colour_viridis_c()+
labs(color="Average Expression",size="Percent expressed",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0),face="bold"),axis.text.x = element_text(size=rel(1.5),face="bold",angle=90),axis.text.y = element_text(size=rel(1.0),face="bold")) 
pdf("Neutrophils.Gene.pdf",height=10)
print(g)
dev.off()

##DEG analysis between tumor cells from SCs and PCs
SCPCDEG=FindMarkers(EPC,ident.1=c(0,10,13),ident.2=c(3),pct.min=0.25)
SCPCDEG=SCPCDEG[SCPCDEG$p_val_adj<0.01,]
SCPCDEGTmp=SCPCDEG[abs(SCPCDEG$avg_log2FC)>1,]
SPHigh=SCPCDEGTmp[SCPCDEGTmp$avg_log2FC>0,]
SPLow=SCPCDEGTmp[SCPCDEGTmp$avg_log2FC<0,]
SPHighId=bitr(rownames(SPHigh),"SYMBOL","ENTREZID",OrgDb="org.Hs.eg.db")[,2]
SPLowId=bitr(rownames(SPLow),"SYMBOL","ENTREZID",OrgDb="org.Hs.eg.db")[,2]
GeneList=list("SPHighId"=SPHighId,"SPLowId"=SPLowId)
SPDEGIdGO <- compareCluster(geneCluster = GeneList, fun = "enrichGO",OrgDb="org.Hs.eg.db",minGSSize = 5,ont="BP")
bp2 <- simplify(SPDEGIdGO, cutoff=0.7, by="p.adjust", select_fun=min)
pdf("SPGOBPbp2.pdf",width=10,height=5)
dotplot(bp2,showCategory=10,label_format=100)
dev.off()

