##Generate the Figure S2B and S2C
setwd("/Projects/deng/Rectal/EPC")
#https://cytotrace.stanford.edu/
devtools::install_local("/Softwares/deng/CytoTRACE_0.3.3.tar.gz")
library(CytoTRACE)
EPC=readRDS("/Projects/deng/Rectal/EPC.rds")
EPC=subset(EPC,idents=c(0:17))

EPCExprMat  <-  as.matrix(EPC@assays$RNA@data)
results <- CytoTRACE(EPCExprMat, ncores = 8, subsamplesize = 1000)
phenotype=data.frame("Cluster"=paste0("C",EPC$seurat_clusters,sep=""))
rownames(phenotype)=rownames(EPC@meta.data)
all(names(results$CytoTRACE)==rownames(phenotype))
diffInfo=data.frame("cellNames"=rownames(phenotype),"Cluster"=phenotype$Cluster,"DifPotential"=results$CytoTRACE)
ClusterOrder=factor(diffInfo$Cluster,levels=paste0("C",c(17,2,8,6,11,16,15,1,7,12,3,14,9,13,4,5,0,10)))
g=ggplot(diffInfo, aes(x=ClusterOrder, y=DifPotential, fill=Cluster))+
geom_boxplot()+labs(title="",x="", y = "")+
theme_classic()+
theme(legend.position="none")+
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))+
theme(axis.title.y = element_text(size=rel(1.5)),axis.text.x = element_text(size=rel(1.5),face="bold",angle=90),axis.text.y = element_text(size=rel(1.5),face="bold"))
pdf("DifPotential8CytoTRACE.pdf",width=15)
print(g)
dev.off()




library(SCENT)
#download from NCBI, for symbol and gene ID convert
geneInfo=read.table("/Projects/deng/Public/NCBI/HumanSymbol2ID.txt",check.names=F)
geneInfo=geneInfo[!duplicated(geneInfo[,3]),]
rownames(geneInfo)=geneInfo[,3]
gene=intersect(rownames(EPCExprMat),rownames(geneInfo))
EPCExprMat=EPCExprMat[gene,]
geneInfo=geneInfo[gene,]
all(rownames(geneInfo)==rownames(EPCExprMat))
rownames(EPCExprMat)=geneInfo[,2]

data(net17Jan16);
range(EPCExprMat);
ccat.v <- CompCCAT(exp = EPCExprMat, ppiA = net17Jan16.m);
all(colnames(EPCExprMat)==rownames(EPC@meta.data))

boxplot(ccat.v ~ colnames(EPCExprMat),xlab="Embryonic stages",ylab="CCAT",col=colorLiver.v);

result=data.frame(cbind("CellName"=colnames(EPCExprMat),"Cluster"=paste0("C",EPC$seurat_clusters),"Entropy"=ccat.v))
write.table(result,file="SCENTResult.txt",sep="\t",quote=F,row.names=F)
result=read.table("SCENTResult.txt",header=T)

g1=ggplot(result, aes(x=ClusterOrder, y=Entropy, fill=Cluster))+
geom_boxplot()+labs(title="",x="", y = "")+
theme_classic()+
theme(legend.position="none")+
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))+
theme(axis.title.y = element_text(size=rel(1.5)),axis.text.x = element_text(size=rel(1.5),face="bold",angle=90),axis.text.y = element_text(size=rel(1.5),face="bold"))
pdf("DifPotential8SCENT.pdf",width=15)
print(g1)
dev.off()



#cell cycle analysis
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
EPC <- CellCycleScoring(EPC, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
pdf("cellCyclingCell.pdf",height=4,width=6)
TSNEPlot(EPC, group.by="Phase")&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
