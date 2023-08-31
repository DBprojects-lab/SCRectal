
setwd("/Projects/deng/Rectal/DEG")
library(org.Hs.eg.db)
library(Seurat)

#Rectal=readRDS("/Projects/deng/Rectal/Rectal.rds")
Rectal$Category=paste0(Rectal$cellType,"_",Rectal$group,sep="")
Idents(Rectal)=Rectal$Category

setwd("/Projects/deng/Rectal/DEG/SP")
for(i in unique(Rectal$cellType)){
   CaseGroup=paste(i,"_SC",sep="");
   CtlGroup=paste(i,"_PC",sep="");
   if(is.na(table(Rectal$Category)[CaseGroup])){
      next
    }
   if(is.na(table(Rectal$Category)[CtlGroup])){
      next
   }
   if(table(Rectal$Category)[CaseGroup]<10){
      next
   }
   if(table(Rectal$Category)[CtlGroup]<10){
      next
   }
   DEG <- FindMarkers(Rectal, ident.1 = CaseGroup, ident.2 = CtlGroup, verbose = FALSE,min.pct = 0.25)
   DEG=DEG[DEG$ p_val_adj<0.01,]
   write.table(DEG,file=paste("allGene/",i,"Deg.txt",sep=""),sep="\t",quote=F)

   UpDEG=DEG[DEG$avg_log2FC>0,]
   write.table(UpDEG,file=paste("upGene/",i,"Deg.txt",sep=""),sep="\t",quote=F)
   degId=bitr(rownames(UpDEG),"SYMBOL","ENTREZID",OrgDb="org.Hs.eg.db")[,2]
   Function=enrichGO(degId,OrgDb="org.Hs.eg.db",minGSSize = 5,ont="BP")
   if(!is.null(Function)){
    y <- setReadable(Function, 'org.Hs.eg.db', 'ENTREZID')
    write.table(data.frame(y),file=paste("upGene/",i,"Function.txt",sep=""),sep="\t",quote=F)
    }

   DownDEG=DEG[DEG$avg_log2FC<0,]
   write.table(DownDEG,file=paste("downGene/",i,"Deg.txt",sep=""),sep="\t",quote=F)
   degId=bitr(rownames(DownDEG),"SYMBOL","ENTREZID",OrgDb="org.Hs.eg.db")[,2]
   Function=enrichGO(degId,OrgDb="org.Hs.eg.db",minGSSize = 5,ont="BP")
   if(!is.null(Function)){
    y <- setReadable(Function, 'org.Hs.eg.db', 'ENTREZID')
    write.table(data.frame(y),file=paste("downGene/",i,"Function.txt",sep=""),sep="\t",quote=F)
 }
}



