
library(infercnv)
library(Seurat)
EPC=readRDS("/Projects/deng/Rectal/EPC.rds")
table(EPC$seurat_clusters)
#EPC=subset(EPC,idents=c(0:17))
counts_matrix = as.matrix(EPC@assays$RNA@counts)
pData=EPC@meta.data$seurat_clusters
names(pData)=rownames(EPC@meta.data)
pData=data.frame(pData)
annotation=read.table("/Projects/deng/refGene/gencode.v38.annotation.removeDupliSymbol.geneType.txt",header=TRUE,row.names=1)
anno=data.frame("Chr"=annotation$Chr,"Start"=annotation$Start,"End"=annotation$End)
rownames(anno)=annotation$Symbol
genes=intersect(rownames(counts_matrix),rownames(anno))
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix[genes,],
                                      annotations_file=pData,
                                      delim="\t",
                                      gene_order_file=anno[genes,],
                                      ref_group_names=c("17","2","8","6","11","16","15"))

infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics   
                               denoise=TRUE,
                               cluster_by_groups = TRUE,
                               HMM=TRUE,
                               out_dir="/Projects/deng/Rectal/InferCNVHmm/"
                               )
#There are 7957 genes and 5637 cells remaining in the expr matrix.
