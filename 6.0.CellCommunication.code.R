library(CellChat)
setwd("/Projects/deng/Rectal/CellChat")
Rectal=readRDS("/Projects/deng/Rectal/Rectal.rds")
EPC=readRDS("/Projects/deng/Rectal/EPC.rds")
RectalImmuneCell=subset(Rectal,idents=c("Macrophages","Mast cells","Neutrophils","T cells","Naive B cells","Plasma B cells"))

new.cluster.ids <- c("T cells","B cells","Macrophages","B cells","Neutrophils","Mast cells")
names(new.cluster.ids) <- levels(RectalImmuneCell)
RectalImmuneCell <- RenameIdents(RectalImmuneCell, new.cluster.ids)
RectalImmuneCell$Type=Idents(RectalImmuneCell)

EPCTmp=subset(EPC,idents=c(0:17))
new.cluster.ids <- c("SC","Normal","Normal","PC","Normal","Normal","Normal","Normal","Normal","Normal","SC","Normal","Normal","SC","Normal","Normal","Normal","Normal")
names(new.cluster.ids) <- levels(EPCTmp)
EPCTmp <- RenameIdents(EPCTmp, new.cluster.ids)
EPCTmp$Type=Idents(EPCTmp)

Rectal.combined <- merge(RectalImmuneCell, y =c(EPCTmp) , add.cell.ids = c("Immune", "EPC"), project = "AllCells")
table(Rectal.combined$group)
table(Rectal.combined$Type)

exprMat  <-  as.matrix(Rectal.combined@assays$RNA@data)
meta = Rectal.combined@meta.data
cellchat <- createCellChat(object = exprMat, meta = meta, group.by = "Type")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)
#dplyr::glimpse(CellChatDB$interaction)

options(stringsAsFactors = FALSE)
interaction_input <- read.table(file = 'Interaction8CellChatDB.txt', row.names = 1,header=T,sep="\t") #add the celltalkDB into database
MyCellChatDB <- list()
MyCellChatDB$interaction <- interaction_input
MyCellChatDB$complex <- CellChatDB$complex
MyCellChatDB$cofactor <- CellChatDB$cofactor
MyCellChatDB$geneInfo <- CellChatDB$geneInfo
#dplyr::glimpse(MyCellChatDB$interaction)
#showDatabaseCategory(MyCellChatDB)
table(MyCellChatDB$interaction$annotation)
Cell-Cell Contact       ECM-Receptor Secreted Signaling
               319                421               4040
#CellChatDB.use <- CellChatDB
CellChatDB.use <- MyCellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = "cellchat8CellChatDB.rds")
