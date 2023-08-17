library(Seurat)
library(tidyverse)
library(ggplot2)
library(scCustomize)
library(gridExtra)
library(DESeq2)
library(celldex)
library(ggpubr)
library(SingleR)
library(scMRMA)
library(SeuratWrappers)
library(Nebulosa)
library(dittoSeq)
library(harmony)
library(cowplot)
library(viridis)
library(scuttle)
library(trqwe)
library(GEOquery)
setwd("D:/UL/dissertation/CD4/code")

GEO_Num = "GSE134520"
project = "GASTRIC"
GEO_Num_D = paste("./", GEO_Num,"/", sep = '')
# getGEOSuppFiles(GEO_Num) # getting the gastric 
untar(paste(GEO_Num_D, GEO_Num, '_RAW.tar', sep = ''), exdir = GEO_Num_D)
files_list = untar(paste(GEO_Num_D, GEO_Num, '_RAW.tar', sep = ''),list=TRUE)
files_list
func <- function(f)
{
  gastric.p1.sparse.m <- read.table(paste(GEO_Num_D, f, sep=''), sep="\t",header=TRUE)
  gastric.p1.seurat.obj <- CreateSeuratObject(counts = gastric.p1.sparse.m, min.cells = 3, min.genes = 300, project = project)
}
seurat.obj.list <- c()
for (file in files_list){
  seurat.obj.list <- append(seurat.obj.list, func(file))
}
b <- c("patient1_", "patient2_", "patient3_", "patient4_", "patient5_", "patient6_",
       "patient7_", "patient8_", "patient9_", "patient10_", "patient11_", "patient12_",
       "patient13_")
merged_seurat <- merge(seurat.obj.list[[1]], y = (seurat.obj.list)[2:13], add.cell.ids = b, project = 'GASTRIC')


merged_seurat





merged_seurat$sample <- rownames(merged_seurat@meta.data)
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Patient', 'Barcode'), sep = '__')
values <- c("NAG", "NAG", "NAG", "CAG", "CAG", "CAG", "IMW", "IMW", "IMS", "IMS", "IMS", "IMS", "EGC")
merged_seurat$type <- values[match(merged_seurat$Patient, b)]
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern="^MT-")
merged_seurat$riboPercent <- PercentageFeatureSet(merged_seurat, pattern="^RPL")
merged_seurat_filtered <- subset(merged_seurat, subset = nFeature_RNA > 400 & nFeature_RNA < 7000 & mitoPercent < 20 & riboPercent < 20)
merged_unprocessed <- merged_seurat_filtered

merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.05, x.high.cutoff = 5, y.cutoff = 0.5, do.plot = T)
merged_seurat_filtered <- ScaleData(merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(merged_seurat_filtered, verbose = FALSE)


merged_seurat_filtered
merged_seurat_Aggregated <- AggregateExpression(  merged_seurat_filtered, slot = "count")





merged_seurat_filtered <- RunTSNE(merged_seurat_filtered, dims = 1:20, reduction.name ="UMAP_BatchEffect_Uncorrected")
merged_seurat_filtered <- RunUMAP(merged_seurat_filtered, dims = 1:20, reduction.name ="UMAP_BatchEffect_Uncorrected")
ElbowPlot(merged_seurat_filtered)
# merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
# merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)

dittoDimPlot(merged_seurat_filtered, var="Patient", reduction.use = "UMAP_BatchEffect_Uncorrected", do.label = TRUE) + ggtitle("Batch effects")








merged_mnn <- merged_unprocessed
merged_mnn <- NormalizeData(object = merged_mnn, normalization.method = "LogNormalize", scale.factor = 10000)
merged_mnn <- FindVariableFeatures(object = merged_mnn, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.05, x.high.cutoff = 5, y.cutoff = 0.5, do.plot = T)
merged_seurat_filtered_mnn <- RunFastMNN(object.list = SplitObject(merged_mnn, split.by = "Patient"))
merged_seurat_filtered_mnn <- RunUMAP(merged_seurat_filtered_mnn, reduction = "mnn", dims = 1:20,  reduction.name = "UMAP_mnnCorrected")
merged_seurat_filtered_mnn <- FindNeighbors(merged_seurat_filtered_mnn, reduction = "mnn", dims = 1:20)
merged_seurat_filtered_mnn <- FindClusters(merged_seurat_filtered_mnn)
dittoDimPlot(merged_seurat_filtered_mnn, "Patient", reduction.use = "UMAP_mnnCorrected", do.label = TRUE) + ggtitle("Batch correction using MNN")







library(SingleR)
library(celldex)





seurat.singleR <- merged_seurat_filtered_mnn
dataset.ref <- BlueprintEncodeData()
seurat.singleR.sce <- as.SingleCellExperiment(seurat.singleR)
singleR.labels <- SingleR(test = seurat.singleR.sce ,assay.type.test = 1, ref = dataset.ref, labels = dataset.ref$label.main)
seurat.singleR@meta.data$singleR.labels <- singleR.labels$pruned.labels
seurat.singleR <- SetIdent(seurat.singleR, value = "singleR.labels")

write.csv(seurat.singleR@meta.data, 'seurat_singleR.csv')

all.cluster.markers <- FindAllMarkers(seurat.singleR, logfc.threshold = log2(1.2), min.pct = 0.5)
write.csv(all.cluster.markers, "ISCBgastric_singleR_labels_mnn.csv")
UMAP_singleR <- dittoDimPlot(seurat.singleR, reduction.use = 'UMAP_mnnCorrected', var = 'singleR.labels') + ggtitle("SingleR Annotations") + theme(legend.text = element_text(size = 10), aspect.ratio = 1)
ggarrange(UMAP_singleR)





merged_seurat_scMRMA <- merged_seurat_filtered_mnn
result <- scMRMA(input = merged_seurat_scMRMA,species = "Hs", db = "panglaodb", p = 0.05, normalizedData = F, selfDB = NULL, selfClusters = NULL, k=20)
merged_seurat_scMRMA[["scMRMA"]] <- result$multiR$annotationResult[colnames(merged_seurat_scMRMA), ncol(result$multiR$annotationResult)]
UMAP_scMRMA <- dittoDimPlot(merged_seurat_scMRMA,reduction.use = "UMAP_mnnCorrected",var = "scMRMA")+ ggtitle("scMRMA Annotations") + theme(legend.text = element_text(size = 10), aspect.ratio = 1)
ggarrange(UMAP_scMRMA)
merged_seurat_scMRMA <- SetIdent(merged_seurat_scMRMA, value = "scMRMA")


write.csv(merged_seurat_scMRMA@meta.data, 'seurat_scMRMA.csv')

all.cluster.markers <- FindAllMarkers(merged_seurat_scMRMA, logfc.threshold = log2(1.2), min.pct = 0.5) 
write.csv(all.cluster.markers, "ISCBgastric_scMRMA_panglaodb_labels_mnn.csv")




result <- scMRMA(input = merged_seurat_scMRMA, species = "Hs", db = "TcellAI", p = 0.05, normalizedData = F, selfDB = NULL, selfClusters = NULL, k=20)
merged_seurat_scMRMA[["scMRMA"]] <- result$multiR$annotationResult[colnames(merged_seurat_scMRMA), ncol(result$multiR$annotationResult)]
UMAP_scMRMA_tcellAI <- dittoDimPlot(merged_seurat_scMRMA,reduction.use = "UMAP_mnnCorrected",var = "scMRMA_TcellAI", do.label = TRUE)+ ggtitle("scMRMA Annotations with TcellAI") + theme(legend.text = element_text(size = 10), aspect.ratio = 1)
ggarrange(UMAP_scMRMA_tcellAI)
merged_seurat_scMRMA <- SetIdent(merged_seurat_scMRMA, value = "scMRMA")
all.cluster.markers <- FindAllMarkers(merged_seurat_scMRMA, logfc.threshold = log2(1.2), min.pct = 0.5) 
write.csv(all.cluster.markers, "ISCBgastric_scMRMA_TcellAI_labels_mnn.csv")












