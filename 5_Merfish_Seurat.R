######################################################################################
################################ Detail of the script ################################
######################################################################################

# This script aims at identifying cDC1 cell from segmented cells identified in the publicly available dataset
# https://www.nature.com/articles/s41591-023-02345-0. This dataset used MERFISH to identify niches and gene modules 
# enriched in specific regions of the tumour. Inputs are matrices of segmented cells and genes expression for 6 samples.
# 1. Create Seurat object for each matrix
# 2. Normalize, scale, Run PCA and UMAP, cluster using different resolutions
# 3. Export DEGs, Featureplots to identify clusters of cDC1 to export in future analyses.

# Follow up analyses: use the metadata of the identified cDC1s and split by CCR7 and CXCL9 expression. 
# Generate masks based on densities of other cell types and calculate distances of the CCR7 and CXCL9 cDC1 population
# to thoses masks

library(Seurat)
options(bitmapType='cairo')


ff <- as.list(list.files(path="../7758080/gene/"))


myfilelist <- lapply(ff, function(l){
  patient_temp <- strsplit(as.character(l), split="_")[[1]][1:3]
  dataset <- paste0(patient_temp[1],patient_temp[2],patient_temp[3],"_")
  print(dataset)
  table <- read.csv(paste0("../7758080/gene/",l), row.names = 1)
  
  row.names(table) <- paste0(dataset, row.names(table))
  table_transpose <- t(table)
  
  test <- CreateSeuratObject(counts = table_transpose)
  test <- subset(test, subset = nFeature_RNA > 10)
  return(test)
  
})
names(myfilelist) <- list.files(path="../7758080/gene/", full.names=FALSE)
saveRDS(myfilelist,"Srt_object_after_creation.RDS")

ifnb.list <- myfilelist


dir.create("Louvain_res_test")
dir.create("Metadata_Seurat")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(pbmc_A4) {
  
  sanmple_in_the_loop <- names(ifnb.list)[parent.frame()$i[]]
  
  pbmc_A4 <- NormalizeData(pbmc_A4)
  pbmc_A4 <- FindVariableFeatures(pbmc_A4, selection.method = "vst",  nfeatures = 2000)
  pbmc_A4 <- ScaleData(pbmc_A4, verbose = FALSE)
  
  pbmc_A4 <- RunPCA(pbmc_A4, npcs = 50, verbose = FALSE)
  pbmc_A4 <- RunUMAP(pbmc_A4, reduction = "pca", dims = 1:50)
  
  pbmc_A4 <- FindNeighbors(pbmc_A4, reduction = "pca", dims = 1:50)
  
  vector_res <- seq(1, 1.4, by=0.2) 
  for(i in vector_res) {
    pbmc_A4 <- FindClusters(pbmc_A4, resolution = i)
  }
  
  #dir.create("Louvain_res_test")
  for(i in vector_res) {
    name_to_use <- paste0("RNA_snn_res.",i)
    Idents(pbmc_A4) <- name_to_use
    pdf(paste0("Louvain_res_test/",sanmple_in_the_loop,"_Louvain_",name_to_use,".pdf"))
    print(DimPlot(pbmc_A4, label = T))
    dev.off()
  }
  
 
  Idents(pbmc_A4) <- "RNA_snn_res.1"
  pbmc_A4_markers <- FindAllMarkers(pbmc_A4, only.pos = TRUE,
                                    min.pct = 0.2, logfc.threshold = 0.2)

  write.csv(pbmc_A4_markers,paste0("DEG_res_1/",sanmple_in_the_loop,"_DEG.csv"))
  write.csv(pbmc_A4@meta.data,paste0("Metadata_Seurat/",sanmple_in_the_loop,"_metadata_Srt.csv"))
  #saveRDS(pbmc_A4, file =paste0("Srt_Obj/",sanmple_in_the_loop,"_SrtObj.rds"))
  return(pbmc_A4)
  
})
saveRDS(ifnb.list, file = "list_all_Srt_object.rds")



ifnb.list <- readRDS("list_all_Srt_object.rds")


dir.create("Featureplots")


#plot specific genes
lapply(X = ifnb.list, FUN = function(pbmc_A4) {
  sanmple_in_the_loop <- names(ifnb.list)[parent.frame()$i[]]
  print(sanmple_in_the_loop <- names(ifnb.list)[parent.frame()$i[]])
  pdf(file=paste0("Featureplots/",sanmple_in_the_loop,"_.pdf"))
  print(FeaturePlot(pbmc_A4,"CD74", order = T))
  dev.off()
  
  
})

dir.create("Dotplot")


lapply(X = ifnb.list, FUN = function(pbmc_A4) {
  sanmple_in_the_loop <- names(ifnb.list)[parent.frame()$i[]]
  print(sanmple_in_the_loop <- names(ifnb.list)[parent.frame()$i[]])
  
  # pdf(file=paste0("DimPlot/",sanmple_in_the_loop,"_DimPlot.pdf"))
  # print(DimPlot(pbmc_A4,group.by ="RNA_snn_res.1.1", label = T ))
  # dev.off()
  
  pdf(file=paste0("Dotplot/",sanmple_in_the_loop,"_.pdf"))
  print(DotPlot(pbmc_A4,features = c("CXCL9","CCR7"),group.by = "RNA_snn_res.1"))
  dev.off()
  
  
})

dir.create("DEG_res_1.4")

lapply(X = ifnb.list, FUN = function(pbmc_A4) {
  sanmple_in_the_loop <- names(ifnb.list)[parent.frame()$i[]]
  Idents(pbmc_A4) <- "RNA_snn_res.1.4"
  pbmc_A4_markers <- FindAllMarkers(pbmc_A4, only.pos = TRUE,
                                    min.pct = 0.2, logfc.threshold = 0.2)
  
  write.csv(pbmc_A4_markers,paste0("DEG_res_1.4/",sanmple_in_the_loop,"_DEG.csv"))
  
})

lapply(X = ifnb.list, FUN = function(pbmc_A4) {
  sanmple_in_the_loop <- names(ifnb.list)[parent.frame()$i[]]
  Idents(pbmc_A4) <- "RNA_snn_res.1"
  pbmc_A4_markers <- FindAllMarkers(pbmc_A4, only.pos = TRUE,
                                    min.pct = 0.1, logfc.threshold = 0.1)
  
  write.csv(pbmc_A4_markers,paste0("DEG_res_1/",sanmple_in_the_loop,"0.1Threshold_DEG.csv"))
  
})

ifnb.list <- readRDS("../Seurat_200K_filterout_take2/list_all_Srt_object.rds")


###########################################################################################
##################### Add CXCL9 and CCR7 expression to the metadata #######################
###########################################################################################


ifnb.list <- readRDS("list_all_Srt_object.rds")


#dir.create("Featureplots")

# Save CCR7 and CXCL9 expression together with metadata
lapply(X = ifnb.list, FUN = function(pbmc_A4) {
  
  sanmple_in_the_loop <- names(ifnb.list)[parent.frame()$i[]]
  
  test <- pbmc_A4@assays$RNA$counts[row.names(pbmc_A4@assays$RNA$counts) %in% c("CXCL9", "CCR7"),]
  test2 <- t(as.matrix(test))
  
  metadata <- pbmc_A4@meta.data
  
  Merged <- merge(metadata,test2,by= "row.names" )
  
  write.csv(Merged, paste0("Metadata_Seurat/",sanmple_in_the_loop,"_metadata_Srt_With_Cxcl9_CCR7_expression.csv"), row.names = F)
  
})




