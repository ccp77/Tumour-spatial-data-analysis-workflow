library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
batch_name <- "IntegrationTake14"

################################################################################################################################
######################################### Details of the script ################################################################
################################################################################################################################

# This script builds on the Seurat pipelines to integrate three samples of dendritic cell type 1 (cDC1s) in tumours.
# More information on the pipeline can be found on the Seurat website: https://satijalab.org/seurat/articles/integration_introduction
# 1. Selects cells that are to be kept for this analysis by building what has been discovered in previous steps.
# 2. Generate seurat objects for the three samples A4, A5 and A6, end keep the cell IDs selected in the step 1
# 3. Proceed to Seurat analysis. Integrate to remove batch effect, scale, run PCA and UMAP, cluster using different resolutions (fine tune Louvain)
#   Select resolution, generate DEGs. Generate Heatmaps, Violoin plots, Feature plots.


################################################################################################################################
################################ 1. Select the cell IDs to keep all the cDC1s of the dataset ################################
################################################################################################################################
# cDC1s tumours: This selects cells that were identified as cDC1 in the previous analysis called Take4
Take4_cDC1 <- readRDS( file = "../Seurat_take4_cDC1/IntegrationTake4_cDC1__afterIntegratedAnalysis_Louvain_test_Object5.rds")
Idents(Take4_cDC1) <-"integrated_snn_res.0.5"
DimPlot(Take4_cDC1, label = T) 
Take4_cDC1_no_cluster8 <- subset(Take4_cDC1, idents = ("8"), invert=T)
DimPlot(Take4_cDC1_no_cluster8)
CellIDs_cDC1 <- colnames(Take4_cDC1_no_cluster8@assays$RNA) #5468 cells

# cDC1s in activated cDCs: : This selects cells that were identified as cDC1 in the previous analysis called Take_activatedT13
Take_activated13 <- readRDS(file = "../Seurat_take13_activated_cDCs/IntegrationTake13_afterIntegratedAnalysis_SrtObject4.rds")
Idents(Take_activated13) <- "integrated_snn_res.0.5_CP_Annot"
DimPlot(Take_activated13)
Take_activated13_cDC1s <- subset(Take_activated13,idents = c("6_cDC1s","4","3","1","2"))
CellIDs_Activated_cDC1 <- colnames(Take_activated13_cDC1s@assays$RNA) #969 cells

CellIDs_Activated_not_doublon <- CellIDs_Activated_cDC1[!CellIDs_Activated_cDC1 %in% CellIDs_cDC1] # 338 cells
DimPlot(Take_activated13,cells.highlight =CellIDs_Activated_not_doublon )

# compile all cell IDs of cDC1s in tumours
CellIDs_cDC1s_and_activated_cDC1s <- c(CellIDs_cDC1,CellIDs_Activated_not_doublon)

Take3Bis <-  readRDS( file = "../Seurat_take3_BIS/IntegrationTake3_BIS__afterIntegratedAnalysis_SrtObject_res0.4.rds")
DimPlot(Take3Bis,cells.highlight=CellIDs_cDC1s_and_activated_cDC1s) # cell Ids are correct

################################################################################################################
################# 2. Generate new seurat object with the cDC1s and activated cDC1s ###########################
################################################################################################################


pbmc.data_A4 <- Read10X(data.dir = "../Mapping_A4/Cellranger_output/merging_eGFP_mapping_SC21125_PIO3801A4_v2/outs/filtered_feature_bc_matrix/")
pbmc.data_A4[1:5,1:5]
pbmc_A4 <- CreateSeuratObject(counts = pbmc.data_A4, project = "A4", min.cells = 3, min.features = 200)

# Generate new cell names adding A4, A5 or A6
head(colnames(pbmc_A4))
New_names <- paste0("A4_",colnames(pbmc_A4))
head(New_names)

# Rename cells in Seurat object
head(colnames(pbmc_A4))
pbmc_A4 <- RenameCells(pbmc_A4,new.names=New_names)
head(colnames(pbmc_A4))

######################### Object A5 #########################

pbmc.data_A5 <- Read10X(data.dir = "../Mapping_A5/Cellranger_output_A5/merging_eGFP_mapping_SC21125_PIO3801A5/outs/filtered_feature_bc_matrix/")
pbmc.data_A5[1:5,1:5]
pbmc_A5 <- CreateSeuratObject(counts = pbmc.data_A5, project = "A5", min.cells = 3, min.features = 200)

# Generate new cell names adding A4, A5 or A6
head(colnames(pbmc_A5))
New_names <- paste0("A5_",colnames(pbmc_A5))
head(New_names)

# Rename cells in Seurat object
head(colnames(pbmc_A5))
pbmc_A5 <- RenameCells(pbmc_A5,new.names=New_names)
head(colnames(pbmc_A5))

######################### Object A6 #########################

pbmc.data_A6 <- Read10X(data.dir = "../Mapping_A6/Cellranger_output_A6/merging_eGFP_mapping_SC21125_PIO3801A6/outs/filtered_feature_bc_matrix/")
pbmc.data_A6[1:5,1:5]
pbmc_A6 <- CreateSeuratObject(counts = pbmc.data_A6, project = "A6", min.cells = 3, min.features = 200)


# Generate new cell names adding A4, A5 or A6
head(colnames(pbmc_A6))
New_names <- paste0("A6_",colnames(pbmc_A6))
head(New_names)

# Rename cells in Seurat object
head(colnames(pbmc_A6))
pbmc_A6 <- RenameCells(pbmc_A6,new.names=New_names)
head(colnames(pbmc_A6))

################################################################################################
#################### Keep cells cDC1s and activated cDC1s based on the cell IDs generated above ###################
################################################################################################

cellIDs_to_keep_for_take2 <- CellIDs_cDC1s_and_activated_cDC1s

##### CHANGE FOR A4 #####
CellIDs_A4_in_cDC_lineage <- colnames(pbmc_A4@assays$RNA)[colnames(pbmc_A4@assays$RNA) %in% cellIDs_to_keep_for_take2]
pbmc_A4 <- subset(pbmc_A4, cells = CellIDs_A4_in_cDC_lineage)
##### CHANGE FOR A5 #####
CellIDs_A5_in_cDC_lineage <- colnames(pbmc_A5@assays$RNA)[colnames(pbmc_A5@assays$RNA) %in% cellIDs_to_keep_for_take2]
pbmc_A5 <- subset(pbmc_A5, cells = CellIDs_A5_in_cDC_lineage)
##### CHANGE FOR A6 #####
CellIDs_A6_in_cDC_lineage <- colnames(pbmc_A6@assays$RNA)[colnames(pbmc_A6@assays$RNA) %in% cellIDs_to_keep_for_take2]
pbmc_A6 <- subset(pbmc_A6, cells = CellIDs_A6_in_cDC_lineage)


################################################################################################
#################### 3. Proceed to Seurat analysis ######################################
################################################################################################

######################### Combine objects in a list #########################

ifnb.list <- list(pbmc_A4,pbmc_A5,pbmc_A6)
#ifnb.list <- SplitObject(ifnb, split.by = "stim")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)


######################### Integrated analysis #########################

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"
saveRDS(immune.combined, file = paste0(batch_name,"_afterIntegrationSrtObject.rds"))

immune.combined <- readRDS("IntegrationTake14_afterIntegrationSrtObject.rds")

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 50, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:50)
remotes::install_version("Matrix", version = "1.6-1.1")


immune.combined <- readRDS("IntegrationTake14_afterIntegrationSrtObject_UMAP.rds")

immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:50)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# Visualization
pdf(file=paste0(batch_name,"UMAP.pdf"))
DimPlot(immune.combined, reduction = "umap", label=T)
dev.off()
DefaultAssay(immune.combined) <- "RNA"
FeaturePlot(immune.combined,"Cd24a")
VlnPlot(immune.combined,"Cxcl9")

Take3Bis <- readRDS("../Seurat_take3_BIS/IntegrationTake3_BIS__afterIntegratedAnalysis_SrtObject_res0.4_Object4.rds")
DimPlot(Take3Bis,cells.highlight = Cells(subset(immune.combined,idents= 5)))
FeaturePlot(Take3Bis,"Cxcl9")


###################################################################################################
################################# Fine Tune Louvain #################################
###################################################################################################

DefaultAssay(immune.combined) <- "integrated"

vector_res <- seq(0.4, 1.2, by=0.1) 
for(i in vector_res) {
  immune.combined <- FindClusters(immune.combined, resolution = i)
}

dir.create("Louvain_res_test")
for(i in vector_res) {
  name_to_use <- paste0("integrated_snn_res.",i)
  Idents(immune.combined) <- name_to_use
  pdf(paste0("Louvain_res_test/",batch_name,"Louvain_,",name_to_use,".pdf"))
  print(DimPlot(immune.combined, label = T))
  dev.off()
}

saveRDS(immune.combined, file = "IntegrationTake14_afterIntegrationSrtObject_UMAP_2.rds")

## CHoose resolution 0.8
Idents(immune.combined) <- "integrated_snn_res.0.8"
DimPlot(immune.combined,label = T)

Take3Bis <- readRDS("../Seurat_take3_BIS/IntegrationTake3_BIS__afterIntegratedAnalysis_SrtObject_res0.4_Object4.rds")
DimPlot(Take3Bis,cells.highlight = Cells(subset(immune.combined,idents= 7)))
VlnPlot(immune.combined,"Cd40")

##############################################################################
########################## DEG for Louvain resolution 0.8 ##########################
##############################################################################
immune.combined <- readRDS("IntegrationTake14_afterIntegrationSrtObject_UMAP_2.rds")
DefaultAssay(immune.combined) <- "RNA"

Idents(immune.combined) <- "integrated_snn_res.0.8"
DimPlot(immune.combined, label = T)
immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE,
                                          min.pct = 0.2, logfc.threshold = 0.2)
dir.create("DEG_all_res0.8")
write.csv(immune.combined.markers,"DEG_all_res0.8/FindAllMarkers_res0.8.csv",
          row.names = T, col.names = T)


DEG <- read.csv("DEG_all_res0.8/FindAllMarkers_res0.8.csv", row.names = 1)

DefaultAssay(immune.combined) <-"integrated"

DEG %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC) -> top10
head(top10)
DoHeatmap(immune.combined, features = top10$gene) + NoLegend()
dev.off()

######################## Save Metadata ########################
Metadata_take14 <- immune.combined@meta.data
head(Metadata_take14)
Metadata_take14$Seurat_cluster_to_use <- immune.combined@meta.data$integrated_snn_res.0.8
UMAP_coord_take_14 <- immune.combined@reductions$umap@cell.embeddings
head(UMAP_coord_take_14)
Merged <- merge(Metadata_take14,UMAP_coord_take_14, by="row.names")
head(Merged)

write.csv(Merged,"Metadata_UMAP_take14_all_cDCs.csv")

##################################################################################################
################################ DPlot heatmap of all cDC1s ################
###############################################################################################


immune.combined <- readRDS("IntegrationTake14_afterIntegrationSrtObject_UMAP_2.rds")
DefaultAssay(immune.combined) <- "RNA"
immune.combined.markers <- read.csv("DEG_all_res0.8/FindAllMarkers_res0.8.csv")
head(immune.combined.markers)

######################################## TOP 25 for Heatmap ########################################
immune.combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC) -> top10
head(top10)
levels(as.factor(top10$cluster))
top10_reordered <- rbind(top10[top10$cluster==4,],
                         top10[top10$cluster==1,],
                         top10[top10$cluster==11,],
                         top10[top10$cluster==5,],
                         top10[top10$cluster==0,],
                         top10[top10$cluster==3,],
                         top10[top10$cluster==6,],
                         top10[top10$cluster==2,],
                         top10[top10$cluster==7,],
                         top10[top10$cluster==9,],
                         top10[top10$cluster==10,],
                         top10[top10$cluster==8,])
tail(top10_reordered)

immune.combined <- ScaleData(immune.combined)
levels(as.factor(immune.combined@meta.data$integrated_snn_res.0.8))
Idents(immune.combined) <- "integrated_snn_res.0.8"

levels(immune.combined) <- c("4","1","11","5","0","3","6","2","7","9","10","8")
DoHeatmap(immune.combined, features = top10_reordered$gene) + NoLegend()
genes <- top10_reordered$gene

pdf("Heatmap_reordered_all_cDC1s_top25.pdf")
DoHeatmap(immune.combined, features = top10_reordered$gene) + NoLegend()
dev.off()

write.csv(top10_reordered,"DEG_all_res0.8/Heatmap_genes_DEG_top25.csv")

######################################## TOP 10 for Vln plot ########################################
immune.combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
head(top10)
levels(as.factor(top10$cluster))
top10_reordered <- rbind(top10[top10$cluster==4,],
                         top10[top10$cluster==1,],
                         top10[top10$cluster==11,],
                         top10[top10$cluster==5,],
                         top10[top10$cluster==0,],
                         top10[top10$cluster==3,],
                         top10[top10$cluster==6,],
                         top10[top10$cluster==2,],
                         top10[top10$cluster==7,],
                         top10[top10$cluster==9,],
                         top10[top10$cluster==10,],
                         top10[top10$cluster==8,])


genesVlnPlot <- top10_reordered$gene
write.csv(top10_reordered,"DEG_all_res0.8/DEG_top10.csv")

Genes_to_plot <- read.csv("DEG_all_res0.8/DEG_top10.csv")
head(Genes_to_plot)

pdf("Dotplot_reordered_all_cDC1s_top10.pdf", width = 7.5, height = 3.5)
plot <- DotPlot(object = immune.combined, features = unique(c(Genes_to_plot$gene)),dot.scale = 1.5,)+ NoLegend()
plot + theme(axis.text.x = element_text(angle = 90), axis.text=element_text(size=6))
dev.off()


write.csv(top10_reordered,"DEG_all_res0.8/Heatmap_genes_DEG_top25.csv")


######################################## TOP 5 for Vln plot ########################################
immune.combined <- readRDS("IntegrationTake14_afterIntegrationSrtObject_UMAP_2.rds")
DefaultAssay(immune.combined) <- "RNA"
immune.combined.markers <- read.csv("DEG_all_res0.8/FindAllMarkers_res0.8.csv")
head(immune.combined.markers)

levels(as.factor(immune.combined@meta.data$integrated_snn_res.0.8))
Idents(immune.combined) <- "integrated_snn_res.0.8"

levels(immune.combined) <- c("4","1","11","5","0","3","6","2","7","9","10","8")



immune.combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
head(top10)
levels(as.factor(top10$cluster))
top10_reordered <- rbind(top10[top10$cluster==4,],
                         top10[top10$cluster==1,],
                         top10[top10$cluster==11,],
                         top10[top10$cluster==5,],
                         top10[top10$cluster==0,],
                         top10[top10$cluster==3,],
                         top10[top10$cluster==6,],
                         top10[top10$cluster==2,],
                         top10[top10$cluster==7,],
                         top10[top10$cluster==9,],
                         top10[top10$cluster==10,],
                         top10[top10$cluster==8,])


genesVlnPlot <- top10_reordered$gene
write.csv(top10_reordered,"DEG_all_res0.8/DEG_top5.csv")

Genes_to_plot <- read.csv("DEG_all_res0.8/DEG_top5.csv")
head(Genes_to_plot)

pdf("Dotplot_reordered_all_cDC1s_top5.pdf", width = 7.5, height = 3)
plot <- DotPlot(object = immune.combined, features = unique(c(Genes_to_plot$gene)),dot.scale = 1.6,)+ NoLegend()
plot + theme(axis.text.x = element_text(angle = 90), axis.text=element_text(size=6.5))
dev.off()


write.csv(top10_reordered,"DEG_all_res0.8/Heatmap_genes_DEG_top25.csv")


#### Make VLn plots

#Maturation genes
VlnPlot(immune.combined,"Cxcl9",idents = c("7","9","10","8"),cols = c("darkturquoise","deeppink1","chocolate1","deeppink4"))
VlnPlot(immune.combined,"Cd40",idents = c("1","2","6_cDC1s","3","4"),cols = c("indianred","rosybrown","chocolate1","deeppink","deeppink4"))
VlnPlot(immune.combined,"Cd86",idents = c("1","2","6_cDC1s","3","4"),cols = c("indianred","rosybrown","chocolate1","deeppink","deeppink4"))
VlnPlot(immune.combined,"H2-K1",idents = c("1","2","6_cDC1s","3","4"),cols = c("indianred","rosybrown","chocolate1","deeppink","deeppink4"))


gene_to_plot <- "Cxcl9"
pdf(paste0("Vln_Plots/Vln_plots_activated_cDC1s_",gene_to_plot,".pdf"),width = 1.2, height = 1.5)
plot<- VlnPlot(immune.combined,"Cxcl9",idents = c("7","9","10","8"),cols = c("darkturquoise","deeppink1","chocolate1","deeppink4")) + NoLegend()
plot + theme(plot.title = element_text(size =(4)), axis.title = element_text(size =(4)),axis.text=element_text(size=6.5))
dev.off()

