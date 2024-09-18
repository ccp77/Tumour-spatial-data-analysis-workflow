# Tumour-spatial-data-analysis-workflow

This repository includes multiple steps that were used to characterize distributions of immune cell subsets with respect to tumour archetypal regions, based on immunofluorescence data and spatial transcriptomics data. Additionally, it includes a pipeline to identify subclusters of immune cells using scRNAseq data to characterize immune cells heterogeneity in tumours. 

Below are the different steps of the workflow, including scripts that were necessary to complete these steps but also including details on the parameters that were used when using user-interfaces (for example for Imaris, CytoMAP).

1.	Cell segmentation: This was achieved using cellpose (https://github.com/MouseLand/cellpose/tree/main) or the Imaris software (https://imaris.oxinst.com/), for which the details of the parameters that were used are detailed. This includes the scripts 1_Imaris_segment.R and 1_segment_Python_cellpose.py
2.	Generate archetypal tumour regions using CytoMAP (https://github.com/DrStoltzfus/CytoMAP/tree/main/CytoMAP). This methos was developed to be used as GUI. The details on the different steps are detailed in the script. This includes the script Two_Regions_CytoMAP.m 
3.	Generate masks based on cell density to identify regions enriched in one or multiple cell types. This includes the script three_CellClouds_to_Masks.m
4.	Characterize the heterogeneity of cells using scRNAseq. This includes the script 4_Seurat_Integration_scRNAseq.R 
5.	Characterize locations of cells within tumours using spatial transcriptomics. This includes the script 5_Merfish_Seurat.R

