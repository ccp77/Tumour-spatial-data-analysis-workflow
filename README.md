# Tumour-spatial-data-analysis-workflow

This repository includes descriptions and scripts for the steps that were used to characterize distributions of immune cell subsets with respect to tumour archetypal regions, based on immunofluorescence data and spatial transcriptomics data. Additionally, it includes a pipeline to identify subclusters of immune cells by scRNAseq that was used to characterize immune cells heterogeneity in tumours. 

Below are the different steps of the workflow, including scripts that were necessary to complete these steps but also details on the parameters that were used with software or GUI interfaces (for example for Imaris and CytoMAP).

1.	Cell segmentation: This was achieved using the Imaris software (https://imaris.oxinst.com/). The details of the parameters that were used are detailed in the script “1_Imaris_segment.R”. Additionally, the script aims at merging the individual csv files for each feature into one metadata file.
2.	Generate archetypal tumour regions using CytoMAP (https://github.com/DrStoltzfus/CytoMAP/tree/main/CytoMAP). This methos was developed to be used through a GUI interface. The details on the different steps are detailed in the script “Two_Regions_CytoMAP.m” 
3.	Generate masks based on cell density to identify regions enriched in one or multiple cell types. This includes the script “three_CellClouds_to_Masks.m”
4.	Characterize cells heterogeneity using scRNAseq: script “4_Seurat_Integration_scRNAseq.R” 
5.	Characterize cells heterogeneity in a spatial transcriptomics dataset (MERFISH): script “5_Merfish_Seurat.R”. Once the cells of interest were identified, masks were generated to identify regions enriched in these cells using script from step 3.

