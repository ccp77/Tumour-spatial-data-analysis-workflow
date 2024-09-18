# Tumour-spatial-data-analysis-workflow

This repository includes multiple steps that were used to characterize distributions of immune cell subsets with respect to tumour archetypal regions based on immunofluorescence data and spatial transcriptomics data. Additionally, it includes a pipeline to identify subclusters of immune cells using scRNAseq data to characterize immune cells heterogeneity in tumours. 

Below are the different steps of the workflow, either as code or detailing the parameters that were used when using a user-interface.

1.	Cell segmentation: This was achieved using cellpose (https://github.com/MouseLand/cellpose/tree/main) or the Imaris software, for which the details of the parameters that were used are detailed.
2.	Generate archetypal tumour regions using CytoMAP (https://github.com/DrStoltzfus/CytoMAP/tree/main/CytoMAP)
3.	Generate masks based on cell density to identify regions enriched in one or multiple cell types.
4.	Characterize the heterogeneity of cells using scRNAseq
5.	Characterize locations of cells within tumours using spatial transcriptomics

