%% This script opens the MATLAB-based tool CytoMAP that has been developed to be used via a GUI
% Original publication: https://doi.org/10.1016/j.celrep.2020.107523

% launch CytoMAP GUI
CytoMAP

% Steps to be performed with the GUI:
% 1. Load the samples: stucture must be one folder per sample, containing
% one csv file per cell type with the same sets of parameters. See R script
% 1_Imaris_segemnt.R to generate such csv files. 

% 2. Define neighborhoods and choose radius. Typically 50 um.

% 3. Cluster neighborhoods into regions. The number of regions can be
% chosen manually based on the number of cell types included and the level
% of granularity required. Typically, with 4 cell types included, 6-8
% clusters are a good place to start. Different models can be generated and
% compared. Number of clusters can also be defined in an unsupervised way
% using the Davis-Bouldin criterion which calculates the ratio of
% within-cluster to between-cluster distances

% 4. Export the classified neighborhoods and cell types for further
% analyses
