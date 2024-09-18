#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This is a pipeline adapted from https://github.com/MouseLand/cellpose/blob/main/notebooks/run_cyto3.ipynb
# The GUI or code runs the generalist model "cyto3" on images for cell segmentation.
# Inputs are confocal IF images with the channel to be used for cell segmentation 

#%% cyto3 can be used via a GUI

# by using 'conda activate cellpose' on the terminal
# and 'cellpose' to open the GUI. Images can be loaded and cell segmentation executed 
# using 'cyto3' or other trained models.
# ROIs can then be exported for ImageJ


#%% Alternatively, cyto3 can be run on python for larger samples size
# Here, done for one image but can be looped for multiple images.

import numpy as np
import time, os, sys
from urllib.parse import urlparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib as mpl
%matplotlib inline
mpl.rcParams['figure.dpi'] = 200
from cellpose import utils, io

# load image
dataPath = '/Users/piotc/Library/CloudStorage/OneDrive-TheFrancisCrickInstitute/2. Flow Cytometry/CP138_XCR1_WT_KO/WT TCF-1 staining/TCF1 tumour/1 Classify CD8 T cells cellpose/T9B_ROI1'
image_path = os.path.join(dataPath,'CD8_C3-CP138_T9B Block 1-Stitching-01-1.tif')
imgs = mpimg.imread(image_path)

masks_true = dat["masks_true"]

plt.imshow(imgs)
plt.axis('off')  # Hide axes
plt.show()

from cellpose import io

io.logger_setup() # run this to get printing of progress

# model_type='cyto3' or model_type='nuclei'
model = models.Cellpose(gpu=True, model_type="cyto3")

# define CHANNELS to run segementation on
# grayscale=0, R=1, G=2, B=3
# channels = [cytoplasm, nucleus]
# if NUCLEUS channel does not exist, set the second channel to 0
# channels = [0,0]
# IF ALL YOUR IMAGES ARE THE SAME TYPE, you can give a list with 2 elements
# channels = [0,0] # IF YOU HAVE GRAYSCALE
# channels = [2,3] # IF YOU HAVE G=cytoplasm and B=nucleus
# channels = [2,1] # IF YOU HAVE G=cytoplasm and R=nucleus
# OR if you have different types of channels in each image
# channels = [[2,3], [0,0], [0,0]]

# choose the diamter in pixels. Can be evaluated using the cellpose GUI 

masks_pred, flows, styles, diams = model.eval(imgs, diameter=30, channels=[0,0],
                                              niter=2000) 



# Save Cellpose output masks 
io.masks_flows_to_seg(image_path, masks_pred, flows, diams, imgs)

# Create a directory for saving the ROI files
roi_output_dir = 'roi_output'
os.makedirs(roi_output_dir, exist_ok=True)

# Iterate through the masks and export them as ImageJ ROIs
for i in np.unique(masks_pred):
    if i == 0:
        continue  # Skip background (value 0)

    masks_pred = masks_pred == i  # Extract the mask for each cell/ROI
    coords = np.column_stack(np.nonzero(masks_pred))  # Get the coordinates of the mask
    
    # Convert the coordinates to ImageJ ROI
    roi = ImagejRoi.frompoints(coords)
    
    # Save the ROI to a file compatible with Fiji/ImageJ
    roi_filename = os.path.join(roi_output_dir, f'roi_{i}.roi')
    roiwrite(roi_filename, roi)

# Optionally, you can zip all the ROI files to import them all at once in Fiji
import zipfile

with zipfile.ZipFile('rois_for_fiji.zip', 'w') as zipf:
    for root, _, files in os.walk(roi_output_dir):
        for file in files:
            zipf.write(os.path.join(root, file), file)




